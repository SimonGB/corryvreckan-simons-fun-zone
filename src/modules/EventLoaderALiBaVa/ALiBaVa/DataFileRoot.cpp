#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <csignal>
#include <cstring>
#include <numeric>
#include <sys/stat.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TProfile2D.h>

#ifdef HAVE_CONFIG_H
#ifdef G_OS_WIN32
#include "config_win.h"
#else
#include <config.h>
#endif
#endif

#include "DataFileRoot.h"
#include "AsciiRoot.h"
#include "HDFRoot.h"
#include "utils.h"

// #ifndef HAVE_HDF5
// #define HAVE_HDF5
// #endif

#ifdef __APPLE__
#define sighandler_t sig_t
#endif
#ifdef __WIN32__
#define sighandler_t __p_sig_fn_t
#endif
#ifdef __MINGW32__
#define sighandler_t __p_sig_fn_t
#endif
#if _WINDOWS
typedef void (*sighandler_t)(int);
#endif

bool _A_do_run = true;
void _A_got_intr(int)
{
    _A_do_run = false;
}

DataFileRoot::DataFileRoot(const char *nam, const char *pedfile, const char *gainfile)
    : _nchan(max_nchan),_seedcut(5.), _neighcut(3.), _average_gain(1.),
      _version(2), _polarity(1), _t1(0.0), _t2(99999.)
{
    int i;

    for (i=0;i<max_nchan;i++)
    {
        _ped[i] = 0.;
        _gain[i] = 1.;
        _noise[i] = 1.;
        _data.data[i] = 0;
        _mask[i] = false;
    }

////  Cannot call open in the constructor since overrides are not
////  available in the constructor
//    if (nam)
//        open(nam);

    if (pedfile)
        load_pedestals(pedfile);

    if (gainfile)
        load_gain(gainfile);

    // TODO: this loads a "fixed" name file. Be more general...
    // load_masking();
}

DataFileRoot::~DataFileRoot()
{
    // TODO Auto-generated destructor stub
}

void DataFileRoot::reset_data()
{
    memset(&_data, 0, sizeof(_data));
}


void DataFileRoot::set_timecut(double t1, double t2)
{
    if (t1>0)
        _t1 = t1;

    if (t2>0)
        _t2 = t2;

    if ( _t1 > _t2 )
    {
        t1 = _t1;
        _t1 = _t2;
        _t2 = t1;
    }
}

bool DataFileRoot::valid_time(double tim) const
{
    return (tim>=_t1 && tim<=_t2);
}

void DataFileRoot::set_data(int nchan, const unsigned short int *data)
{
    int i;
    _nchan = nchan;
    for (i=0;i<_nchan;i++)
        _data.data[i] = data[i];
}


TH1 *DataFileRoot::show_pedestals()
{
    int ic;
    TH1 *hst = create_h1("hPed","Pedestals",nchan(),-0.5, nchan()-0.5);
    hst->SetYTitle("ADCs");
    hst->SetXTitle("Channel no.");
    for (ic=0; ic<nchan(); ic++)
        hst->SetBinContent(ic+1, _ped[ic]);

    return hst;
}

TH1 *DataFileRoot::show_noise()
{
    int ic;
    TH1 *hst = create_h1("hNoise","Noise",nchan(),-0.5, nchan()-0.5);
    if (gain()==1)
    {
        hst->SetYTitle("ADCs");
    }
    else
    {
        hst->SetYTitle("e^{-} ENC");
    }
    hst->SetXTitle("Channel no.");
    for (ic=0; ic<nchan(); ic++)
        hst->SetBinContent(ic+1, noise(ic));

    return hst;
}

void DataFileRoot::compute_pedestals_fast(int mxevts, double wped, double wnoise)
{
    int i, ievt;

    if (!valid())
        return;

    if (mxevts<0)
        mxevts = 100000000;

    for (i=0;i<max_nchan;i++)
        _ped[i] = _noise[i] = 0.;

    // std::cout << "Computing fast pedestals..." << std::endl;
    for (ievt=0; read_data()==0 && ievt<mxevts; ievt++)
    {
        if (!(ievt%1000))
        {
            // std::cout << "event " << ievt << ": " << _ped[i] << "\n";
            // std::cout << "\revent " << std::setw(10) << ievt << std::flush;
        }
        common_mode();
        for (i=0; i<nchan(); i++)
        {
            // TODO: figure out how to determine the chip number when
            //       Plugin::filter_event has been called
            int ichip = i/128;
            // IF noise is 0, set it arbitrarily to 1.
            if (_noise[i]==0.)
                _noise[i] = 1.;

            if (_ped[i]==0.)
            {

                // If pedestal is not yet computed we assume the current
                // channel value should not be too far
                _ped[i] = _data.data[i];
            }
            else
            {
                // Do the pedestal and noise correction
                double corr;
                double xs;

                _signal[i] = _data.data[i] - _ped[i];
                corr = _signal[i] * wped;

                xs = (_signal[i]-_cmmd[ichip])/_noise[i];
                if (corr > 1.)
                    corr = 1.;

                if (corr < -1)
                    corr = -1.;

                _ped[i] += corr;

                if (fabs(xs) < 3.)
                {
                    _noise[i] = _noise[i]*(1.0-wnoise) + xs*xs*wnoise;
                }
            }
        }
    }
    // std::cout << "\nDone" << std::endl;
    rewind();
}

void DataFileRoot::compute_pedestals_alternative()
{
    int mxevts = 100000000;
    int max_nchan = 128;

    int i, ievt;
    std::vector <double> pedestal_data[max_nchan];
    double pedestal_average[max_nchan];
    double pedestal_stdev[max_nchan];

    if (!valid())
        return;

    for (i=0;i<max_nchan;i++)
        _ped[i] = _noise[i] = 0.;

    for (ievt=0; read_data()==0 && ievt<mxevts; ievt++)
    {
        for (i=0; i<max_nchan; i++)
        {
            pedestal_data[i].push_back(_data.data[i]);
        }
    }
    for (i=0;i<max_nchan;i++){

      pedestal_average[i] = std::accumulate(pedestal_data[i].begin(), pedestal_data[i].end(), 0.0);
      pedestal_average[i] = pedestal_average[i] / pedestal_data[i].size();

      pedestal_stdev[i] = std::inner_product(pedestal_data[i].begin(), pedestal_data[i].end(), pedestal_data[i].begin(), 0.0);
      pedestal_stdev[i] = std::sqrt(pedestal_stdev[i] / pedestal_data[i].size() - pedestal_average[i] * pedestal_average[i]);

      _ped[i] = pedestal_average[i];
      _noise[i] = pedestal_stdev[i];

    }
    rewind();
}

void DataFileRoot::compute_cmmd_alternative()
{
    int mxevts = 100000000;
    int max_nchan = 128;
    int nEvents;

    int i, ievt;

    double event_bias = 0; //common mode noise per pedestal events
    double cmn = 0; //common mode noise averaged over all pedestal events

    std::vector <double> pedestal_data[max_nchan];
    std::vector <double> corrected_pedestal_data[max_nchan];
    double pedestal_average[max_nchan];
    double pedestal_stdev[max_nchan];

    if (!valid())
        return;

    for (ievt=0; read_data()==0 && ievt<mxevts; ievt++)
    {
        for (i=0; i<max_nchan; i++)
        {
            pedestal_data[i].push_back(_data.data[i]);
        }
    }

    nEvents = pedestal_data[1].size();

    for (ievt=0; ievt<nEvents; ievt++){
        double event_sum = 0;
        for (i=0; i<max_nchan; i++)
        {
            event_sum += (pedestal_data[i].at(ievt) - _ped[i]);
        }
        event_bias = event_sum/max_nchan;
        cmn += event_bias;
        for (i=0; i<max_nchan; i++)
        {
            corrected_pedestal_data[i].push_back(pedestal_data[i].at(ievt) - event_bias);
        }
    }

    _cmmd[0] = cmn/nEvents;

    for (i=0;i<max_nchan;i++){

        pedestal_average[i] = std::accumulate(corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), 0.0);
        pedestal_average[i] = pedestal_average[i] / corrected_pedestal_data[i].size();

        pedestal_stdev[i] = std::inner_product(corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), corrected_pedestal_data[i].begin(), 0.0);
        pedestal_stdev[i] = std::sqrt(pedestal_stdev[i] / corrected_pedestal_data[i].size() - pedestal_average[i] * pedestal_average[i]);

        _ped[i] = pedestal_average[i];
        _noise[i] = pedestal_stdev[i];
    }

    rewind();
}

TH2 *DataFileRoot::compute_pedestals(int mxevts, bool do_cmmd)
{
    if (!valid())
        return 0;

    if (mxevts<0)
        mxevts = 100000000;

    int ievt, ichan;
    TH2 *hst = create_h2("hRaw","Raw data",nchan(), -0.5,nchan()-0.5, 256, -0.5,1023.5);
    TH2 *hsts = create_h2("hSig","Signal",nchan(), -0.5,nchan()-0.5,256, -127.5,127.5);


    std::cout << "Computing pedestals..." << std::endl;
    for (ievt=0; read_data()==0 && ievt<mxevts; ievt++)
    {
        process_event(do_cmmd);
        for (ichan=0; ichan<nchan(); ichan++)
            // TODO: get right chip number in all situations (after calling set_data)
            hst->Fill(ichan, data(ichan)-get_cmmd(ichan/128));

        if (!(ievt%1000))
        {
            std::cout << "\revent " << std::setw(10) << ievt << std::flush;
        }
    }
    std::cout << "\nDone" << std::endl;
    rewind();

    // TODO: _nchan can be updated in an event by event basis
    //       while here we are assuming that it is the same
    //       for all the events
    for (ichan=0; ichan<nchan(); ichan++)
    {
        TF1 *g = new TF1("g1", "gaus");
        TH1 *h1 = hst->ProjectionY("__hx__", ichan+1, ichan+1);
        g->SetParameters(h1->GetSumOfWeights(), h1->GetMean(), h1->GetRMS());
        g->SetRange(h1->GetMean()-2.5*h1->GetRMS(), h1->GetMean()+2.5*h1->GetRMS());
        h1->Fit("g1", "q0wr");
        _ped[ichan] = h1->GetFunction("g1")->GetParameter(1);
        _noise[ichan] = h1->GetFunction("g1")->GetParameter(2);
        delete h1;
        delete g;
    }

    rewind();
    for (ievt=0; read_data()==0 && ievt<mxevts; ievt++)
    {
        process_event(do_cmmd);
        for (ichan=0; ichan<nchan(); ichan++)
            hsts->Fill(ichan, signal(ichan));

        if (!(ievt%1000))
        {
            std::cout << "\revent " << std::setw(10) << ievt << std::flush;
        }
    }
    std::cout << "\nDone" << std::endl;
    rewind();

    return hst;
}




void DataFileRoot::find_clusters(int ichip)
{
    int chan0=0;
    int chan1=255;
    if (ichip>=0 && ichip<2)
        {
            chan0 = ichip*128;
            chan1 = (ichip+1)*128 -1;
        }

    std::ostringstream ostr;
    ostr << chan0 << '-' << chan1;
    ChanList C(ostr.str().c_str());

    clear();
    find_clusters(C);
    _hits = C.hit_list();
}

void DataFileRoot::find_clusters(ChanList &C)
{
    // TODO: figure out how to determine the chip number in
    //       all the situations
    int i, j, imax=-1, left, right;
    double mxsig=-1.e20, sg, val;
    bool *used = new bool[C.Nch()];

    for (i=0;i<C.Nch();i++)
    {
        used[i]= _mask[C[i]] ? true : false;
    }



    while (true)
    {
        /*
         * Find the highest
         */
        imax = -1;
        for (j=0; j<C.Nch(); j++)
        {
            i = C[j];
            if (used[j] || _signal[i]*polarity()<0.)
                continue;

            if ( polarity()*sn(i) > _seedcut)
            {
                val = fabs(signal(i));
                if (mxsig<val)
                {
                    mxsig = val;
                    imax = j;
                }
            }
        }

        if (imax<0 || imax >= C.Nch() )
            break;

        sg = signal(C[imax]);
        used[imax]=true;
        // Now look at neighbors
        // first to the left
        left = C[imax];
        for (j=imax-1;j>=0;j--)
        {
            i = C[j];
            if ( used[j] || _signal[i]*polarity()<0.)
                break;

            if ( fabs(sn(i)) > _neighcut )
            {
                used[j] = true;
                sg += signal(i);
                left = i;
            }
        }

        // now to the right
        right = C[imax];
        for (j=imax+1;j<C.Nch();j++)
        {
            i = C[j];
            if ( used[j] || _signal[i]*polarity()<0.)
                break;
            if ( fabs(sn(i))>_neighcut )
            {
                used[j] = true;
                sg += signal(i);
                right = i;
            }
        }
        C.add_hit(Hit(C[imax], left, right, sg));
    }
    delete [] used;
}

void DataFileRoot::save_pedestals(const char *fnam)
{
    std::ofstream ofile(fnam);
    if (!ofile)
    {
        std::cout << "Could not open " << fnam << " to save pedestals." << std::endl;
        return;
    }

    // TODO: _nchan can be updated in an event by event basis
    //       while here we are assuming that it is the same
    //       for all the events
    ofile << _cmmd[0] << "\n";
    int i;
    for (i=0; i<nchan(); i++)
    {
        ofile << _ped[i] << "\t" << _noise[i] << "\n";
    }
    ofile.close();
}

void DataFileRoot::load_pedestals(const char *fnam, bool show)
{
    std::ifstream ifile(fnam);
    if (!ifile)
    {
        std::cout << "Could not open " << fnam << " to load pedestals." << std::endl;
        return;
    }

    ifile >> _cmmd[0] >> std::ws;

    int i;
    for (i=0; i<max_nchan; i++)
    {
        if (ifile.eof())
            break;

        ifile >> _ped[i] >> std::ws >> _noise[i] >> std::ws;
        // _mask[i] = (_noise[i]>20. || _noise[i]<=0.);
    }
    ifile.close();
    if (show)
    {
        TCanvas *pedcanvas = create_canvas("Pedcanvas", "Pedestal Values", 600, 400);
        TH1 *pedestalhisto = create_h1("pedestalhisto", "Pedestal Values", 128, -0.5, 127.5);
        for (i=0; i<128; i++)
        {
            pedestalhisto->Fill(i, _ped[i]);
        }
        pedcanvas->cd(1);
        pedestalhisto->Draw();

        // for (i=0; i<256; i++)
        // {
          // std::cout << "channel " << i << ": " << _ped[i] << "\n";
        // }
    }
}

void DataFileRoot::load_masking(const char *fnam)
{
    std::ifstream ifile(fnam);
    if (!ifile)
    {
        std::cout << "Could not open masked.txt. " << std::endl;
        return;
    }
    int val;
    for (int i=0; i<500; i++)
    {
        ifile >> val >> std::ws;
        if (ifile.eof())
            break;
        if (val>255)
        {
            std::cout << "A value is greater than 255, causing an overflow crash. Please check the text file again. It has been set to 1 for continuation purposes. " << std::endl;
            val = 1;
        }
        _mask[val] = true;
    }
}

void DataFileRoot::load_gain(const char *fnam)
{
    std::ifstream ifile(fnam);
    if (!ifile)
    {
        std::cout << "Could not open " << fnam << " to load the gain." << std::endl;
        return;
    }
    int i;
    int ichan;
    double val, xn, xm;
    xn=xm=0.;
    for (i=0; i<max_nchan; i++)
    {
        ifile >> ichan >> std::ws;
        if (ifile.eof())
            break;
        ifile >> val;
        if (ifile.eof())
            break;

        xn++;

        xm += val;
        _gain[ichan] = val;

        ifile >> std::ws;
        if (ifile.eof())
            break;
    }
    if (xn>0)
    {
        _average_gain = xm/xn;
    }
    ifile.close();
}

// This function processes the event, it subtracts pedestals from data
// and fills in the signal/noise ratio

void DataFileRoot::process_event()
{
    int i;
    for (i=0; i<nchan(); i++){
      _signal[i] = _data.data[i]-_ped[i];
      _sn[i] = _signal[i]/_noise[i];
    }
}

void DataFileRoot::add_channel_list(const ChanList &C)
{
    chan_list.push_back(C);
}


void DataFileRoot::common_mode()
{
    ChanList C("0-127");
    // ChanList C("18-122");
    common_mode(C);

    _cmmd[0] = C.CommonMode();
    _cnoise[0] = C.Noise();


    C.Set("128-255");
    common_mode(C);

    _cmmd[1] = C.CommonMode();
    _cnoise[1] = C.Noise();
}

void DataFileRoot::common_mode(ChanList &C, bool correct)
{
    int ip, i, j;

    double mean, sm, xn, xx, xs, xm, tmp;
    bool use_it;
    mean = sm = 0.;

    for (ip=0;ip<3;ip++)
    {
        xn = xs = xm = 0.;
        for (j=0; j<C.Nch(); j++)
        {
            i = C[j];
            if (_mask[i])
                continue;

            use_it = true;
            xx = data(i) - _ped[i];
            if (ip)
            {
                tmp = fabs((xx-mean)/sm);
                use_it = (tmp<2.5);
            }
            if (use_it)
            {
                xn++;
                xm += xx;
                xs += xx * xx;
            }
        }
        if (xn>0.)
        {
            mean = xm / xn;
            sm = sqrt( xs/xn - mean*mean);
        }
        //  std::cout << "...iter " << ip << ": xm " << mean << " xs: " << sm << std::endl;
    }
    C.CommonMode(mean);
    C.Noise(sm);

    if (correct)
    {
        for ( j=0; j<C.Nch(); j++ )
        {
            i = C[j];
            _signal[i] = _data.data[i]-_ped[i] - C.CommonMode();
            _sn[i] = (_noise[i] >1. && !_mask[i] ? _signal[i]/_noise[i] : 0.);
        }
    }
}




void DataFileRoot::spy_data(bool with_signal, double t0, double t1, int nevt)
{
    TVirtualPad *pad;
    if (!valid())
        return;

    sighandler_t old_handler = ::signal(SIGINT, _A_got_intr);
    _A_do_run = true;

    TCanvas *cnvs = (TCanvas *)gROOT->FindObject("cnvs");
    if (cnvs)
    {
        cnvs->Clear();
    }
    else
       cnvs = new TCanvas("cnvs","cnvs", 700, 800);

    cnvs->Divide(2,3);


    TH1 *hsignal = create_h1("hsignal","signal (ADC)",256, -0.5, 255.0);
    hsignal->SetXTitle("Channel");
    hsignal->SetYTitle("ADC");
    hsignal->SetMinimum(-300);
    hsignal->SetMaximum(300);

    TH1 *helec = create_h1("helec","signal (elec)", 256, -0.5, 255.5);
    helec->SetXTitle("Channel");
    helec->SetYTitle("electrons");
    helec->SetMinimum(-300/gain());
    helec->SetMaximum(300/gain());

    TH1 *hraw = create_h1("hraw","Raw Data (around 512.)",256, 0., 256.);
    hraw->SetXTitle("Channel");
    hraw->SetYTitle("ADC");
    hraw->SetMinimum(-300);
    hraw->SetMaximum(+300);

    TH1 *hrawc = create_h1("hrawc","Raw Data (no commd)",256, 0., 256.);
    hrawc->SetXTitle("Channel");
    hrawc->SetYTitle("ADC");
    hrawc->SetMinimum(-300);
    hrawc->SetMaximum(+300);


    TH1 *hcmmd[2];
    hcmmd[0] = create_h1("hcmmd0","Common mode (Chip 0)",50,-100.,100.);
    hcmmd[0]->SetXTitle("Common mode");
    hcmmd[1] = create_h1("hcmmd1","Common mode (Chip 1)",50,-100.,100.);
    hcmmd[1]->SetXTitle("Common mode");

    int ievt,jevt;
    for (ievt=jevt=0; read_event()==0 && _A_do_run && ievt<nevt;jevt++)
    {
        std::vector<ChanList>::iterator icl;
        process_event();
        if ( ! (t0==t1&& t0==0 ) )
        {
            if (t0>t1)
            {
                double tmp = t0;
                t0 = t1;
                t1 = tmp;
            }
            double tval = time();
            if (tval<t0 || tval>t1)
                continue;
        }
        if ( chan_list.empty() )
            find_clusters();
        else
        {
            for (icl=chan_list.begin(); icl!=chan_list.end(); ++icl)
            {
                icl->clear_hits();
                find_clusters(*icl);
                for (int ii=0; ii<icl->nhits(); ii++)
                    _hits.push_back( icl->get_hit(ii) );
            }
        }
        if ( with_signal && empty())
            continue;

        int i,ichip=-1;
        for (i=0; i<nchan(); i++)
        {
            // TODO: figure out chip number
            if (!(i%128))
                ichip++;

            hsignal->SetBinContent(i+1, _signal[i]);
            helec->SetBinContent(i+1, signal(i));
            hraw->SetBinContent(i+1,data(i)-512.);
            hrawc->SetBinContent(i+1, data(i)-_ped[i]);
            // TODO: why we draw the signal + common mode ?
            //       May be cause signal should be ~0...
            hcmmd[ichip]->Fill(_signal[i]+get_cmmd(ichip));
        }
        pad = cnvs->cd(1);
        pad->SetGrid(1,1);
        hsignal->Draw();
        pad = cnvs->cd(2);
        pad->SetGrid(1,1);
        helec->Draw();

        pad = cnvs->cd(3);
        pad->SetGrid(1,1);
        hraw->Draw();

        pad = cnvs->cd(4);
        pad->SetGrid(1,1);
        hrawc->Draw();

        pad = cnvs->cd(5);
        pad->SetGrid(1,1);
        hcmmd[0]->Draw();

        pad = cnvs->cd(6);
        pad->SetGrid(1,1);
        hcmmd[1]->Draw();

        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << "*** Event " << jevt << " *****" << std::endl;
        std::cout << "Common Mode:" << std::endl
                  << "   Chip 0 " << std::setw(6) << std::setprecision(1) << get_cmmd(0) << " noise: " << get_cnoise(0)
                  << std::endl
                  << "   Chip 1 " << std::setw(6) << std::setprecision(1) << get_cmmd(1) << " noise: " << get_cnoise(1)
                  << std::endl;

        std::cout << "Time: " << time() << " ns" << std::endl;
        std::cout << "Clusters: " << std::endl;

        HitList::iterator ip;
        for (ip=begin(); ip!=end(); ++ip)
        {
            std::cout << "   Seed: " << ip->center()
                      << " sig: "
                      << std::setw(6) << std::setprecision(1) << ip->signal()
                      << " left: " << ip->left() << " right: " << ip->right()
                      << std::endl;
            std::cout << '\t' << "channels: " << std::endl;
            int j;
            for (j=ip->left();j<=ip->right();j++)
                std::cout << "\t   " << j << " sn: " << _sn[j] << " signal: " << _signal[j] << " noise: " << _noise[j] << '\n';
            std::cout << std::endl;
        }

        cnvs->Update();
        ievt++;
    }
    std::cout << std::endl;
    _A_do_run= true;
    ::signal(SIGINT, old_handler);
}



bool is_text(const char *fnam)
{
    int nc;
    char buffer[1024];
    std::ifstream ifile(fnam);
    if (!fnam)
        return false;

    ifile.read(buffer, sizeof(buffer));
    nc = ifile.gcount();
    ifile.close();
    if (!nc) // empty files are text
    {
        return true;
    }

    std::string ss(buffer, nc);
    ifile.close();

    if ( ss.find('\0') != ss.npos )
        return false;

    double nontext = 0.;
    double ntotal = 0.;
    std::string::iterator ip;
    for (ip=ss.begin(); ip!=ss.end(); ++ip)
    {
        ntotal++;
        char c = *ip;
        if ( (c<' ' || c >'~') && !strchr("\n\t\r\b", c) )
            nontext++;
    }
    if ( nontext/ntotal > 0.3 )
        return false;

    return true;
}

unsigned int DataFileRoot::clock_counter() const
{
    return _data.clock;
}

double DataFileRoot::time() const
{
    return 0.0;
}
double DataFileRoot::temp() const
{
    return 0.0;
}


DataFileRoot *DataFileRoot::OpenFile(const char *nam, const char *pedfile, const char *gainfile)
{
	struct stat sb;
	if ( stat(nam, &sb) == -1 )
		return 0;

#ifdef HAVE_HDF5
    std::ifstream ifile(nam);
    if (!ifile)
        return 0;

    char buf[5]={'\0'};
    ifile.read(buf, 4);
    ifile.close();
    std::string idf(buf+1);

    if (idf=="HDF")
        return new HDFRoot(nam, pedfile, gainfile);
    else
        return new AsciiRoot(nam, pedfile, gainfile);
#else
    return new AsciiRoot(nam, pedfile, gainfile);
#endif
}
