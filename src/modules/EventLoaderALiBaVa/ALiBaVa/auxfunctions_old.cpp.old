#include "auxfunctions.h"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <fstream>
#include <TROOT.h>

#include "Data.h"
#include "Tracer.h"
#include "utils.h"
#include "core/utils/log.h"

bool debug_event=false;
void set_debug(bool x)
{
    debug_event = x;
}
bool get_debug()
{
    return debug_event;
}

TH1 *draw_gain_hist(DataFileRoot &A, const char *hname, double factor, bool flip, int chip,
               int odd, int mxevts)
{
    double s = A.step()/2.;
    if (mxevts<0)
        mxevts = 100000000;

    int chan0=0;
    int chan1=A.nchan();
    int step=1;
    if (chip>=0 && chip<2)
    {
        //chip = -1
        // chan0 = 0
        // chan1 = 256
        //chip = 0
        // chan0 = 0
        // chan1 = 128
        //chip = 1
        // chan0 = 128
        // chan1 = 256
        chan0 = chip*128;
        chan1 = (chip+1)*128;
    }
    if (odd>0)
    {
        step = 2;
        if (odd>1)
        {
            chan0 ++;
            chan1 ++;
        }
    }

    int nbin = (chan1-chan0)/step;

    std::string name1 = hname + std::string("_all");
    TH2 *hst = create_profile2d(name1.c_str(), "Scan", nbin, chan0-step/2.0, chan1-step/2.0, A.npts(), A.from()-s, A.to()-s);
    hst->SetXTitle("Channel number");
    hst->SetYTitle(A.type()==1 ? "x10^{3} electrons" : "ns");
    bool delay_scan = A.scan_type() == DataFileRoot::Time;
    if ( delay_scan )
        hst->SetYTitle("Strobe Delay (ns)");

    hst->SetZTitle("Signal (ADC units)");
    int ichan, ievt;
    short delay, charge;
    for (ievt=0; A.read_event()==0 && ievt<mxevts; ievt++)
    {
        A.process_event();
        A.get_scan_values(delay, charge);
        std::cout << "Delay " << delay << " charge " << charge << std::endl;
        //short delay = int(A.value()) >> 16;
        //short charge = int(A.value()) & 0xffff;
        double val;
        if (delay_scan)
        {
            val = delay;
        }
        else
        {
            val = factor*charge;
        }
        for (ichan=chan0; ichan<chan1; ichan +=step)
        {
            double ff=1.0;
            double xxx = A.signal(ichan);
            if (delay_scan)
            {
                if (fabs(xxx)>200.)
                    continue;

                double ss = ichan % 2 ? -1.0 : 1.0 ;
                if (ss*xxx < 0.0 )
                    continue;
            }
            if (flip && xxx<0 )
            {    xxx = -xxx;
            }


            hst->Fill(ichan, val, ff*xxx);
        }
    }

    TH1 *hout=0;
    if (A.type()==1)
    {
        hout = create_h1(hname, "Gain Scan", nbin, chan0-step/2.0, chan1-step/2.0);
        hout->SetXTitle("Channel number");
        hout->SetYTitle("ADC/electron");
        TF1 *gpol = (TF1 *)gROOT->GetFunction("pol1");
        if (A.type()==1)
        {
            for (ichan=chan0; ichan<chan1; ichan +=step)
            {
                TH1 * h1 = (TH1 *)gROOT->FindObject("__h__");
                if (h1)
                    delete h1;

                h1 = hst->ProjectionY("__h__", ichan+1, ichan+1);
                if (h1->GetSumOfWeights()>0)
                {
                    h1->Fit("pol1", debug_event ? "w" : "qw");
                    if (debug_event)
                        gPad->Update();

                    double gain = gpol->GetParameter(1);
                    if (gain<0)
                        gain = -gain;

                    hout->Fill(ichan, gain);
                }
                delete h1;
            }
        }
    }
    return hout;
}


void save_text_file(TH1 *h1, const char *name)
{
    std::ofstream ofile(name);
    if (!ofile)
    {
        std::cout << "Could not open " << name << " for writing" << std::endl;
        return;
    }
    int ib, nb=h1->GetNbinsX();
    for (ib=1; ib<=nb; ib++)
    {
        ofile << h1->GetBinCenter(ib) << '\t' << h1->GetBinContent(ib) << std::endl;
    }
    ofile.close();
}

/*
 * A:          Pointer to a user supplied DataFileRoot (or descendant) object.
 * data_file:  The path of the data file.
 * cal_file:   The path of a calibration file: Either an Alibava data file
 *             produced during a calibration run, or an ASCII text file with as many
 *             lines as channels with gain and offset in each line.
 * ped_file:   The path of a pedestal file. It can be an Alibava data file from where
 *             pedestals can be computed, or an ascii text file with as many lines as
 *             channels, and pedestal and noise for each channel.
 *             If no file is given the data file will be used.
 */
int ALiBaVa_loader(DataFileRoot *A,
                       const char *data_file, const char *cal_file, const char *ped_file)
{

    const char *ped_f = "";
    const char *cal_f = "";

    LOG(DEBUG) << "Computing pedestals";
// ----------------------
    // if (!ped_file || !is_text(ped_file))
    // {
    //         LOG(DEBUG) << "1";
    //         if (A->valid())
    //         {
    //             LOG(DEBUG) << "2";
    //             A->save();
    //             A->rewind();
    //         }
    //         else{
    //             LOG(DEBUG) << "2.5";
    //             A->open(data_file);
    //         }
    //         LOG(DEBUG) << "3";
    //         A->compute_pedestals_fast();
    //         LOG(DEBUG) << "4";
    //         A->save_pedestals(ped_f);
    //         LOG(DEBUG) << "5";
    //         A->restore();
    //         }
    //         else
    //         ped_f = ped_file;
// ----------------
if (!is_text(ped_file))
        LOG(DEBUG) << "Pedestal file is HDF5 format";
        DataFileRoot * PedestalPointer = DataFileRoot::OpenFile(0,ped_file,0);
        LOG(DEBUG) << "1";
        // PedestalPointer->open(ped_file);



        LOG(DEBUG) << "2";
        PedestalPointer->compute_pedestals_fast();
        LOG(DEBUG) << "3";
        PedestalPointer->save_pedestals(ped_f);
        LOG(DEBUG) << "4";
        PedestalPointer->close();
        delete PedestalPointer;
        LOG(DEBUG) << "HDF5 pedestal succesfully converted";
    }
    else
        ped_f = ped_file;
    LOG(DEBUG) << "Pedestals finished";

    // Get the calibration (i.e. if necessary, convert HDF5 to ASCII)
    // If none given, no calibration will be applied
    if (cal_file)
    {
        LOG(DEBUG) << "Calibration file detected";
        if (is_text(cal_file))
        {
            LOG(DEBUG) << "Calibration file is ASCII";
            cal_f = cal_file;
        }
        else
        {
            LOG(DEBUG) << "Calibration file is HDF5";
            DataFileRoot * B = DataFileRoot::OpenFile(cal_file, ped_f, 0);
            draw_gain_hist(*B, "hGain");
            save_text_file((TH1 *)gDirectory->Get("hGain"), cal_f);
            delete B;
            LOG(DEBUG) << "succesfully converted HDF5 calibration file";
        }
    }
    else
        cal_f = 0;
    LOG(DEBUG) << "Calibration loaded succesfully";

    // Load the data, pedestal, and calibration
    A->open(data_file);

    if (cal_f)
        A->load_gain(cal_f);

    A->load_pedestals(ped_f);

    LOG(DEBUG) << "Loader finished";
}
