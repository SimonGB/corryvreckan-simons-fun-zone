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

TH1 *draw_gain_hist(DataFileRoot &A, const char *hname, int mxevts)
{
    LOG(DEBUG) << "##############in drawgainhist";
    double s = A.step()/2.;
    LOG(DEBUG) << s;
    if (mxevts<0)
        mxevts = 999999999;

    int chan0=0;
    int chan1=A.nchan();
    int step=1;
    LOG(DEBUG) << chan1;
    int nbin = (chan1-chan0)/step;
    LOG(DEBUG) << nbin;
    std::string name1 = hname + std::string("_all");
    LOG(DEBUG) << name1;
    TH2 *hst = create_profile2d(name1.c_str(), "Scan", nbin, chan0-step/2.0, chan1-step/2.0, A.npts(), A.from()-s, A.to()-s);
    LOG(DEBUG) << chan0-step/2.0;
    LOG(DEBUG) << chan1-step/2.0;
    LOG(DEBUG) << A.npts();
    LOG(DEBUG) << A.from()-s;
    LOG(DEBUG) << A.to()-s;
    hst->SetXTitle("Channel number");
    hst->SetYTitle(A.type()==1 ? "x10^{3} electrons" : "ns");

    hst->SetZTitle("Signal (ADC units)");
    int ichan, ievt;
    short delay, charge;
    for (ievt=0; A.read_event()==0 && ievt<mxevts; ievt++)
    {
        A.process_event();
        A.get_scan_values(delay, charge);
        // std::cout << "Delay " << delay << " charge " << charge << std::endl;
        //short delay = int(A.value()) >> 16;
        //short charge = int(A.value()) & 0xffff;
        double val;
        // val = factor*charge;
        val = 1024.*charge;
        for (ichan=chan0; ichan<chan1; ichan +=step)
        {
            double ff=1.0;
            double xxx = A.signal(ichan);

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
    LOG(DEBUG) << "saving text file";
    std::ofstream ofile(name);
    LOG(DEBUG) << "ofstream created";
    if (!ofile)
    {
        LOG(DEBUG) << "no outfile";
        std::cout << "Could not open " << name << " for writing" << std::endl;
        return;
    }
    LOG(DEBUG) << "starting loop";
    int nb = h1->GetNbinsX();
    LOG(DEBUG) << "info read from hist";
    for (int ib=1; ib<=nb; ib++)
    {
        LOG(DEBUG) << "reading histogram entry";
        ofile << h1->GetBinCenter(ib) << '\t' << h1->GetBinContent(ib) << std::endl;
    }
    LOG(DEBUG) << "closing file";
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
  const char *ped_f = "/tmp/alibava_ped.ped";
  const char *cal_f = "/tmp/alibava_cal.cal";
  LOG(DEBUG) << "inside alibava loader";

  // Compute the pedestals
  // If not given, use the data file to compute pedestals
  if (!ped_file || !is_text(ped_file))
  //if there is no pedestal file -> TRUE or TRUE -> TRUE
  //if there is an HDF5 or Binary pedestal file -> FALSE or TRUE -> TRUE
  //if there is an ASCII pedestal file -> FALSE or FALSE -> FALSE
  {
      LOG(DEBUG) << "no pedestal file or non-ascii file";
      if (A->valid())
      {
          LOG(DEBUG) << "A is valid";
          A->save();
          LOG(DEBUG) << "A saved";
          A->rewind();
          LOG(DEBUG) << "A rewound";
      }
      else{
          LOG(DEBUG) << "A not valid";
          A->open(data_file);
          LOG(DEBUG) << "Opened data file in A";
      }

      LOG(DEBUG) << "Out of if statement";
      A->compute_pedestals_fast();
      LOG(DEBUG) << "pedestals computed fast";
      A->save_pedestals(ped_f);
      LOG(DEBUG) << "pedestals saved";
      A->restore();
      LOG(DEBUG) << "A restored";
  }
  else{
      LOG(DEBUG) << "Ascii pedestal file!";
      ped_f = ped_file;
  }

  // Get calibration file
  if (cal_file)
  {
      LOG(DEBUG) << "Cal file found";
      if (is_text(cal_file))
      {
          LOG(DEBUG) << "Cal file is text file";
          cal_f = cal_file;
      }
      else
      {
          LOG(DEBUG) << "Cal file is not text file";
          DataFileRoot *B = DataFileRoot::OpenFile(cal_file, ped_f);
          // DataFileRoot *B = DataFileRoot::OpenFile(ped_f,cal_file);
          LOG(DEBUG) << "Datafileroot object for cal created";
          LOG(DEBUG) << B->scan_type();
                                                draw_gain_hist(*B, "hGain");
          LOG(DEBUG) << "gain hist created";
          save_text_file((TH1 *)gDirectory->Get("hGain"), cal_f);
          LOG(DEBUG) << "cal file saved";
          delete B;
          LOG(DEBUG) << "B deleted";
      }
  }
  else
      LOG(DEBUG) << "No cal file";
      cal_f = 0;

  // Analyze the data file
  LOG(DEBUG) << "Analyzing file";
  if (data_file && !A->valid())
      LOG(DEBUG) << "Data file found and valid";
      A->open(data_file);
      LOG(DEBUG) << "Data file opened";

  if (cal_f)
      LOG(DEBUG) << "Loading cal file";
      A->load_gain(cal_f);

  if (ped_f)
      LOG(DEBUG) << "Loading ped file";
      A->load_pedestals(ped_f);

    LOG(DEBUG) << "Loader finished";
}
