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
    // LOG(DEBUG) << "TYPE OF A IN DRAW_GAIN_HIST (&) .: " << A.type();
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
    // hst->SetYTitle(A.type()==1 ? "x10^{3} electrons" : "ns");
    hst->SetYTitle("x10^{3} electrons");
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
        // std::cout << "Delay " << delay << " charge " << charge << std::endl;
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
    // if (A.type()==1)
    if (true)
    {
        hout = create_h1(hname, "Gain Scan", nbin, chan0-step/2.0, chan1-step/2.0);
        hout->SetXTitle("Channel number");
        hout->SetYTitle("ADC/electron");
        TF1 *gpol = (TF1 *)gROOT->GetFunction("pol1");
        // if (A.type()==1)
        if (true)
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
        LOG(ERROR) << "Could not open " << name << " for writing" << std::endl;
        return;
    }
    int ib, nb=h1->GetNbinsX();
    for (ib=1; ib<=nb; ib++)
    {
        ofile << h1->GetBinCenter(ib) << '\t' << h1->GetBinContent(ib) << std::endl;
    }
    ofile.close();
}

int ALiBaVa_loader(DataFileRoot *A,
                       const char *data_file, const char *ped_file, const char *cal_file,
                       TH1F* & pedestalValues, TH1F* & noiseValues, TH1F* & correctedPedestalValues, TH1F* & correctedNoiseValues,
                       const int lowerChannel, const int upperChannel)
{
    const char *ped_f = "alibava_ped.ped";
    const char *cal_f = "alibava_cal.cal";

    // Create a pointer with the pedestal file
    DataFileRoot * PedestalPointer = DataFileRoot::OpenFile(ped_file);

    // Calculate the pedestals, and compute and apply the common mode noise correction
    PedestalPointer->compute_pedestals_alternative(lowerChannel, upperChannel);
    pedestalValues = PedestalPointer->show_pedestals(lowerChannel, upperChannel, false);
    noiseValues = PedestalPointer->show_noise(lowerChannel, upperChannel, false);
    PedestalPointer->compute_cmmd_alternative(lowerChannel, upperChannel);
    correctedPedestalValues = PedestalPointer->show_pedestals(lowerChannel, upperChannel, true);
    correctedNoiseValues = PedestalPointer->show_noise(lowerChannel, upperChannel, true);

    // Save the calculated pedestal information in a temporary file
    PedestalPointer->save_pedestals(ped_f);
    PedestalPointer->close();
    delete PedestalPointer;







    // INSERT CALIBRATION
    // Still in testing phase, non functional!
    DataFileRoot * CalibrationPointer = DataFileRoot::OpenFile(cal_file);
    CalibrationPointer->load_pedestals(ped_f, kTRUE);
    double channelCalVal[upperChannel-lowerChannel+1];
    double channelCalValNeg[upperChannel-lowerChannel+1];
    int posChargeEventCount[upperChannel-lowerChannel+1];
    int negChargeEventCount[upperChannel-lowerChannel+1];
    int nCalEvents=0;
    int DebugCOUNT = 0;
    int DebugCHARGE = -1;

    for(int chan = 0; chan <= upperChannel-lowerChannel; chan++){
        channelCalVal[chan] = 0.;
        channelCalValNeg[chan] = 0.;
        posChargeEventCount[chan] = 0;
        negChargeEventCount[chan] = 0;
    }

    while(CalibrationPointer->read_event()==1)
    {
        short delay, charge; // delay is pointless in this context
        CalibrationPointer->process_event();
        CalibrationPointer->get_scan_values(delay,charge);

        if(DebugCHARGE != (int)charge){
          DebugCHARGE = (int)charge;
          std::cout << DebugCHARGE << "\n";
        //   DebugCOUNT++;
        }

        if(charge <= 5) continue;
        double actualCharge = (double)charge * 1025 * 2;
        nCalEvents++;
        // if(nCalEvents == 10) break;

        // std::cout << "Event " << nCalEvents-1 << " with charge " << charge << "\n";
        for(int chan = 0; chan <= upperChannel-lowerChannel; chan++){
          double CalADC = CalibrationPointer->ADC_signal(chan);

          if(CalADC>0){
            // if(charge == 10 && chan == 0) DebugCOUNT++;
            channelCalVal[chan] = channelCalVal[chan] + actualCharge/CalADC;
            // std::cout << "      " << " Pos Channel " << chan << " ADC " << CalADC << "  summed_constant = " << channelCalVal[chan] << "\n";
            posChargeEventCount[chan]++;
            // std::cout << "      " << "             " << "count = " << posChargeEventCount[chan] << "\n";
          }
          else if(CalADC<0){
            channelCalValNeg[chan] = channelCalValNeg[chan] - actualCharge/CalADC;
            // std::cout << "      " << " Neg Channel " << chan << " ADC " << CalADC << "  summed_constant = " << channelCalValNeg[chan] << "\n";
            negChargeEventCount[chan]++;
            // std::cout << "      " << "             " << "  count = " << negChargeEventCount[chan] << "\n";

          }
          // std::cout << "          Channel " << chan+lowerChannel << " ADC = " << CalADC << "\n";
        }
    }

    double average = 0;
    for(int chan = 0; chan <= upperChannel-lowerChannel; chan++){
      channelCalVal[chan] = channelCalVal[chan]/posChargeEventCount[chan];
      channelCalValNeg[chan] = channelCalValNeg[chan]/negChargeEventCount[chan];
      std::cout << "Channel " << lowerChannel+chan << " cal-val= " << (channelCalVal[chan]+channelCalValNeg[chan])/2 << "\n";
      average+=(channelCalVal[chan]+channelCalValNeg[chan])/2;
    }
    average = average/(upperChannel-lowerChannel+1);
    std::cout << "Average " << average << "\n" << DebugCOUNT<< "\n";
      //   std::cout << "Event " << nCalEvents-1 << " with charge " << charge << "\n";
        // for(int chan = 0; chan <= upperChannel-lowerChannel; chan++){
        //   std::cout << "Channel " << chan+lowerChannel << " ADC = " << CalibrationPointer->ADC_signal(chan) << "\n";
        // }









    // Load the calculated pedestal info into the original datafile
    A->load_pedestals(ped_f, kTRUE);
    //A->load_gain(cal_f);

    return -1;
}

void crosstalk_correction(){

}
