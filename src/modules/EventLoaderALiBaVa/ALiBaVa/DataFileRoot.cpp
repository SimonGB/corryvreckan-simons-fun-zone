#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <cmath>
#include <csignal>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <sys/stat.h>

#ifdef HAVE_CONFIG_H
#ifdef G_OS_WIN32
#include "config_win.h"
#else
#include <config.h>
#endif
#endif

#include "AsciiRoot.h"
#include "DataFileRoot.h"
#include "HDFRoot.h"
//#include "utils.h"

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
void _A_got_intr(int) {
    _A_do_run = false;
}

DataFileRoot::DataFileRoot(const char* nam, const char* pedfile, const char* gainfile)
    : _nchan(max_nchan), _seedcut(5.), _neighcut(3.), _average_gain(1.), _version(2), _polarity(1), _t1(0.0), _t2(99999.),
      _roi({}) {
    int i;

    for(i = 0; i < max_nchan; i++) {
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

    if(pedfile)
        load_pedestals(pedfile);

    // if (gainfile)
    // load_gain(gainfile);

    // TODO: this loads a "fixed" name file. Be more general...
    // load_masking();
}

DataFileRoot::~DataFileRoot() {
    // TODO Auto-generated destructor stub
}

void DataFileRoot::set_data(int nchan, const unsigned short int* data) {
    int i;
    _nchan = nchan;
    for(i = 0; i < _nchan; i++)
        _data.data[i] = data[i];
}

void DataFileRoot::reset_data() {
    memset(&_data, 0, sizeof(_data));
}

void DataFileRoot::set_ROI(std::vector<unsigned int> bounds) {
    int n_arrays;
    int i, j;
    std::vector<unsigned int> temp_array;

    n_arrays = bounds.size() / 2;

    for(i; i <= n_arrays - 1; i++) {
        temp_array = {};
        for(j = bounds[2 * i]; j <= bounds[2 * i + 1]; j++) {
            temp_array.push_back(j);
        }
        _roi.insert(_roi.end(), temp_array.begin(), temp_array.end());
    }
    /*
    std::cout << "In set_ROI" << std::endl;
    std::cout << _roi.size() << std::endl;
    for (auto k: _roi)
        std::cout << k << std::endl;
    std::cout << _roi.size() << std::endl;
    std::cout << "end set_ROI" << std::endl;
     */
}

void DataFileRoot::set_timecut(double t1, double t2) {
    if(t1 > 0)
        _t1 = t1;

    if(t2 > 0)
        _t2 = t2;

    if(_t1 > _t2) {
        t1 = _t1;
        _t1 = _t2;
        _t2 = t1;
    }
}

bool DataFileRoot::valid_time(double tim) const {
    return (tim >= _t1 && tim <= _t2);
}

void DataFileRoot::compute_pedestals_alternative() {
    int mxevts = 10000000;
    // int max_nchan = 128;

    int i, ievt;
    std::vector<double> pedestal_data[max_nchan];
    double pedestal_average[max_nchan];
    double pedestal_stdev[max_nchan];

    if(!valid())
        return;

    for(unsigned int i : _roi)
        _ped[i] = _noise[i] = 0.;

    for(ievt = 0; read_data() == 1 && ievt < mxevts; ievt++) {
        // std::cout << "IT DOESN'T GET TO THIS LOOP?";
        // now it does, read_data() return value needs to be checked, different for HDF5 and binary
        for(unsigned int i : _roi) {
            pedestal_data[i].push_back(_data.data[i]);
            // std::cout << "Pedestal data: " << _data.data[i] << "\n";
        }
    }

    for(unsigned int i : _roi) {
        // for(int test=0;test<101;test++){ std::cout << pedestal_data[i][test] << "\n";}
        pedestal_average[i] = std::accumulate(pedestal_data[i].begin(), pedestal_data[i].end(), 0.0);
        // std::cout << pedestal_average[i] << "\n";

        pedestal_average[i] = pedestal_average[i] / pedestal_data[i].size();
        // std::cout << pedestal_average[i] << "\n";

        pedestal_stdev[i] =
            std::inner_product(pedestal_data[i].begin(), pedestal_data[i].end(), pedestal_data[i].begin(), 0.0);
        pedestal_stdev[i] =
            std::sqrt(pedestal_stdev[i] / pedestal_data[i].size() - pedestal_average[i] * pedestal_average[i]);

        _ped[i] = pedestal_average[i];
        _noise[i] = pedestal_stdev[i];
        // std::cout << "Before cmmd" << _ped[i] << "\n";
    }

    rewind();
}

void DataFileRoot::compute_cmmd_alternative() {
    int mxevts = 10000000;
    int max_nchan = _nchan; // was 128 originally
    int nEvents = 0;

    int i, ievt;

    double event_bias = 0; // common mode noise per pedestal events
    double cmn = 0;        // common mode noise averaged over all pedestal events

    std::vector<double> corrected_pedestal_data[max_nchan];
    double pedestal_average[max_nchan];
    double pedestal_stdev[max_nchan];

    if(!valid())
        return;

    for(ievt = 0; read_data() == 1 && ievt < mxevts; ievt++) {
        nEvents++;
        double event_sum = 0;
        for(unsigned int i : _roi) {
            // std::cout << _data.data[i] << "\n";
            event_sum += (_data.data[i] - _ped[i]);
        }
        event_bias = event_sum / (_roi.size());
        cmn += event_bias;
        // std::cout << cmn << "\n";
        for(unsigned int i : _roi) {
            corrected_pedestal_data[i].push_back(_data.data[i] - event_bias);
        }
    }

    _cmmd[0] = cmn / nEvents;

    for(unsigned int i : _roi) {

        pedestal_average[i] = std::accumulate(corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), 0.0);
        pedestal_average[i] = pedestal_average[i] / corrected_pedestal_data[i].size();

        pedestal_stdev[i] = std::inner_product(
            corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), corrected_pedestal_data[i].begin(), 0.0);
        pedestal_stdev[i] =
            std::sqrt(pedestal_stdev[i] / corrected_pedestal_data[i].size() - pedestal_average[i] * pedestal_average[i]);

        _ped[i] = pedestal_average[i];
        _noise[i] = pedestal_stdev[i];

        // std::cout << "After cmmd" << _ped[i] << "\n";
        // std::cout << _noise[i] << "\n";
    }

    rewind();
}

void DataFileRoot::save_pedestals(const char* fnam) {
    std::ofstream ofile(fnam);
    if(!ofile) {
        std::cout << "Could not open " << fnam << " to save pedestals." << std::endl;
        return;
    }

    // TODO: _nchan can be updated in an event by event basis
    //       while here we are assuming that it is the same
    //       for all the events
    ofile << _cmmd[0] << "\n";

    int i;
    for(i = 0; i < nchan(); i++) {
        ofile << _ped[i] << "\t" << _noise[i] << "\n";
    }
    ofile.close();
}

void DataFileRoot::load_pedestals(const char* fnam, bool show) {
    std::ifstream ifile(fnam);
    if(!ifile) {
        std::cout << "Could not open " << fnam << " to load pedestals." << std::endl;
        return;
    }

    ifile >> _cmmd[0] >> std::ws;

    int i;
    for(i = 0; i < max_nchan; i++) {
        if(ifile.eof())
            break;

        ifile >> _ped[i] >> std::ws >> _noise[i] >> std::ws;
        // std::cout << _ped[i] << std::endl;
        // _mask[i] = (_noise[i]>20. || _noise[i]<=0.);
    }
    ifile.close();
}

// This function processes the event, it subtracts pedestals and common mode from data
// and fills in the signal/noise ratio

void DataFileRoot::process_event(bool do_cmmd) {
    for(int i : _roi) {
        _signal[i] = (_data.data[i] - _ped[i] - _cmmd_roi) * _polarity;
        _sn[i] = _signal[i] / _noise[i];
    }
}

void DataFileRoot::calc_common_mode_signal() {
    int ip, i, n;
    double mean, st_dev, sig, signal_square, signal_sum, tmp;
    bool use_it;

    // Iterate common mode calculation three times to get better result
    for(ip = 0; ip < 3; ip++) {
        n = 0;
        signal_square = signal_sum = 0;
        // Use only channels in ROI
        for(int i : _roi) {
            use_it = true;
            sig = _data.data[i] - _ped[i];
            // In first iteration of calculation, mean ist not defined -> Use
            // all roi channels for cmmd calculation
            if(ip) {
                // Filter out all channels for which the deviation from the mean
                // signal is larger than 2.5 standard deviations
                tmp = fabs((sig - mean) / st_dev);
                use_it = (tmp < 2.5);
            }
            if(use_it) {
                // Counter how many channels have been used for the mean calculation
                n++;
                signal_sum += sig;
                signal_square += sig * sig;
            }
        }
        // Check if any channels have been used for calculation
        if(n > 0) {
            mean = signal_sum / n;
            st_dev = sqrt(signal_square / n - mean * mean);
        }
        // std::cout << "Iteration " << ip << ": Mean: " << mean << " Deviation: " << st_dev << std::endl;
    }
    _cmmd_roi = mean;
    _cnoise_roi = st_dev;
}

unsigned int DataFileRoot::clock_counter() const {
    return _data.clock;
}

double DataFileRoot::time() const {
    return 0.0;
}
double DataFileRoot::temp() const {
    return 0.0;
}

DataFileRoot* DataFileRoot::OpenFile(const char* nam, const char* pedfile, const char* gainfile) {
    struct stat sb;
    if(stat(nam, &sb) == -1)
        return 0;
        
        
        

std::ifstream ifile(nam);
if(!ifile)
    return 0;

char buf[5] = {'\0'};
ifile.read(buf, 4);
ifile.close();
std::string idf(buf + 1);

if(idf == "HDF")
    return new HDFRoot(nam, pedfile, gainfile);
else
    return new AsciiRoot(nam, pedfile, gainfile);

}
