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
#include <config.h>
#endif

#include "AsciiRoot.h"
#include "DataFileRoot.h"
#include "HDFRoot.h"

#ifdef __APPLE__
#define sighandler_t sig_t
#endif

bool _A_do_run = true;
void _A_got_intr(int);
void _A_got_intr(int) { _A_do_run = false; }
std::string DataFileRoot::_idf;

DataFileRoot::DataFileRoot(const char*, const char* pedfile, const char*)
    : _nchan(MAX_NCHAN), _seedcut(5.), _neighcut(3.), _average_gain(1.), _version(2), _polarity(1), _t1(0.0), _t2(99999.),
      _roi({}) {
    int i;

    for(i = 0; i < MAX_NCHAN; i++) {
        _ped[i] = 0.;
        _gain[i] = 1.;
        _noise[i] = 1.;
        _data.data[i] = 0;
        _mask[i] = false;
    }

    //  Cannot call open in the constructor since overrides are not
    //  available in the constructor

    if(pedfile)
        load_pedestals(pedfile);

    // TODO: this loads a "fixed" name file. Be more general...
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

void DataFileRoot::reset_data() { memset(&_data, 0, sizeof(_data)); }

void DataFileRoot::set_ROI(std::vector<unsigned int> roi) { _roi = roi; }

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

bool DataFileRoot::valid_time(double tim) const { return (tim >= _t1 && tim <= _t2); }

void DataFileRoot::compute_pedestals_alternative() {
    int mxevts = 10000000;

    int ievt;
    std::vector<double> pedestal_data[MAX_NCHAN];
    double pedestal_average[MAX_NCHAN];
    double pedestal_stdev[MAX_NCHAN];

    if(!valid())
        return;

    for(unsigned int i : _roi)
        _ped[i] = _noise[i] = 0.;

    for(ievt = 0; read_data() == 1 && ievt < mxevts; ievt++) {
        for(unsigned int i : _roi) {
            pedestal_data[i].push_back(_data.data[i]);
        }
    }

    for(unsigned int i : _roi) {
        pedestal_average[i] = std::accumulate(pedestal_data[i].begin(), pedestal_data[i].end(), 0.0);

        pedestal_average[i] = pedestal_average[i] / static_cast<double>(pedestal_data[i].size());

        pedestal_stdev[i] =
            std::inner_product(pedestal_data[i].begin(), pedestal_data[i].end(), pedestal_data[i].begin(), 0.0);
        pedestal_stdev[i] = std::sqrt(pedestal_stdev[i] / static_cast<double>(pedestal_data[i].size()) -
                                      pedestal_average[i] * pedestal_average[i]);

        _ped[i] = pedestal_average[i];
        _noise[i] = pedestal_stdev[i];
    }

    rewind();
}

void DataFileRoot::compute_cmmd_alternative() {
    int mxevts = 10000000;
    int nEvents = 0;

    int ievt;

    double event_bias = 0; // common mode noise per pedestal events
    double cmn = 0;        // common mode noise averaged over all pedestal events

    std::vector<double> corrected_pedestal_data[MAX_NCHAN];
    double pedestal_average[MAX_NCHAN];
    double pedestal_stdev[MAX_NCHAN];

    double temp_sum = 0;
    double temperature;

    if(!valid())
        return;

    for(ievt = 0; read_data() == 1 && ievt < mxevts; ievt++) {
        nEvents++;
        double event_sum = 0;

        // Get temp
        if(_idf == "HDF") {
            temperature = _data.temp;
        } else {
            temperature = 0.12 * _data.temp - 39.8;
        }
        temp_sum += temperature;

        for(unsigned int i : _roi) {
            event_sum += (_data.data[i] - _ped[i]);
        }
        event_bias = event_sum / static_cast<double>(_roi.size());
        cmn += event_bias;
        for(unsigned int i : _roi) {
            corrected_pedestal_data[i].push_back(_data.data[i] - event_bias);
        }
    }

    _cmmd[0] = cmn / nEvents;
    _mean_temp_pedestal = temp_sum / nEvents;

    for(unsigned int i : _roi) {

        pedestal_average[i] = std::accumulate(corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), 0.0);
        pedestal_average[i] = pedestal_average[i] / static_cast<double>(corrected_pedestal_data[i].size());

        pedestal_stdev[i] = std::inner_product(
            corrected_pedestal_data[i].begin(), corrected_pedestal_data[i].end(), corrected_pedestal_data[i].begin(), 0.0);
        pedestal_stdev[i] = std::sqrt(pedestal_stdev[i] / static_cast<double>(corrected_pedestal_data[i].size()) -
                                      pedestal_average[i] * pedestal_average[i]);

        _ped[i] = pedestal_average[i];
        _noise[i] = pedestal_stdev[i];
    }

    rewind();
}

void DataFileRoot::save_pedestals(const char* fnam) {
    std::ofstream ofile(fnam);
    if(!ofile) {
        std::cout << "Could not open " << fnam << " to save pedestals." << std::endl;
        return;
    }

    ofile << _cmmd[0] << "\n";

    int i;
    for(i = 0; i < nchan(); i++) {
        ofile << _ped[i] << "\t" << _noise[i] << "\n";
    }
    ofile.close();
}

void DataFileRoot::load_pedestals(const char* fnam, bool) {
    std::ifstream ifile(fnam);
    if(!ifile) {
        std::cout << "Could not open " << fnam << " to load pedestals." << std::endl;
        return;
    }

    ifile >> _cmmd[0] >> std::ws;

    int i;
    for(i = 0; i < MAX_NCHAN; i++) {
        if(ifile.eof())
            break;

        ifile >> _ped[i] >> std::ws >> _noise[i] >> std::ws;
    }
    ifile.close();
}

// This function processes the event, it subtracts pedestals and common mode from data
// and fills in the signal/noise ratio

void DataFileRoot::process_event(bool) {
    for(auto i : _roi) {
        _signal[i] = (_data.data[i] - _ped[i] - _cmmd_roi) * _polarity;
        _sn[i] = _signal[i] / _noise[i];
    }
}

void DataFileRoot::calc_common_mode_signal() {
    int ip, n;
    double mean = 0, st_dev = 0, sig, signal_square, signal_sum, tmp;
    bool use_it;

    // Iterate common mode calculation three times to get better result
    for(ip = 0; ip < 3; ip++) {
        n = 0;
        signal_square = signal_sum = 0;
        // Use only channels in ROI
        for(auto i : _roi) {
            use_it = true;
            sig = _data.data[i] - _ped[i];
            // In first iteration of calculation, mean is not defined -> Use
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
    }
    _cmmd_roi = mean;
    _cnoise_roi = st_dev;
}

unsigned int DataFileRoot::clock_counter() const { return _data.clock; }

double DataFileRoot::time() const { return 0.0; }
double DataFileRoot::temp() const { return 0.0; }

DataFileRoot* DataFileRoot::OpenFile(const char* nam, const char* pedfile, const char* gainfile) {
    struct stat sb;
    if(stat(nam, &sb) == -1)
        return nullptr;

    std::ifstream ifile(nam);
    if(!ifile)
        return nullptr;

    char buf[5] = {'\0'};
    ifile.read(buf, 4);
    ifile.close();
    _idf = std::string(buf + 1);

    if(_idf == "HDF")
        return new HDFRoot(nam, pedfile, gainfile);
    else
        return new AsciiRoot(nam, pedfile, gainfile);
}
