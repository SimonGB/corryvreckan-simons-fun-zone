#ifndef DATAFILEROOT_H_
#define DATAFILEROOT_H_

#include <vector>
#include "Data.h"
// #include "Hit.h"
// #include "ChanList.h"
#include <TH1.h>
#include <TH2.h>
#include <ctime>

// maximum number of channels
#define MAX_NCHAN 256

/**
 * class DataFileRoot
 *
 * This class defines the interface to read the Alibava data from file
 * and also defines a number of tools to help analyzing the data
 */
class DataFileRoot {
public:
    enum ScanType { UnknownScan, Charge, Time, LaserScan };
    enum RunType { UnknownRun = 0, Calibration = 1, LaserSync, Laser, RadSource, Pedestal, ChargeScan, LastRType };

protected: // This is ugly but comfortable
    RunType _type;
    time_t _t0;
    int _nchips;
    int _chip_mask;
    int _firmware;
    ScanType _scantype;
    unsigned int _npoints;
    unsigned int _from;
    unsigned int _to;
    unsigned int _step;
    unsigned int _nevts;
    int _charge;
    int _delay;
    int _nchan; // current number of channels
    double _seedcut;
    double _neighcut;
    unsigned short _header[2][16];
    double _ped[MAX_NCHAN];
    double _noise[MAX_NCHAN];
    double _signal[MAX_NCHAN];
    double _sn[MAX_NCHAN];
    double _cmmd[2];
    double _cnoise[2];
    double _gain[MAX_NCHAN];
    double _average_gain;
    bool _mask[MAX_NCHAN];
    int _version;
    int _polarity;
    double _t1, _t2;

    // not happy with this, would be nicer as a static
    std::vector<unsigned int> _roi;
    double _cmmd_roi;
    double _cnoise_roi;

    EventDataBlock _data;
    static std::string _idf;
    double _mean_temp_pedestal = std::numeric_limits<double>::quiet_NaN();

protected:
    void reset_data();

public:
    /**
     * The constructor needs three parameters
     *
     * nam      - the path of the input file
     * pedfile  - the path of a pedestal file
     * gainfile - the path of a gain file (translates from ADC to charge)
     *
     * The pedestal file is an ASCII file where each line contains
     * the pedestal and noise value of a channel. There should be
     * as many lines as channels.
     *
     * A gain file is an ASCII file where each line contains the
     * channel number and the gain (conversion from ADC to charge
     * in electrons).
     *
     */
    DataFileRoot(const char* nam = nullptr, const char* pedfile = nullptr, const char* gainfile = nullptr);

    /**
     * Destructor
     */
    virtual ~DataFileRoot();

    /**
     * Opens a file. It can be either an Alibava binary file or a HDF file.
     * It will return the corresponding handler. The parameters are the same as
     * for the constructor:
     */
    static DataFileRoot* OpenFile(const char* nam = nullptr, const char* pedfile = nullptr, const char* gainfile = nullptr);

    /**
     * These are the file operation methods that need to be
     * implemented
     *
     */

    /**
     * Sets the number of channels and the data in the case
     * of non "standard" values. If data==0, then only the number
     * of channels is changed
     */
    void set_data(int nchan, const unsigned short* data = nullptr);

    /**
     * Tells if the data stream is valid.
     *
     * @return a boolean
     */
    virtual bool valid() const = 0;

    // all of the below functions are implemented in AsciiRoot.h. I dont know why they work here

    /*
     * Scan information
     */
    // return the scan type
    int scan_type() const { return _scantype; }

    /**
     * Opens a new file.
     * \param name the path of the file to open
     *
     * The success of the operation can be checked with valid
     */
    virtual void open(const char* name) = 0;

    /// Closes the file.
    virtual void close() = 0;

    /// Moves the read pointer to the beginning of the file
    virtual void rewind() = 0;

    /// Save the current state of the file
    virtual void save(){};

    /// Restore the last saved state of the file
    virtual void restore() { rewind(); }

    /**
     * Reads all the information about a given event.
     *
     * \returns 0 in case of success.
     */
    virtual int read_event() = 0;

    /**
     *  Reads only the ASIC data of a given event
     *  \returns 0 in case of success
     */
    virtual int read_data() = 0;

    virtual void check_point(unsigned int, const char*){};
    virtual void new_file(unsigned int, const char*) {}
    virtual void start_of_run(unsigned int, const char*) {}
    virtual void end_of_run(unsigned int, const char*) {}
    virtual void new_data_block(unsigned int, const char*){};

    // The data format version
    int version() const { return _version; }

    int get_polarity() const { return _polarity; }
    void set_polarity(int x) { _polarity = (x < 0 ? -1 : 1); }

    /*
     * General information
     */
    // returns the number of channels
    int nchan() const { return _nchan; }

    // returns the run type
    int type() const { return _type; }

    // returns the date of the run
    char* date() const { return ctime(&_t0); }

    // returns the pedestal value of channel i
    double ped(unsigned int i) const { return _ped[i]; }

    // returns the noise value of channel i

    double noise(unsigned int i) const { return _noise[i]; }

    // returns the signal over noise ratio of channel i
    double sn(unsigned int i) const { return _sn[i]; }

    // returns the common mode  of chip i
    double get_cmmd(unsigned int i) const { return _cmmd[i]; }

    // returns the common mode noise of chip i
    double get_cnoise(unsigned int i) const { return _cnoise[i]; }

    /*
     * Event specific information
     */

    // returns the raw data of channel i
    unsigned short data(unsigned int i) const { return _data.data[i]; }

    // returns the signal value of channel i
    double signal(unsigned int i) const { return _signal[i] / _gain[i]; }

    double ADC_signal(unsigned int i) const { return _signal[i]; }

    // return the scan value
    double value() const { return _data.value; }

    // return the clock counter
    virtual unsigned int clock_counter() const;

    // return the TDC time
    virtual double time() const;

    // return the temperature
    virtual double temp() const;

    unsigned short get_header(int ichip, int ibit) { return _header[ichip][ibit]; }

    // returns the gain of channel i
    double get_gain(int i) const { return _gain[i]; }

    // returns the average gain
    double gain() const { return _average_gain; }

    /*
     * Analysis
     *
     * To find clusters we follow a very simple algorithm in which
     * we look among the "unused" channels a possible cluster seed
     * by requiring that it has a signal over noise above the
     * seed_cut value. If more than one channel fulfills this
     * requirement we take the one with the highest value. We then
     * add the neighbor channels while their signal over noise is
     * above the neighbor cut. Those channels are marked as
     * "used" and the process starts again until we do not find a
     * cluster seed.
     *
     */

    // returns the seed cut
    double seed_cut() const { return _seedcut; }

    // returns the neighbour cut
    double neigh_cut() const { return _neighcut; }

    // Set the seed and neighbor cuts
    void set_cuts(double seed, double neighbor) {
        _seedcut = seed;
        _neighcut = neighbor;
    }

    // Define the time interval where the pulse has the peak
    void set_timecut(double t1, double t2);
    bool valid_time(double time) const;

    // Alternative functions
    virtual void compute_pedestals_alternative();
    virtual void compute_cmmd_alternative();

    // Subtract pedestals and compute common mode and s/n
    void process_event(bool do_cmmd = true);

    // Save/load pedestals in/from a file
    void save_pedestals(const char* fnam);
    void load_pedestals(const char* fnam, bool show = false);

    int chip_mask() const { return _chip_mask; }

    int firmware() const { return _firmware; }

    int nchips() const { return _nchips; }

    void set_ROI(std::vector<unsigned int> bounds);

    std::vector<unsigned int> get_ROI() const { return _roi; }

    void calc_common_mode_signal();

    double get_mean_pedestal_temp() const { return _mean_temp_pedestal; }
};

#endif /* DATAFILEROOT_H_ */
