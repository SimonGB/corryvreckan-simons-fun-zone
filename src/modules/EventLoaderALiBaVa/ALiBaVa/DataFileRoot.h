#ifndef DATAFILEROOT_H_
#define DATAFILEROOT_H_

#include <vector>
#include "Data.h"
#include "Hit.h"
#include "ChanList.h"
#include <ctime>
#include <TH1.h>
#include <TH2.h>

/**
 * class DataFileRoot
 *
 * This class defines the interface to read the Alibava data from file
 * and also defines a number of tools to help analyzing the data
 */
class DataFileRoot
{
    public:
        enum ScanType { UnknownScan, Charge, Time, LaserScan };
        enum RunType
        {
            UnknownRun=0,
            Calibration=1,
            LaserSync,
            Laser,
            RadSource,
            Pedestal,
            ChargeScan,
            LastRType
        };

    protected: // This is ugly but comfortable
        static const int max_nchan=256;

        RunType _type;
        time_t _t0;
        int _nchips;
        int _chip_mask;
        int _firmware;
        ScanType _scantype;
        int _npoints;
        int _from;
        int _to;
        int _step;
        int _nevts;
        int _charge;
        int _delay;
        int _nchan; // current number of channels
        double _seedcut;
        double _neighcut;
        unsigned short _header[2][16];
        double _ped[max_nchan];
        double _noise[max_nchan];
        double _signal[max_nchan];
        double _sn[max_nchan];
        double _cmmd[2];
        double _cnoise[2];
        double _gain[max_nchan];
        double _average_gain;
        bool   _mask[max_nchan];
        int     _version;
        int     _polarity;
        double _t1, _t2;
        HitList _hits;
        std::vector<ChanList> chan_list;
        EventDataBlock _data;

    protected:
        void reset_data();

    public:
        void set_data(int i, unsigned short x) { _data.data[i] = x; }

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
        DataFileRoot(const char *nam=0, const char *pedfile=0, const char *gainfile=0);

        /**
         * Destructor
         */
        virtual ~DataFileRoot();

        /**
         * Opens a file. It can be either an Alibava binary file or a HDF file.
         * It will return the corresponding handler. The parameters are the same as
         * for the constructor:
         */
        static DataFileRoot *OpenFile(const char *nam=0, const char *pedfile=0, const char *gainfile=0);

        /**
         * These are the file operation methods that need to be
         * implemented
         *
         */

         /**
          * Gives the number of events.
          *
          * @return an int
          */
         virtual int nevents() const =0;

        /**
         * Tells if the data stream is valid.
         *
         * @return a boolean
         */
        virtual bool valid() const =0;

        /**
         * Opens a new file.
         * \param name the path of the file to open
         *
         * The success of the operation can be checked with valid
         */
        virtual void open(const char *name)=0;

        /// Closes the file.
        virtual void close()=0;

        /// Moves the read pointer to the beginning of the file
        virtual void rewind()=0;

        /// Save the current state of the file
        virtual void save() {};

        /// Restore the last saved state of the file
        virtual void restore() { rewind(); }

        /**
         * Reads all the information about a given event.
         *
         * \returns 0 in case of success.
         */
        virtual int read_event()=0;

        /**
         *  Reads only the ASIC data of a given event
         *  \returns 0 in case of success
         */
        virtual int read_data()=0;

        virtual void check_point(int, const char *) {};
        virtual void new_file(int, const char *) {}
        virtual void start_of_run(int, const char *) {}
        virtual void end_of_run(int, const char *) {}
        virtual void new_data_block(int, const char *) {};

        // The data format version
        int version() const { return _version; }


        int polarity() const { return _polarity; }
        void polarity(int x) { _polarity = ( x<0 ? -1 : 1); }

        /*
         * Sets the number of channels and the data in the case
         * of non "standard" values. If data==0, then only the number
         * of channels is changed
         */
        void set_data(int nchan, const unsigned short *data=0);


        /*
         * General information
         */
        // returns the number of channels
        int nchan() const { return _nchan; }

        // returns the run type
        int type() const
        {
            return _type;
        }

        // returns the date of the run
        char *date() const
        {
            return ctime(&_t0);
        }

        // returns the pedestal value of channel i
        double ped(int i) const
        {
            return _ped[i];
        }

        // double ped(int i) const
        // {
        //     return _ped[i]/_gain[i];
        // }

        // returns the noise value of channel i
        
        double noise(int i) const
        {
            return _noise[i];
        }

        // double noise(int i) const
        // {
        //     return _noise[i]/_gain[i];
        // }

        // returns the signal over noise ratio of channel i
        double sn(int i) const
        {
            return _sn[i];
        }

        // returns the common mode  of chip i
        double get_cmmd(int i) const
        {
            return _cmmd[i];
        }

        // returns the common mode noise of chip i
        double get_cnoise(int i) const
        {
            return _cnoise[i];
        }


        /*
         * Event specific information
         */

        // returns the raw data of channel i
        unsigned short data(int i) const
        {
            return _data.data[i];
        }

        // returns the signal value of channel i
        double signal(int i) const
        {
            return _signal[i]/_gain[i];
        }

        double ADC_signal(int i) const
        {
            return _signal[i];
        }

        // return the scan value
        double value() const
        {
            return _data.value;
        }

        // return the clock counter
        virtual unsigned int clock_counter() const;

        // return the TDC time
        virtual double time() const;

        // return the temperature
        virtual double temp() const;


        /*
         * Scan information
         */
        // return the scan type
        int scan_type() const { return _scantype; }

        // number of points
        int npts() const
        {
            return _npoints;
        }

        // first value of the scan
        int from() const
        {
            return _from;
        }

        // last value of the scan
        int to() const
        {
            return _to;
        }

        // the step from one point to the next
        int step() const
        {
            return _step;
        }

        // Number of events acquired per scan point
        int nevts() const
        {
            return _step;
        }


        // returns the value of the scan variables
        virtual void get_scan_values(short &delay, short &charge) = 0;

        unsigned short get_header(int ichip, int ibit) { return _header[ichip][ibit]; }


        /*
         * Event hit information
         */

        // Adds a new hit
        void add_hit(const Hit &h)
        {
            _hits.push_back(h);
        }

        /*
         * Iterators
         */
        HitList::iterator begin()
        {
            return _hits.begin();
        }
        HitList::iterator end()
        {
            return _hits.end();
        }

        // number of hits
        int nhits() const
        {
            return _hits.size();
        }

        // check if the hit list is empty
        bool empty() const
        {
            return _hits.empty();
        }

        // retrieve i-th hit
        const Hit &hit(int i) const
        {
            return _hits[i];
        }

        // Empty the hit list
        void clear()
        {
            _hits.clear();
        }

        // replace the hit list by the given
        void set_hit_list(const HitList &L) { _hits = L; }

        // returns the gain of channel i
        double get_gain(int i) const
        {
            return _gain[i];
        }

        // returns the average gain
        double gain() const
        {
            return _average_gain;
        }



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
        double seed_cut() const
        {
            return _seedcut;
        }

        // returns the neightbour cut
        double neigh_cut() const
        {
            return _neighcut;
        }

        // Set the seed and neighbor cuts
        void set_cuts(double seed, double neighbor)
        {
            _seedcut = seed;
            _neighcut = neighbor;
        }

        // Define the time interval where the pulse has the peak
        void set_timecut(double t1, double t2);
        bool valid_time(double time) const;

        // Returns a histogram with the pedestals
        TH1F *show_pedestals(const int lowerChannel, const int upperChannel, bool corrected);

        // Returns a histogram with the noise
        TH1F *show_noise(const int lowerChannel, const int upperChannel, bool corrected);

        // Computes pedestals from the data of the file
        virtual TH2 *compute_pedestals(int mxevts=-1, bool do_cmmd=true);
        virtual void compute_pedestals_fast(int mxevts = -1, double ped_weight=0.01, double noise_weight=0.001);

        // Alternative functions
        virtual void compute_pedestals_alternative(const int lowerChannel, const int upperChannel);
        virtual void compute_cmmd_alternative(const int lowerChannel, const int upperChannel);

        // Computes common mode,
        void common_mode();
        void common_mode(ChanList &C, bool correct=false);


        // Subtract pedestals and compute common mode and s/n
        void process_event(bool do_cmmd=true);

        // find clusters
        void find_clusters(int ichip=-1);
        void find_clusters(ChanList &C);

        // Save/load pedestals in/from a file
        void save_pedestals(const char *fnam);
        void load_pedestals(const char *fnam, bool show=false);

        // Loads a gain file
        void load_gain(const char *fnam);

        // Loads a channel mask
        void load_masking(const char *fnam);
        void spy_data(bool with_signal=false, double t0=0.0, double t1=0.0, int nevt=1);

        int  n_channel_list() const { return chan_list.size(); }
        void add_channel_list(const ChanList &C);
        void clear_channel_lists() { chan_list.clear(); }
        ChanList get_channel_list(int i) const { return chan_list[i]; }

        int chip_mask() const
        {
            return _chip_mask;
        }

        int firmware() const
        {
            return _firmware;
        }

        int nchips() const
        {
            return _nchips;
        }
};

// Return true if file is an ASCII text file
bool is_text(const char *);

#endif /* DATAFILEROOT_H_ */
