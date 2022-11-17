#ifndef __Alibava_AsciiRoot_h__
#define __Alibava_AsciiRoot_h__

#include <string>
#include "DataFileRoot.h"

/**
 * This is the class that reads the data files
 */

struct AsciiRootPriv;
class AsciiRoot : public DataFileRoot {
public:
    typedef std::vector<std::string> XtraValues;
    enum BlockType { NewFile = 0, StartOfRun, DataBlock, CheckPoint, EndOfRun };

private:
    AsciiRootPriv* priv;
    long data_start;
    long saved_state;
    XtraValues _xtra; // extra values from header

public:
    AsciiRoot(const char* nam = nullptr, const char* pedfile = nullptr, const char* gainfile = nullptr);
    virtual ~AsciiRoot();

    bool valid() const;
    int nevents() const;
    void open(const char* name);
    void save();
    void restore();
    void close();
    void rewind();
    int read_event();
    /**
     * AsciiRoot does not implement an independent access of ASIC data
     * and therefore needs to read the whole event data.
     */
    int read_data() { return read_event(); }

    double time() const;
    double temp() const;

    void get_scan_values(short& delay, short& charge);

    size_t nxtra() const { return _xtra.size(); }
    const std::string xtra(size_t i) const { return _xtra[i]; }
    void add_xtra(const std::string& x) { _xtra.push_back(x); }
    void add_xtra(const char* x) { _xtra.push_back(x); }
};

#endif
