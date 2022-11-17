/*
 * HDFRoot.cc
 *
 *  Created on: Mar 5, 2014
 *      Author: lacasta
 */

#ifdef HAVE_CONFIG_H
#ifdef G_OS_WIN32
#include "config_win.h"
#else
#include <config.h>
#endif
#endif

#include <ctime>
#include <hdf5.h>
#include <iostream>
#include "HDFRoot.h"
#include "core/utils/log.h"

struct ScanDef {
    enum ScanTypes { Unknown, Charge, Time, Laser };
    unsigned int type;
    unsigned int from;
    unsigned int to;
    unsigned int npts;
    unsigned int nevts;
    unsigned int charge;
    unsigned int delay;
};

struct Header {
    int run_type;
    unsigned short nchips;
    unsigned short chip_mask;
    unsigned short firmware;
    float pedestals[256];
    float noise[256];
    struct ScanDef scan;
    struct tm time;
};

struct DataScanPoint {
    float value; // Value of current scan value
    int start;   // Event index where this scan point starts
    int end;     // Event index where this scan point finishes
};

struct AlibavaData {
    float temp;
    float time;
    unsigned int clock;
    unsigned short header[32];
    unsigned short data[256];
};

struct HDFRootPrivate {
    hid_t fileid;
    hid_t headerID; /* The DataSet ID for the event headers */
    hid_t dataID;   /* The DataSet ID for data */
    hid_t tempID;   /* */
    hid_t timeID;
    hid_t clockID;
    hid_t scvalID;
    hid_t startID;
    hid_t endID;
    unsigned int nevts;
    unsigned int npoints;
    unsigned int ievt;
    unsigned int ipoint;
    unsigned int saved_evt, saved_point;
    AlibavaData data;
    DataScanPoint scan;
    HDFRootPrivate()
        : fileid(H5I_BADID), headerID(H5I_BADID), dataID(H5I_BADID), tempID(H5I_BADID), timeID(H5I_BADID),
          scvalID(H5I_BADID), startID(H5I_BADID), endID(H5I_BADID), nevts(0), npoints(0), ievt(0), ipoint(0), saved_evt(0),
          saved_point(0) {}
};

herr_t read_file_data(hid_t data_set, unsigned int offset, int ncols, void* dest);
herr_t read_file_data(hid_t data_set, unsigned int offset, int ncols, void* dest) {
    /*
     * Create the memory data space.
     * The one we want to add
     */
    int rank = (ncols > 1 ? 2 : 1);
    hsize_t mdims[2] = {1, static_cast<hsize_t>(ncols)};
    hid_t m_sid = H5Screate_simple(rank, mdims, nullptr);

    /*
     * ...and select the hyperslab in the file dataset
     */
    hsize_t start[2] = {static_cast<hsize_t>(offset), 0};
    hsize_t count[2] = {1, static_cast<hsize_t>(ncols)};
    hid_t sid = H5Dget_space(data_set); /* A copy of data space for writing */
    H5Sselect_hyperslab(sid, H5S_SELECT_SET, start, nullptr, count, nullptr);

    /*
     * Now read
     */
    hid_t tid = H5Dget_type(data_set);
    herr_t ret = H5Dread(data_set, tid, m_sid, sid, H5P_DEFAULT, dest);

    H5Sclose(m_sid);
    H5Sclose(sid);
    H5Tclose(tid);

    return ret;
}

HDFRoot::HDFRoot(const char* nam, const char* pedfile, const char* gainfile)
    : DataFileRoot(nam, pedfile, gainfile), priv(nullptr) {
    priv = new HDFRootPrivate;
    if(nam)
        open(nam);
}

HDFRoot::~HDFRoot() {
    if(valid())
        close();
}

bool HDFRoot::valid() const {
    return priv->fileid != H5I_BADID;
}

unsigned int HDFRoot::nevents() const {
    return priv->nevts;
}

void HDFRoot::open(const char* name) {
    Header H;
    priv->fileid = H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT);
    priv->ievt = 0;
    priv->ipoint = 0;

    if(!valid())
        return;

    memset(&H, 0, sizeof(H));

    // Get the stuff in the header
    hid_t header = H5Gopen2(priv->fileid, "/header", H5P_DEFAULT);
    hid_t setup = H5Aopen(header, "setup", H5P_DEFAULT);
    hid_t dtype = H5Aget_type(setup);
    H5Aread(setup, dtype, &H);
    H5Aclose(setup);
    H5Tclose(dtype);

    _type = static_cast<RunType>(H.run_type);
    _t0 = mktime(&H.time);

    _nchips = H.nchips;
    _nchan = 128 * H.nchips;
    _chip_mask = H.chip_mask;
    _firmware = H.firmware;

    // Now get the scan definition (if any)
    hid_t scan = H5Gopen2(priv->fileid, "/scan", H5P_DEFAULT);
    if(H5Aexists(scan, "scan_definition")) {
        setup = H5Aopen(scan, "scan_definition", H5P_DEFAULT);
        dtype = H5Aget_type(setup);
        H5Aread(setup, dtype, &H.scan);
        H5Aclose(setup);
        H5Tclose(dtype);
    }
    _npoints = H.scan.npts;
    _from = H.scan.from;
    _to = H.scan.to;
    _nevts = H.scan.nevts;
    _charge = static_cast<int>(H.scan.charge);
    _delay = static_cast<int>(H.scan.delay);
    _scantype = static_cast<DataFileRoot::ScanType>(H.scan.type);
    if(_npoints > 0)
        _step = (_from - _to) / _npoints;
    else
        _step = 0;

    LOG(DEBUG) << "Run type: " << static_cast<int>(_type) << " acquired on " << ctime(&_t0);
    LOG(DEBUG) << _t0;

    // Get pedestals and noise
    hid_t dset = H5Dopen2(priv->fileid, "/header/pedestal", H5P_DEFAULT);
    read_file_data(dset, 0, _nchan, H.pedestals);
    for(int ic = 0; ic < _nchan; ic++)
        _ped[ic] = H.pedestals[ic];
    H5Dclose(dset);

    dset = H5Dopen2(priv->fileid, "/header/noise", H5P_DEFAULT);
    read_file_data(dset, 0, _nchan, H.noise);
    for(int ic = 0; ic < _nchan; ic++)
        _noise[ic] = H.noise[ic];
    H5Dclose(dset);

    /*
     * Get all the relevant datasets
     */
    priv->headerID = H5Dopen2(priv->fileid, "/events/header", H5P_DEFAULT);
    priv->dataID = H5Dopen2(priv->fileid, "/events/signal", H5P_DEFAULT);
    priv->tempID = H5Dopen2(priv->fileid, "/events/temperature", H5P_DEFAULT);
    priv->timeID = H5Dopen2(priv->fileid, "/events/time", H5P_DEFAULT);
    priv->clockID = H5Dopen2(priv->fileid, "/events/clock", H5P_DEFAULT);
    priv->scvalID = H5Dopen2(priv->fileid, "/scan/value", H5P_DEFAULT);
    priv->startID = H5Dopen2(priv->fileid, "/scan/start", H5P_DEFAULT);
    priv->endID = H5Dopen2(priv->fileid, "/scan/end", H5P_DEFAULT);

    hsize_t dims[2];
    hid_t dspace = H5Dget_space(priv->timeID);
    H5Sget_simple_extent_dims(dspace, dims, nullptr);
    H5Sclose(dspace);
    priv->nevts = static_cast<unsigned int>(dims[0]);

    dspace = H5Dget_space(priv->startID);
    H5Sget_simple_extent_dims(dspace, dims, nullptr);
    H5Sclose(dspace);
    priv->npoints = static_cast<unsigned int>(dims[0]);

    LOG(DEBUG) << "Nr. of evts. " << priv->nevts << " - Nr. of points " << priv->npoints << std::endl;

    /*
     * Now get the very first scan point
     */
    next_scan_point();

    priv->saved_evt = priv->ievt;
    priv->saved_point = priv->ipoint;
}

void HDFRoot::close() {
    if(valid()) {
        H5Dclose(priv->headerID);
        H5Dclose(priv->dataID);
        H5Dclose(priv->tempID);
        H5Dclose(priv->timeID);
        H5Dclose(priv->clockID);
        H5Dclose(priv->scvalID);
        H5Dclose(priv->startID);
        H5Dclose(priv->endID);
        H5Fclose(priv->fileid);
        priv->fileid = H5I_BADID;
    }
}

void HDFRoot::next_scan_point() {
    if(priv->ipoint < priv->npoints) {
        read_file_data(priv->scvalID, priv->ipoint, 1, &priv->scan.value);
        read_file_data(priv->startID, priv->ipoint, 1, &priv->scan.start);
        read_file_data(priv->endID, priv->ipoint, 1, &priv->scan.end);
        _data.value = priv->scan.value;
        priv->ipoint++;
    }
}
void HDFRoot::rewind() {
    priv->ievt = 0;
    priv->ipoint = 0;
    next_scan_point();
}

void HDFRoot::save() {
    if(priv->fileid) {
        priv->saved_evt = priv->ievt;
        priv->saved_point = priv->ipoint;
    }
}

void HDFRoot::restore() {
    if(priv->fileid) {
        priv->ievt = priv->saved_evt;
        priv->ipoint = priv->saved_point;
    }
}

int HDFRoot::read_data() {
    herr_t rc;
    if(priv->ievt >= priv->nevts)
        return -1;

    rc = read_file_data(priv->dataID, priv->ievt, _nchan, &_data.data); // priv->data.data);
    if(rc < 0)
        return rc;

    priv->ievt++;

    return 0;
}

int HDFRoot::read_event() {
    herr_t rc;
    if(priv->ievt >= priv->nevts)
        return -1;
    // End of file reached

    rc = read_file_data(priv->clockID, priv->ievt, 1, &priv->data.clock);
    if(rc < 0)
        // return rc;
        return 4;

    rc = read_file_data(priv->tempID, priv->ievt, 1, &priv->data.temp);
    if(rc < 0)
        // return rc;
        return 4;

    rc = read_file_data(priv->timeID, priv->ievt, 1, &priv->data.time);
    if(rc < 0)
        // return rc;
        return 4;

    rc = read_file_data(priv->headerID, priv->ievt, 16 * _nchips, &priv->data.header);
    if(rc < 0)
        // return rc;
        return 4;

    rc = read_file_data(priv->dataID, priv->ievt, _nchan, &_data.data); // priv->data.data);
    if(rc < 0)
        // return rc;
        return 4;

    int ij = 0;
    for(int ichip = 0; ichip < _nchips; ichip++) {
        for(int ih = 0; ih < 16; ih++, ij++)
            _header[ichip][ih] = priv->data.header[ij];
    }
    priv->ievt++;
    if(priv->ievt >= static_cast<unsigned int>(priv->scan.end))
        next_scan_point();

    return 1;
}

void HDFRoot::get_scan_values(short& delay, short& charge) {
    if(scan_type() == DataFileRoot::Charge) {
        charge = static_cast<short>(value() / 1024);
        delay = static_cast<short>(_delay);
    } else {
        delay = static_cast<short>(value());
        charge = static_cast<short>(_charge / 1024);
    }
}

double HDFRoot::time() const {
    return priv->data.time;
}
double HDFRoot::temp() const {
    return priv->data.temp;
}
unsigned int HDFRoot::clock_counter() const {
    return priv->data.clock;
}
