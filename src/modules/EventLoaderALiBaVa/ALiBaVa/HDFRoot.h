/*
 * HDFRoot.h
 *
 *  Created on: Mar 5, 2014
 *      Author: lacasta
 */

#ifndef HDFROOT_H_
#define HDFROOT_H_

#include "DataFileRoot.h"

struct HDFRootPrivate;

class HDFRoot : public DataFileRoot {
private:
    HDFRootPrivate* priv;

    void next_scan_point();

public:
    HDFRoot(const char* nam = nullptr, const char* pedfile = nullptr, const char* gainfile = nullptr);
    virtual ~HDFRoot();

    bool valid() const;
    unsigned int nevents() const;
    void open(const char* name);
    void close();
    void rewind();
    void save();
    void restore();
    int read_event();
    int read_data();

    unsigned int clock_counter() const;
    double time() const;
    double temp() const;
    void get_scan_values(short& delay, short& charge);
};

#endif /* HDFROOT_H_ */
