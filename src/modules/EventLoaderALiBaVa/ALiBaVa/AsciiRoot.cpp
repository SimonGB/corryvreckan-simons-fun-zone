#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "AsciiRoot.h"
#include "core/utils/log.h"
#include "utils.h"

struct AsciiRootPriv {
    std::ifstream* ifile;

    AsciiRootPriv() : ifile(nullptr) {}
};

// decodes the header and returns a vector with the integers found
std::vector<unsigned int> decode_header(const std::string& h, AsciiRoot::XtraValues& xtra);
std::vector<unsigned int> decode_header(const std::string& h, AsciiRoot::XtraValues& xtra) {
    std::vector<unsigned int> vout;
    std::istringstream istr(h);
    char* endptr;
    char buf[256];
    long val;

    xtra.clear();

    while(istr) {
        istr.getline(buf, sizeof(buf), ';');
        if(!istr)
            break;

        errno = 0;
        val = strtol(buf, &endptr, 0);

        if((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0)) {
            std::string sval(buf), sout;
            sout = trim_str(sval);
            if(!sout.empty())
                xtra.push_back(sout);
        } else if(endptr == buf || *endptr != '\0') {
            std::string sval(buf), sout;
            sout = trim_str(sval);
            if(!sout.empty())
                xtra.push_back(sout);
        } else {
            vout.push_back(static_cast<unsigned int>(atoi(buf)));
        }
    }
    return vout;
}

AsciiRoot::AsciiRoot(const char* nam, const char* pedfile, const char* gainfile)
    : DataFileRoot(nam, pedfile, gainfile), priv(nullptr) {
    priv = new AsciiRootPriv;

    if(nam)
        open(nam);
}
AsciiRoot::~AsciiRoot() {
    if(priv->ifile) {
        priv->ifile->close();
    }
    delete priv;
}

int AsciiRoot::nevents() const {
    LOG(ERROR) << "Total number of events hasn't been properly implemented yet for binary files.";
    return -1;
}

bool AsciiRoot::valid() const {
    return (priv->ifile != nullptr);
}

void AsciiRoot::open(const char* name) {
    priv->ifile = new std::ifstream(name);
    if(!(*priv->ifile)) {
        std::cout << "Could not open data file: " << name << std::endl;
        delete priv->ifile;
        priv->ifile = nullptr;
        return;
    }
    std::string header;
    unsigned int ic, lheader;
    char c;
    priv->ifile->read(reinterpret_cast<char*>(&_t0), sizeof(time_t));
    priv->ifile->read(reinterpret_cast<char*>(&_type), sizeof(int));
    if(_type > LastRType) {
        priv->ifile->seekg(0, std::ios::beg);
        priv->ifile->read(reinterpret_cast<char*>(&_t0), sizeof(int));
        priv->ifile->read(reinterpret_cast<char*>(&_type), sizeof(int));
    }

    priv->ifile->read(reinterpret_cast<char*>(&lheader), sizeof(unsigned int));
    for(ic = 0; ic < 80; ic++) {
        priv->ifile->read(&c, sizeof(char));
        header.append(1, c);
    }
    header = trim_str(header);

    if(header[0] != 'V' && header[0] != 'v') {
        _version = 0;
    } else {
        _version = header[1] - '0';
        header = header.substr(5);
    }

    auto param = decode_header(header, _xtra);
    priv->ifile->read(reinterpret_cast<char*>(_ped), MAX_NCHAN * sizeof(double));
    priv->ifile->read(reinterpret_cast<char*>(_noise), MAX_NCHAN * sizeof(double));
    switch(_type) {
    case 1: // calibration
    case 2: // laser sync
        _npoints = param[0];
        _from = param[1];
        _to = param[2];
        _step = param[3];
        break;
    case 3: // laser run
    case 4: // source run
    case 5: // pedestal run
        if(param.empty())
            _nevts = 100000;
        else
            _nevts = param[0];
        _npoints = _from = _to = _step = 0;
        break;
    default:
        break;
    }
    _scantype = UnknownScan;
    if(type() == Calibration) {
        if(nxtra())
            _scantype = Time;
        else
            _scantype = Charge;
    } else if(type() == LaserSync)
        _scantype = LaserScan;

    data_start = priv->ifile->tellg();
    saved_state = data_start;
}

void AsciiRoot::save() {
    if(priv->ifile)
        saved_state = priv->ifile->tellg();
}

void AsciiRoot::restore() {
    if(priv->ifile) {
        priv->ifile->seekg(saved_state, std::ios::beg);
    }
}

void AsciiRoot::rewind() {
    if(priv->ifile) {
        priv->ifile->clear();
        priv->ifile->seekg(data_start, std::ios::beg);
    }
}

void AsciiRoot::close() {
    if(priv->ifile) {
        priv->ifile->close();
        delete priv->ifile;
        priv->ifile = nullptr;
    }
}

void AsciiRoot::get_scan_values(short& delay, short& charge) {
    delay = static_cast<short>(int(value()) >> 16);
    charge = static_cast<short>(int(value()) & 0xffff);
}

int AsciiRoot::read_event() {
    if(priv->ifile) {
        unsigned int header, size, user = 0, code = 0;
        char* block_data = nullptr;
        if(_version) {
            do {
                do {
                    priv->ifile->read(reinterpret_cast<char*>(&header), sizeof(unsigned int));
                    if(priv->ifile->bad() || priv->ifile->eof())
                        return -1;

                    code = (header >> 16) & 0xFFFF;
                } while(code != 0xcafe);

                code = header & 0x0fff;
                user = header & 0x1000;
                switch(code) {
                case NewFile:
                    priv->ifile->read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                    block_data = new char[size];
                    priv->ifile->read(block_data, size);
                    new_file(size, block_data);
                    break;
                case StartOfRun:
                    priv->ifile->read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                    block_data = new char[size];
                    priv->ifile->read(block_data, size);
                    start_of_run(size, block_data);
                    break;
                case DataBlock:
                    priv->ifile->read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                    if(user) {
                        reset_data();
                        block_data = new char[size];
                        priv->ifile->read(block_data, size);
                        new_data_block(size, block_data);
                    } else {
                        if(_version == 1) {
                            priv->ifile->read(reinterpret_cast<char*>(&_data), sizeof(EventData));
                            for(int ii = 0; ii < 2; ii++)
                                memset(_header[ii], 0, 16 * sizeof(unsigned short));

                        } else {
                            priv->ifile->read(reinterpret_cast<char*>(&_data.value), sizeof(double));
                            if(_version > 2)
                                priv->ifile->read(reinterpret_cast<char*>(&_data.clock), sizeof(unsigned int));
                            priv->ifile->read(reinterpret_cast<char*>(&_data.time), sizeof(unsigned int));
                            priv->ifile->read(reinterpret_cast<char*>(&_data.temp), sizeof(unsigned short));
                            for(int ii = 0; ii < 2; ii++) {
                                priv->ifile->read(reinterpret_cast<char*>(_header[ii]), 16 * sizeof(unsigned short));
                                priv->ifile->read(reinterpret_cast<char*>(&_data.data[ii * 128]),
                                                  128 * sizeof(unsigned short));
                            }
                        }
                    }

                    break;
                case CheckPoint:
                    priv->ifile->read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                    block_data = new char[size];
                    priv->ifile->read(block_data, size);
                    check_point(size, block_data);
                    break;
                case EndOfRun:
                    priv->ifile->read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
                    block_data = new char[size];
                    priv->ifile->read(block_data, size);
                    end_of_run(size, block_data);
                    break;
                default:
                    std::cout << "Unknown block data type: " << std::hex << header << " - " << code << std::dec << std::endl;
                }
                if(block_data) {
                    delete[] block_data;
                    block_data = nullptr;
                }

            } while(code != DataBlock && !(priv->ifile->bad() || priv->ifile->eof()));
        } else {
            priv->ifile->read(reinterpret_cast<char*>(&_data), sizeof(EventData));
            for(int ii = 0; ii < 2; ii++)
                memset(_header[ii], 0, 16 * sizeof(unsigned short));
        }

        if(priv->ifile->eof()) {
            return -1;
        } else if(priv->ifile->bad()) {
            return 0;
        } else {
            return 1;
        }
    } else
        return 2;
}

double AsciiRoot::time() const {
    short fpart = static_cast<short>(_data.time & 0xffff);
    short ipart = static_cast<short>((_data.time & 0xffff0000) >> 16);
    if(ipart < 0)
        fpart = static_cast<short>(fpart * -1);
    double tt = 100.0 * (ipart + (static_cast<unsigned short>(fpart) / 65535.));
    return tt;
}

double AsciiRoot::temp() const {
    return 0.12 * _data.temp - 39.8;
}
