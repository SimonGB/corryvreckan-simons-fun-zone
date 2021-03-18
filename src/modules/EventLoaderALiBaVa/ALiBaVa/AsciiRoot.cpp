#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cerrno>
#include <climits>

#include "AsciiRoot.h"
#include "Tracer.h"
#include "Hit.h"
#include "ChanList.h"
#include "utils.h"
#include "core/utils/log.h"

struct AsciiRootPriv
{
	std::ifstream *ifile;

	AsciiRootPriv() : ifile(0) {}
};

// decodes the header and returns a vector with the integers found
std::vector<int> decode_header(const std::string &h, AsciiRoot::XtraValues &xtra)
{
    std::vector<int> vout;
    std::istringstream istr(h);
    char *endptr;
    char buf[256];
    long val;

    xtra.clear();

    while (istr)
    {
        istr.getline(buf, sizeof(buf), ';');
        if (!istr)
            break;

        errno = 0;
        val = strtol(buf, &endptr, 0);

        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0))
        {
            std::string sval(buf), sout;
            sout = trim_str(sval);
            if (!sout.empty())
                xtra.push_back( sout );
        }
        else if ( endptr == buf || *endptr != '\0' )
        {
            std::string sval(buf), sout;
            sout = trim_str(sval);
            if (!sout.empty())
                xtra.push_back( sout );
        }
        else
        {
            vout.push_back(atoi(buf) );
        }
    }
    return vout;
}

AsciiRoot::AsciiRoot(const char *nam, const char *pedfile, const char *gainfile)
    : DataFileRoot(nam, pedfile, gainfile), priv(0)
{
    priv = new AsciiRootPriv;

    if (nam)
        open(nam);
}
AsciiRoot::~AsciiRoot()
{
    if (priv->ifile)
    {
        priv->ifile->close();
    }
    delete priv;
}

bool AsciiRoot::valid() const
{
	LOG(DEBUG) << "Testing if file is Valid";
	return (priv->ifile!=0);
	LOG(DEBUG) << "Tested if file is Valid finished";
}


void AsciiRoot::open(const char *name)
{
    priv->ifile = new std::ifstream(name);
    if (!(*priv->ifile))
    {
        std::cout << "Could not open data file: " << name << std::endl;
        delete priv->ifile;
        priv->ifile = 0;
		return;
    }
    std::string header;
    unsigned int ic, lheader;
    char c;
    priv->ifile->read((char *)&_t0, sizeof(time_t));
    //std::cout << "64b: " << ctime(&_t0) << std::endl;
    priv->ifile->read((char *)&_type, sizeof(int));
    //std::cout << "type_ " << _type << std::endl;
    if ( _type > LastRType )
    {
        priv->ifile->seekg(0, std::ios::beg);
        priv->ifile->read((char *)&_t0, sizeof(int));
        priv->ifile->read((char *)&_type, sizeof(int));
        //std::cout << "32b: " << ctime(&_t0) << std::endl;
    }


    priv->ifile->read((char *)&lheader, sizeof(unsigned int));
    for (ic=0; ic<80; ic++)
    {
        priv->ifile->read(&c, sizeof(char));
        header.append(1, c);
    }
    header = trim_str(header);

    if (header[0]!='V' && header[0]!='v')
    {
        _version = 0;
    }
    else
    {
        _version = int(header[1]-'0');
        header = header.substr(5);
    }

    std::cout << "type: " << _type << " header: " << header << std::endl;
    std::vector<int> param = decode_header(header, _xtra);
    priv->ifile->read((char *)_ped, max_nchan*sizeof(double));
    priv->ifile->read((char *)_noise, max_nchan*sizeof(double));
    switch (_type)
    {
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
            if (param.empty())
                _nevts = 100000;
            else
                _nevts = param[0];
            _npoints = _from = _to = _step = 0;
            break;
        default:
        	break;
    }
    _scantype = UnknownScan;
    if ( type() == Calibration )
    {
        if ( nxtra() )
            _scantype = Time;
        else
            _scantype = Charge;
    }
    else if ( type() == LaserSync )
        _scantype = LaserScan;

    data_start = priv->ifile->tellg();
    saved_state = data_start;
}

void AsciiRoot::save()
{
    if (priv->ifile)
        saved_state = priv->ifile->tellg();
}

void AsciiRoot::restore()
{
    if (priv->ifile)
    {
        priv->ifile->seekg(saved_state, std::ios::beg);
    }
}

void AsciiRoot::rewind()
{
    if (priv->ifile)
    {
        priv->ifile->clear();
        priv->ifile->seekg(data_start, std::ios::beg);
    }
}

void AsciiRoot::close()
{
    if (priv->ifile)
    {
        priv->ifile->close();
        delete priv->ifile;
        priv->ifile = 0;
    }
}

void AsciiRoot::get_scan_values(short &delay, short &charge)
{
    delay = int(value()) >> 16;
    charge = int(value()) & 0xffff;
}

int AsciiRoot::read_event()
{
    if (priv->ifile)
    {
        unsigned int header, size, user=0, code=0;
        char *block_data=0;
        if (_version)
        {
            do
            {
                do
                {
                    priv->ifile->read((char *)&header, sizeof(unsigned int));
                    if (priv->ifile->bad() || priv->ifile->eof())
                        return -1;

                    code = (header>>16) & 0xFFFF;
                } while ( code != 0xcafe );

                code = header & 0x0fff;
                user = header & 0x1000;
                switch (code)
                {
                    case NewFile:
                        priv->ifile->read((char *)&size, sizeof(unsigned int));
                        block_data = new char[size];
                        priv->ifile->read(block_data, size);
                        new_file(size, block_data);
                        break;
                    case StartOfRun:
                        priv->ifile->read((char *)&size, sizeof(unsigned int));
                        block_data = new char[size];
                        priv->ifile->read(block_data, size);
                        start_of_run(size, block_data);
                        break;
                    case DataBlock:
                        priv->ifile->read((char *)&size, sizeof(unsigned int));
                        if (user)
                        {
                            reset_data();
                            block_data = new char[size];
                            priv->ifile->read(block_data, size);
                            new_data_block(size, block_data);
                        }
                        else
                        {
                            if ( _version == 1 )
                            {
                                priv->ifile->read((char *)&_data, sizeof(EventData));
                                for (int ii=0; ii<2; ii++)
                                    memset(_header[ii], 0, 16*sizeof(unsigned short));

                            }
                            else
                            {
                                priv->ifile->read((char *)&_data.value, sizeof(double));
                                if (_version > 2)
                                    priv->ifile->read((char *)&_data.clock, sizeof(unsigned int));
                                priv->ifile->read((char *)&_data.time, sizeof(unsigned int));
                                priv->ifile->read((char *)&_data.temp, sizeof(unsigned short));
                                for (int ii=0; ii<2; ii++)
                                {
                                    priv->ifile->read((char *)_header[ii], 16*sizeof(unsigned short));
                                    priv->ifile->read((char *)&_data.data[ii*128], 128*sizeof(unsigned short));
                                }
                            }
                        }

                        break;
                    case CheckPoint:
                        priv->ifile->read((char *)&size, sizeof(unsigned int));
                        block_data = new char[size];
                        priv->ifile->read(block_data, size);
                        check_point(size, block_data);
                        break;
                    case EndOfRun:
                        priv->ifile->read((char *)&size, sizeof(unsigned int));
                        block_data = new char[size];
                        priv->ifile->read(block_data, size);
                        end_of_run(size, block_data);
                        break;
                    default:
                        std::cout << "Unknown block data type: " << std::hex << header << " - " << code << std::dec << std::endl;
                }
                if (block_data)
                {
                    delete [] block_data;
                    block_data = 0;
                }

            } while ( code != DataBlock && !(priv->ifile->bad() || priv->ifile->eof()) );
        }
        else
        {
            priv->ifile->read((char *)&_data, sizeof(EventData));
            for (int ii=0; ii<2; ii++)
                memset(_header[ii], 0, 16*sizeof(unsigned short));
        }

        if (priv->ifile->eof())
        {
            std::cout << "End of file" << std::endl;
            return -1;
        }
        else if (priv->ifile->bad())
        {
            std::cout << "Problems with data file" << std::endl;
            return -1;
        }
        else
        {
            //process_event();
            return 0;
        }
    }
    else
        return -1;
}

double AsciiRoot::time() const
{
    unsigned short fpart = _data.time & 0xffff;
    short ipart = (_data.time & 0xffff0000)>>16;
    if (ipart<0)
        fpart *= -1;
    //double tt = 100.*(1. -(ipart + (fpart/65535.)));
    double tt = 100.0*(ipart + (fpart/65535.));
    return tt;
}

double AsciiRoot::temp() const
{
    return 0.12*_data.temp - 39.8;
}
