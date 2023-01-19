#ifndef DATA_H_
#define DATA_H_

/**
 *  Data is received from the USB port with this format
 *
 */
struct EventBlock {
    unsigned int clock;
    unsigned int time;
    unsigned short temp;
    unsigned short data[256];
};

/**
 * We add value as an extra parameter when moving the
 * data within the application. It tells the value of
 * the variables which are scanned
 */
struct EventData : public EventBlock {
    double value;
};

struct EventDataBlock : public EventData {
    unsigned short header[32];
};

#endif /*DATA_H_*/
