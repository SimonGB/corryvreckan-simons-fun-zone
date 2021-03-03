#ifndef __HIT_H__
#define __HIT_H__

/**
 * A class representing a hit
 */

#include <vector>

class Hit
{
    private:
        int _center;
        int _left;
        int _right;
        double _sig;
        
        void cpy(const Hit &h);
    public:
        Hit(int c=0, int l=0, int r=0, double s=0);
        Hit(const Hit &h);
        ~Hit();
        
        Hit &operator=(const Hit &h);
        
        int center() const { return _center; }
        int left() const { return _left; }
        int right() const { return _right; }
        double signal() const { return _sig; }
        int width() const { return _right - _left + 1; }
};
typedef std::vector< Hit> HitList;


#endif /*__HIT_H__*/
