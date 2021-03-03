#include "Hit.h"

Hit::Hit(int c, int l, int r, double s) :
    _center(c), _left(l), _right(r), _sig(s)
{
}

Hit::Hit(const Hit &h)
{
    cpy(h);
}
Hit::~Hit()
{
}

void Hit::cpy(const Hit &h)
{
    _center = h._center;
    _left = h._left;
    _right = h._right;
    _sig = h._sig;
}
Hit &Hit::operator=(const Hit &h)
{
    if (&h!=this)
        cpy(h);

    return *this;
}
