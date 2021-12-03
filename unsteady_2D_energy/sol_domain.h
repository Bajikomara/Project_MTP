/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_DOMAIN_H_
#define SOL_DOMAIN_H_

using namespace std;

class Domain
{
public:
    double dx;
    double dy;
    vector <double> x;
    vector <double> y;

};

void create_domain(Domain &D)
{
    D.dx = (x_end - x_start) / (imax-1);
    D.dy = (y_end - y_start) / (jmax-1);
    D.x = linspace_interval(x_start, x_end, D.dx, imax);
    D.y = linspace_interval(y_start, y_end, D.dy, jmax);
}

#endif

