/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_PRESSURE_CORRECT_H
#define SOL_PRESSURE_CORRECT_H

#include "sol_PDMA.h"
#include "sol_get_diagonals.h"

using namespace std;
vector<vector<double>> pressure_correct(int imax, int jmax, vector<double> &rhsp, vector<vector<double>> &Ap, vector<vector<double>> p, vector<vector<double>> &p_prime, double alpha)
{
    int i,j,z = 0;
    vector<vector<double>> pressure = p;
    vector <double> p_prime_interior;
    vector<vector <double>> diagonals(5, vector<double> (imax*jmax));
    get_diagonals(imax,jmax,Ap,diagonals);

    p_prime_interior = PDMA(Ap,rhsp,diagonals);
    //convert pressure correction in to a matrix
    //update pressure values

    for(j=0;j<jmax;j++)
    {
        for(i=0;i<imax;i++)
        {
            p_prime[i][j] = p_prime_interior[z];
            z = z + 1;
            pressure[i][j] = p[i][j] + alpha*p_prime[i][j];
        }
    }
    pressure[0][0] = 0;

    return pressure;

}

#endif
