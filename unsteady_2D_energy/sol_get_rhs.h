/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_GET_RHS_H
#define SOL_GET_RHS_H

using namespace std;

void get_rhs(int imax, int jmax, double dx, double dy, double rho, vector<vector<double>> u_star, vector<vector<double>> v_star, vector<double> &rhsp)
{
    int i,j;
    int stride = jmax;
    int position;


    // RHS is the same for all nodes except the p_prime(1,1)
    // Because p(1,1) is set to be zero, it has no pressure correction

    for(j=0;j<jmax;j++)
    {
        for(i=0;i<imax;i++)
        {
            position = i+(j-0)*stride;
            rhsp[position] = rho * (u_star[i][j]*dy - u_star[i+1][j]*dy  + v_star[i][j]*dx - v_star[i][j+1]*dx)  ;
        }
    }

    // modify for p_prime(1,1)
    rhsp[0] = 0.0;
}

#endif
