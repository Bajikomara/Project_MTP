/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_COEFFICIENT_MATRIX_H
#define SOL_COEFFICIENT_MATRIX_H

using namespace std;

void coefficient_matrix(int imax, int jmax, double dx, double dy, double rho, vector<vector<double>> d_u, vector<vector<double>> d_v, vector<vector<double>> &Ap)
{
    int i,j;
    int N = imax*jmax;
    int stride = jmax;
    int position;
    double aE,aW,aN,aS,aP;
    //P-prime for the boundary nodes is set to zero
    //the interior nodes imax-2*jmax-2 solved implicitly for P-prime

    for(j=0;j<jmax;j++)
    {
        for(i=0;i<imax;i++)
        {
            position = i + (j-0)*stride;
            aE = 0.0;
            aW = 0.0;
            aN = 0.0;
            aS = 0.0;

            //set BSc for four corners
            if(i == 0 && j == 0)
            {
                 Ap[position][position] = 1.0;                  //pressure correction at the first node is zero
            }

            else if(i == imax-1 && j == 0)
            {
                Ap[position][position-1] = -rho*d_u[i][j]*dy;
                aW = -Ap[position][position-1];
                Ap[position][position+stride] = -rho*d_v[i][j+1]*dx;
                aN = -Ap[position][position+stride];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            else if(i == 0 && j == jmax-1)
            {
                Ap[position][position+1] = -rho*d_u[i+1][j]*dy;
                aE = -Ap[position][position+1];
                Ap[position][position-stride] = -rho*d_v[i][j]*dx;
                aS = -Ap[position][position-stride];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            else if(i == imax-1 && j == jmax-1)
            {
                Ap[position][position-1] = -rho*d_u[i][j]*dy;
                aW = -Ap[position][position-1];
                Ap[position][position-stride] = -rho*d_v[i][j]*dx;
                aS = -Ap[position][position-stride];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            //set four boundaries
            else if (i == 0)
            {
                Ap[position][position+1] = -rho*d_u[i+1][j]*dy;
                aE = -Ap[position][position+1];
                Ap[position][position+stride] = -rho*d_v[i][j+1]*dx;
                aN = -Ap[position][position+stride];
                Ap[position][position-stride] = -rho*d_v[i][j]*dx;
                aS = -Ap[position][position-stride];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            else if (j == 0)
            {
                Ap[position][position+1] = -rho*d_u[i+1][j]*dy;
                aE = -Ap[position][position+1];
                Ap[position][position+stride] = -rho*d_v[i][j+1]*dx;
                aN = -Ap[position][position+stride];
                Ap[position][position-1] = -rho*d_u[i][j]*dy;
                aW = -Ap[position][position-1];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            else if (i == imax-1)
            {
                Ap[position][position+stride] = -rho*d_v[i][j+1]*dx;
                aN = -Ap[position][position+stride];
                Ap[position][position-stride] = -rho*d_v[i][j]*dx;
                aS = -Ap[position][position-stride];
                Ap[position][position-1] = -rho*d_u[i][j]*dy;
                aW = -Ap[position][position-1];
                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;
            }

            else if (j == jmax-1)
            {
                Ap[position][position+1] = -rho*d_u[i+1][j]*dy;
                aE = -Ap[position][position+1];
                Ap[position][position-stride] = -rho*d_v[i][j]*dx;
                aS = -Ap[position][position-stride];
                Ap[position][position-1]  = -rho*d_u[i][j]*dy;
                aW = -Ap[position][position-1] ;
                aP = aE + aN + aW + aS;
                Ap[position][position]  = aP;
            }

            else
            {
                Ap[position][position-1] = -rho*d_u[i][j]*dy;                          //sub diagonal
                aW = -Ap[position][position-1];

                Ap[position][position+1] = -rho*d_u[i+1][j]*dy;                        //upper diagonal
                aE = -Ap[position][position+1];

                Ap[position][position-stride] = -rho*d_v[i][j]*dx;                     //sub sub diagonal
                aS = -Ap[position][position-stride];

                Ap[position][position+stride] = -rho*d_v[i][j+1]*dx;                   //upper upper diagonal
                aN = -Ap[position][position+stride];

                aP = aE + aN + aW + aS;
                Ap[position][position] = aP;                                           //main diagonal
            }
        }
    }
}

#endif
