/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_ENERGY_EQUATION_H
#define SOL_ENERGY_EQUATION_H

using namespace std;

double AT(double F, double Diffusion)
{
    return fmax(0, pow((1 - 0.1 * abs(F/Diffusion)),5));
}

void energy_equation(int imax,int jmax,double dx,double dy, vector<vector<double>> u, vector<vector<double>> v, vector<vector<double>> &T, vector<vector<double>> Told, vector<vector<double>> T_int, double D, double alpha,double T_top, double T_bottom,double dt,vector<vector<double>> T_o)
{

    int i,j;
    double Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,ap_o,pressure_term;
    double De,Dw,Dn,Ds,U,V;

    for(i=0;i<T.size();i++)
        {
            for(j=0;j<T[0].size();j++)
            {
                Told[i][j] = T[i][j] + T_int[i][j];
            }
        }


    De  = D*dy / dx;  //convective coefficients
    Dw  = D*dy / dx;
    Dn  = D*dx / dy;
    Ds  = D*dx / dy;

    for(i=1;i<imax-1;i++)
    {
        for(j=1;j<jmax-1;j++)
        {
            Fe  = dy*(u[i+1][j]);
            Fw  = dy*(u[i][j]);
            Fn  = dx*(v[i][j+1]);
            Fs  = dx*(v[i][j]);

            U = 0.5*(u[i+1][j] + u[i][j]);
            V = 0.5*(v[i][j+1] + v[i][j]);

            aE = De * AT(Fe,De) + fmax(-Fe,0);
            aW = Dw * AT(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AT(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AT(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            ap_o=1*dx*dy/dt;
            T[i][j] = (alpha/aP) * ( (aE*Told[i+1][j] + aW*Told[i-1][j] + aN*Told[i][j+1] + aS*Told[i][j-1] +ap_o*T_o[i][j] )) + (1-alpha)*Told[i][j];

        }
    }

    //Apply BCs
    for(j=0;j<jmax;j++)
    {
        //T[0][j] = T[1][j];           //left wall
        //T[imax-1][j] = T[imax-2][j];   //right wall
        T[0][j] = 0.0;           //left wall
        T[imax-1][j] = 0.0;   //right wall

    }
    for(i=0;i<imax;i++)
    {
        T[i][0] = 0.0;               //bottom wall
        T[i][jmax-1] = 0.0;            //top wall
        //T[i][0] = T[i][1];               //bottom wall
        //T[i][jmax-1] = T[i][jmax-2];            //top wall
    }

}

#endif
