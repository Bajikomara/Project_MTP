/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_V_MOMENTUM_H
#define SOL_V_MOMENTUM_H

using namespace std;

double AV(double F, double D)
{
    return fmax(0, pow((1 - 0.1 * abs(F/D)),5) );
}

void v_momentum_new(int imax,int jmax,double dx,double dy,double rho, double mu, vector<vector<double>> u, vector<vector<double>> v, vector<vector<double>> &v_star, vector<vector<double>> &d_v, vector<vector<double>> p, vector<vector<double>> T, double alpha,double dt,vector<vector<double>> v_o)
{
    int i,j;
    double Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,ap_o,pressure_term;
    double Boussinesq_source = 0.0;
    double De  = mu*dy / dx;  //convective coefficients
    double Dw  = mu*dy / dx;
    double Dn  = mu*dx / dy;
    double Ds  = mu*dx / dy;

    for(i=1;i<imax-1;i++)
    {
        for(j=1;j<jmax;j++)
        {
            Fe  = 0.5*rho*dy*(u[i+1][j]+u[i+1][j-1]);
            Fw  = 0.5*rho*dy*(u[i][j]+u[i][j-1]);
            Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
            Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

            aE = De * AV(Fe,De) + fmax(-Fe,0);
            aW = Dw * AV(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            ap_o=1*dx*dy/dt;
            pressure_term = (p[i][j-1] - p[i][j]) * dx;
            Boussinesq_source = 0.5*(T[i][j-1] + T[i][j]) * dx*dy ;

            v_star[i][j] = alpha/aP * ( (aE*v[i+1][j] + aW*v[i-1][j] + aN*v[i][j+1] + aS*v[i][j-1]) +ap_o*v_o[i][j]  + pressure_term + Boussinesq_source) + (1-alpha)*v[i][j];

            d_v[i][j] = alpha * dx / aP;   //refer to Versteeg CFD book
        }
    }

    //set d_v for left and right BCs
    //they will be later used by the pressure correction equation
    //they should not be zero, or BCs of pressure correction will get messed up
    //Apply BCs

    i = 0; //Left
    for(j=1;j<jmax;j++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i+1][j-1]);
        Fw  = 0.0;
        Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
        Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

        aE = De * AV(Fe,De) + fmax(-Fe,0);
        aW = 0.0;
        aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
        aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_v[i][j] = alpha * dx / aP;
    }

    i = imax-1; //top
    for(j=1;j<imax;j++)
    {
        Fe  = 0.0;
        Fw  = 0.5*rho*dy*(u[i][j]+u[i][j-1]);
        Fn  = 0.5*rho*dx*(v[i][j]+v[i][j+1]);
        Fs  = 0.5*rho*dx*(v[i][j-1]+v[i][j]);

        aE = 0.0;
        aW = Dw * AV(Fw,Dw) + fmax(Fw,0);
        aN = Dn * AV(Fn,Dn) + fmax(-Fn,0);
        aS = Ds * AV(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_v[i][j] = alpha * dx / aP;
    }

    //Apply BCs
    for(j=0;j<jmax+1;j++)
    {
        v_star[0][j] = 0.0;      //left wall
        v_star[imax-1][j] = 0.0;  //right wall
        //v_star[0][j] =  -v_star[1][j];      //left wall
        //v_star[imax-1][j] =  -v_star[imax-2][j];  //right wall
    }
    for(i=0;i<imax;i++)
    {
        v_star[i][0] =  -v_star[i][1];                     //bottom wall
        v_star[i][jmax] =  -v_star[i][jmax-1] ;          //top wall
        //v_star[i][0] =  0.0;                     //bottom wall
        //v_star[i][jmax] =  0.0 ;          //top wall
    }

}

#endif
