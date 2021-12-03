/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_U_MOMENTUM_H
#define SOL_U_MOMENTUM_H

using namespace std;

double AU(double F, double D)
{
    return fmax(0, pow((1 - 0.1 * abs(F/D)),5));
}

void u_momentum_new(int imax,int jmax,double dx,double dy,double rho, double mu, vector<vector<double>> u, vector<vector<double>> &u_star, vector<vector<double>> &d_u, vector<vector<double>> v, vector<vector<double>> p, double velocity, double alpha,double dt,vector<vector<double>> u_o)
{
    int i,j;
    double Fe,Fw,Fn,Fs,aE,aW,aN,aS,aP,ap_o,pressure_term;
    double De,Dw,Dn,Ds;

    De  = mu*dy / dx;  //convective coefficients
    Dw  = mu*dy / dx;
    Dn  = mu*dx / dy;
    Ds  = mu*dx / dy;

    for(i=1;i<imax;i++)
    {
        for(j=1;j<jmax-1;j++)
        {
            Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
            Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
            Fn  = 0.5*rho*dx*(v[i][j+1]+v[i-1][j+1]);
            Fs  = 0.5*rho*dx*(v[i][j]+v[i-1][j]);

            aE = De * AU(Fe,De) + fmax(-Fe,0);
            aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
            aN = Dn * AU(Fn,Dn) + fmax(-Fn,0);
            aS = Ds * AU(Fs,Ds) + fmax(Fs,0);
            aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
            ap_o=1*dx*dy/dt;

            pressure_term = (p[i-1][j] - p[i][j]) * dy;

            u_star[i][j] = (alpha/aP) * ( (aE*u[i+1][j] + aW*u[i-1][j] + aN*u[i][j+1] + aS*u[i][j-1]) + ap_o*u_o[i][j] + pressure_term ) + (1-alpha)*u[i][j];

            d_u[i][j] = alpha * dy / aP;   //refer to Versteeg CFD book
        }
    }

    //set d_u for top and bottom BCs
    //they will be later used by the pressure correction equation
    //they should not be zero, or BCs of pressure correction will get messed up

    j = 0; //bottom
    for(i=1;i<imax;i++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
        Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
        Fn  = 0.5*rho*dx*(v[i][j+1]+v[i-1][j+1]);
        Fs  = 0;

        aE = De * AU(Fe,De) + fmax(-Fe,0);
        aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
        aN = Dn * AU(Fn,Dn) + fmax(-Fn,0);
        aS = 0.0;
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_u[i][j] = alpha * dy / aP;
    }

    j = jmax-1; //top
    for(i=1;i<imax;i++)
    {
        Fe  = 0.5*rho*dy*(u[i+1][j]+u[i][j]);
        Fw  = 0.5*rho*dy*(u[i-1][j]+u[i][j]);
        Fn  = 0;
        Fs  = 0.5*rho*dx*(v[i][j]+v[i-1][j]);

        aE = De * AU(Fe,De) + fmax(-Fe,0);
        aW = Dw * AU(Fw,Dw) + fmax(Fw,0);
        aN = 0.0;
        aS = Ds * AU(Fs,Ds) + fmax(Fs,0);
        aP = aE + aW + aN + aS + (Fe-Fw) + (Fn-Fs);
        d_u[i][j] = alpha * dy / aP;
    }

    //Apply BCs
    for(j=0;j<jmax;j++)
    {
        u_star[0][j] = -u_star[1][j];           //left wall
        u_star[imax][j] = -u_star[imax-1][j];   //right wall
        //u_star[0][j] = 0.0;           //left wall
        //u_star[imax][j] = 0.0;   //right wall
    }
    for(i=0;i<imax+1;i++)
    {
        u_star[i][0] = 0.0;                      //bottom wall
        u_star[i][jmax-1] = velocity;            //top wall
        //u_star[i][0] =  -u_star[i][1];                      //bottom wall
        //u_star[i][jmax-1] = -u_star[i][jmax-2];            //top wall
    }

}

#endif
