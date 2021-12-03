/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_PDMA_H_
#define SOL_PDMA_H_

using namespace std;

vector<double> PDMA(vector<vector<double>> &A, vector<double> &b, vector<vector<double>> &diagonals)
{
    int N = A.size();
    vector<double> x(A.size(),0);
    vector<double> alpha(A.size(),0);
    vector<double> gam(A.size()-1,0);
    vector<double> delta(A.size()-2,0);
    vector<double> bet(A.size(),0);

    vector<double> c(A.size(),0);
    vector<double> z(A.size(),0);

    vector<double> d = diagonals[0];
    vector<double> e = diagonals[1];
    vector<double> f = diagonals[2];
    vector<double> h = diagonals[3];
    vector<double> g = diagonals[4];

    // Factor A=LR
    alpha[0]=d[0];
    gam[0]=e[0]/alpha[0];
    delta[0]=f[0]/alpha[0];
    bet[1] = h[1];
    alpha[1]=d[1]-bet[1]*gam[0];
    gam[1]=( e[1]-bet[1]*delta[0] )/alpha[1];
    delta[1]=f[1]/alpha[1];

    for(int k=2; k<A.size()-2; k++)
    {
        bet[k]=h[k]-g[k]*gam[k-2];
        alpha[k]=d[k]-g[k]*delta[k-2]-bet[k]*gam[k-1];
        gam[k]=( e[k]-bet[k]*delta[k-1] )/alpha[k];
        delta[k]=f[k]/alpha[k];
    }

    bet[N-2]=h[N-2]-g[N-2]*gam[N-4];
    alpha[N-2]=d[N-2]-g[N-2]*delta[N-4]-bet[N-2]*gam[N-3];
    gam[N-2]=( e[N-2]-bet[N-2]*delta[N-3] )/alpha[N-2];
    bet[N-1]=h[N-1]-g[N-1]*gam[N-3];
    alpha[N-1]=d[N-1]-g[N-1]*delta[N-3]-bet[N-1]*gam[N-2];

    // Update b=Lc
    c[0]=b[0]/alpha[0];
    c[1]=(b[1]-bet[1]*c[0])/alpha[1];

    for(int k=2; k<A.size(); k++)
    {
        c[k]=( b[k]-g[k]*c[k-2]-bet[k]*c[k-1] )/alpha[k];
    }

    // Back substitution Rx=c
    x[N-1]=c[N-1];
    x[N-2]=c[N-2]-gam[N-2]*x[N-1];

    for(int k=A.size()-3;k>=0;k--)
    {
       x[k]=c[k]-gam[k]*x[k+1]-delta[k]*x[k+2];
    }

    return x;
}

#endif
