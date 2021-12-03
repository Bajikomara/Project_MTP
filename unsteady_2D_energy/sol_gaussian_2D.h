/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_GAUSSIAN_2D_H
#define SOL_GAUSSIAN_2D_H

using namespace std;

// Function to create Gaussian filter
vector<vector<double>> gaussian_distribution_2D(int imax_original, int jmax_original,int imax_local, int jmax_local, vector<double> x, vector<double> y, double x_center, double y_center)
{
    // Initializing standard deviation to 1.0
    double sigma = 1.0/sqrt(2*M_PI);
    double r, s = 2.0 * sigma * sigma;
    vector<vector<double>> GKernel(imax_original, vector<double>(jmax_original,0));

    // sum is for normalization
    double sum = 0.0;

    // generating 5x5 kernel
    for (int i = (imax_original - imax_local)/2; i < (imax_original - imax_local)/2+imax_local ; i++)
    {
        for (int j = (jmax_original - jmax_local)/2; j < (jmax_original - jmax_local)/2 + jmax_local ; j++)
        {
            r = sqrt((x[i] - x_center) * (x[i] - x_center) + (y[j] - y_center) * (y[j] - y_center));
            GKernel[i][j] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[i][j];
        }
    }

    // Normalizing the Kernel
    for (int i = (imax_original - imax_local)/2; i < (imax_original - imax_local)/2 + imax_local ; i++)
        for (int j = (jmax_original - jmax_local)/2; j < (jmax_original - jmax_local)/2 + jmax_local ; j++)
            GKernel[i][j] /= sum;

    return GKernel;
}

/*
// Function to create Gaussian filter
vector<vector<double>> gaussian_2D(int imax_original, int jmax_original,int imax_local, int jmax_local, vector<double> x, vector<double> y, double x_center, double y_center)
{
    // Initializing standard deviation to 1.0
    vector<vector<double>> kernal;
    double A, sigx, sigy;
    A = 1.0;
    sigx = 1.0;
    sigy = 1.0;

    for(int i=0; i<imax_local; i++)
    {
        for(int j=0; j<jmax_local; j++)
        {
            kernal[i][j] = exp(-(((x[i]-x_center)*(x[i]-x_center)/(2*sigx*sigx)) + ((y[i]-y_center)*(y[i]-y_center)/(2*sigy*sigy))));
        }
    }


}
*/

#endif

