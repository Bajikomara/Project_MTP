/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_GET_DIAGONALS_H
#define SOL_GET_DIAGONALS_H

using namespace std;

void get_diagonals(int imax, int jmax, vector<vector<double>> &Ap, vector<vector <double>> &diagonals)
{
    int N = imax*jmax;
    int i,j;
    //vector<vector <double>> diagonals(5, vector<double> (imax*jmax));
    //vector<vector <double>> diagonals;
    for(i=0; i<Ap.size(); i++)
    {
        for(j=0; j<Ap[0].size(); j++)
        {
            if(i==j)  //main diagonal
            {
                //diagonals[0].push_back(Ap[i][j]);
                diagonals[0][i] = (Ap[i][j]);
            }
            else if(i == j-1) // upper diagonal
            {
                //diagonals[1].push_back(Ap[i][j]);
                diagonals[1][i] = (Ap[i][j]);
            }
            else if(i == j-2)  //upper upper diagonal
            {
                //diagonals[2].push_back(Ap[i][j]);
                diagonals[2][i] = (Ap[i][j]);
            }
            else if(j == i-1)  // sub diagonal
            {
                //diagonals[3].push_back(Ap[i][j]);
                diagonals[3][i] = (Ap[i][j]);
            }
            else if(j == i-2)  // sub sub diagonal
            {
                //diagonals[4].push_back(Ap[i][j]);
                diagonals[4][i] = (Ap[i][j]);
            }
        }
    }
}

#endif
