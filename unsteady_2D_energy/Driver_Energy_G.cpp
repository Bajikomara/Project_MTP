#include <bits/stdc++.h>

using namespace std;

double x_start, x_end, y_start, y_end;
int imax, jmax;
int max_iteration;
double mu, rho, velocity, alphaU, alphaV, alphaT, alphaP, tol;
double Pr, Ra, Di, T_top, T_bottom;
double maxRes = 1000;
int iteration ;
string mesh;

#include "sol_vectorlib.h"
#include "sol_readinput.h"
#include "sol_domain.h"
#include "sol_gaussian_2D.h"
#include "sol_u_momentum.h"
#include "sol_v_momentum.h"
#include "sol_energy_equation.h"
#include "sol_get_rhs.h"
#include "sol_coefficient_matrix.h"
#include "sol_pressure_correct.h"
#include "sol_update_velocity.h"
#include "sol_check_divergence.h"

int main()
{
    double Re;
    int i,j;
    Domain D;
    read_geometry_file();
    read_simulation_file();
    create_domain(D);

    mu = sqrt(Pr/Ra);
    Di = 1 / sqrt(Ra*Pr); 
    Re = (rho * velocity * (x_end - x_start))/mu;
    string filename;
	

    //Variable declaration
    vector<vector<double>> p(imax, vector<double> (jmax,0));                // p = Pressure
    vector<vector<double>> p_star(imax, vector<double> (jmax,0));           // Guessed Pressure
    vector<vector<double>> p_prime(imax, vector<double> (jmax,0));          // pressure correction
    vector<vector<double>> divergence(imax, vector<double> (jmax,0));
    vector<double> rhsp(imax + (jmax-1)*jmax, 0);                           // Right hand side vector of pressure correction equation
    vector<vector<double>> Ap(imax*jmax, vector<double> (imax*jmax,0));
    vector<vector <double>> diagonals(5, vector<double> (imax*jmax));
    vector<vector <double>> div(imax, vector<double> (jmax,0));
    vector<vector<double>> p_o(imax, vector<double> (jmax,0));

    vector < vector <double>> T_int(imax,vector<double>(jmax,0));
    vector < vector <double>> T(imax,vector<double>(jmax,0));
    vector < vector <double>> Told(imax,vector<double>(jmax,0));
    vector < vector <double>> TRes(imax,vector<double>(jmax,0));
    vector<vector<double>> T_o(imax, vector<double> (jmax,0));

    // Vertical velocity
    vector<vector<double>> v_star(imax, vector<double> (jmax+1,0));
    vector<vector<double>> vold (imax, vector<double> (jmax+1,0));
    vector<vector<double>> vRes(imax, vector<double> (jmax+1,0));
    vector<vector<double>> v (imax, vector<double> (jmax+1,0));
    vector<vector<double>> d_v(imax, vector<double> (jmax+1,0));    //velocity correction coefficient
    vector<vector<double>> v_o(imax, vector<double> (jmax,0));

    // Horizontal Velocity
    vector<vector<double>> u_star(imax+1, vector<double> (jmax,0));
    vector<vector<double>> uold(imax+1, vector<double> (jmax,0));
    vector<vector<double>> uRes(imax+1, vector<double> (jmax,0));
    vector<vector<double>> u(imax+1, vector<double> (jmax,0));
    vector<vector<double>> d_u(imax+1, vector<double> (jmax,0));  //velocity correction coefficient
    vector<vector<double>> u_o(imax, vector<double> (jmax,0));

    

    double laser_center_x = (x_end + x_start)/2;
    double laser_center_y = (y_end + y_start)/2;
    int d;
    d = 2;
    double radius = d*D.dx;
    vector<double> x_laser_vec = linspace_interval(laser_center_x - radius, laser_center_x + radius , D.dx, d+2);
    vector<double> y_laser_vec = linspace_interval(laser_center_y - radius , laser_center_y + radius , D.dy, d+2);


    for(i=0;i<u_star.size();i++)
    {
        u_star[i][jmax-1] = velocity;
        u[i][jmax-1] = velocity;
    }

    T_int = gaussian_distribution_2D(imax, jmax, x_laser_vec.size(), y_laser_vec.size(), D.x, D.y, laser_center_x, laser_center_y);
    T = T_int;
    Told = T_int;

    /*
    for(i=0;i<T.size();i++)
    {
        T[i][jmax-1] = T_top;     //top temperature
        Told[i][jmax-1] = T_top;
        T[i][0] = T_bottom;       //bottom temperature
        Told[i][0] = T_bottom;
    }
    */
   double dt=2,t=0;
   int cnt=0;
    while( t < 50.0){
        cnt++;
        cout<<"Time step = "<<t<<" count = "<<cnt<<endl;
        t+=dt;
        
        u_o=u;
        v_o=v;
        T_o=T;
        p_star = p_o;
        double rho_o= rho;
        iteration=0;
        maxRes=0.5;
        // ---------- iterations -------------------//
        while ( (iteration <= max_iteration) && (maxRes >= tol) )
        {
            ofstream write;
            ostringstream filenamestream;
            filenamestream << mesh << "_" << imax << "x" << jmax << "_result_Ra_"  << Ra << "_time_step_" << t << ".dat";
            filename = filenamestream.str();
            //cout << filename << endl;

            if(iteration == 0 && t==2)
            {
                write.open(filename.c_str());
                write << fixed << setprecision(8);
                write << "VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\t\"vres\"\t\"T\"\t\"p\"\n";
                write << "ZONE T=" << " \" "  << iteration << " \" " << ", I=" << imax << ", J=" << jmax << ", F = point\n";
                for (i=0 ; i < imax ; i++)
                {
                    for (j=0 ; j < jmax ; j++)
                    {
                        write << D.x[i] << "\t" << D.y[j] << "\t" << u[i][j] << "\t" <<v[i][j] << "\t" <<sqrt(pow(u[i][j],2)+ pow(v[i][j],2)) << "\t" << T[i][j] << "\t" << p[i][j] <<endl ;
                        
                    }
                }
                write << "\n\n";
            }

            u_momentum_new(imax,jmax,D.dx,D.dy,rho,mu,u,u_star,d_u,v,p_star,velocity,alphaU,dt,u_o);
            v_momentum_new(imax,jmax,D.dx,D.dy,rho,mu,u,v,v_star,d_v,p_star,T,alphaV,dt,v_o);
            uold = u;
            vold = v;
            get_rhs(imax,jmax,D.dx,D.dy,rho,u_star,v_star,rhsp);
            coefficient_matrix(imax,jmax,D.dx,D.dy,rho,d_u,d_v,Ap);
            p = pressure_correct(imax,jmax,rhsp,Ap,p_star,p_prime,alphaP);
            update_velocity(imax,jmax,u,u_star,v,v_star,p_prime,d_u,d_v,velocity);
            check_divergence(imax,jmax,D.dx,D.dy,u,v,div);
            energy_equation(imax,jmax,D.dx,D.dy,u,v,T,Told,T_int,Di,alphaT,T_top,T_bottom,dt,T_o);
            p_star = p;
            p_o = p;
            
            //transform(T.begin(), T.end(), T_int.begin(), back_inserter(Told), plus<double>());
            //transform(T.begin(), T.end(), T_int.begin(), back_insert_iterator, plus<double>());
            /*
            for(i=0;i<T.size();i++)
            {
                for(j=0;j<T[0].size();j++)
                {
                    Told[i][j] = T[i][j] + T_int[i][j];
                }
            }
            */
            Told = T;


            //find maximum residual in the domain
            for(i=0;i<v.size();i++)
            {
                for(j=0;j<v[0].size();j++)
                {
                    vRes[i][j] = abs(v[i][j] - vold[i][j]);
                }
            }

            for(i=0;i<u.size();i++)
            {
                for(j=0;j<u[0].size();j++)
                {
                    uRes[i][j] = abs(u[i][j] - uold[i][j]);
                }
            }
            for(i=0;i<T.size();i++)
            {
                for(j=0;j<T[0].size();j++)
                {
                    TRes[i][j] = abs(T[i][j] - Told[i][j]);
                }
            }

            vector<double> oneDimU,oneDimV,oneDimT;
            for(int i = 0; i < imax; i++)
            {
                for(int j = 0; j < jmax; j++)
                {
                    oneDimU.push_back(uRes[i][j]);
                    oneDimV.push_back(vRes[i][j]);
                    oneDimT.push_back(TRes[i][j]);
                }
            }
            double maxRes_u = *max_element(oneDimU.begin(),oneDimU.end());
            double maxRes_v = *max_element(oneDimV.begin(),oneDimV.end());
            double maxRes_T = *max_element(oneDimT.begin(),oneDimT.end());
            vector<double> vec = {maxRes_u, maxRes_v, maxRes_T};
            maxRes = *max_element(vec.begin(),vec.end());

            

            cout << "It: "<< iteration << " " <<"Res: "<< maxRes <<endl;

            if (maxRes > 2)
            {
                cout << "Solution is not converging!" << endl;
                break;
            }
            iteration = iteration + 1;

        }
    
        ofstream write;
        write.open(filename.c_str());
        write << fixed << setprecision(8);
        write << "VARIABLES = \"x\"\t\"y\"\t\"u\"\t\"v\"\t\"vres\"\t\"T\"\t\"p\"\n";
        write << "ZONE T=" << " \" "  << iteration << " \" " << ", I=" << imax << ", J=" << jmax << ", F = point\n";
        for (i=0 ; i < imax ; i++)
        {
            for (j=0 ; j < jmax ; j++)
            {
                write << D.x[i] << "\t" << D.y[j] << "\t" << u[i][j] << "\t" <<v[i][j] << "\t" <<sqrt(pow(u[i][j],2)+ pow(v[i][j],2)) << "\t" << T[i][j] << "\t" << p[i][j] <<endl ;
                
            }
        }
        if(maxRes>2){
            cout<<"time_step"<<t<<endl;
            break;
        }
    }
    return 0;
}
