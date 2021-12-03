/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_READINPUT_H_
#define SOL_READINPUT_H_

using namespace std;
void read_geometry_file()
{
    ifstream ifs("input_geometry.txt");

    if (!ifs){
        cout << "Input file not found !!" << endl;
        exit(0);
    }

    string param;
    string value;

    while(!ifs.eof())
    {
        ifs >> param >> value;
        if (param.compare("x_start")==0)
        {
            x_start = atof(value.c_str());
        }
        else if (param.compare("x_end")==0)
        {
            x_end = atof(value.c_str());
        }
        else if (param.compare("y_start")==0)
        {
            y_start = atof(value.c_str());
        }
        else if (param.compare("y_end")==0)
        {
            y_end = atof(value.c_str());
        }
        else if (param.compare("imax")==0)
        {
            imax = atof(value.c_str());
        }
        else if (param.compare("jmax")==0)
        {
            jmax = atof(value.c_str());
        }
        else
        {
            cout << "Invalid parameter " << param << " check input file " << endl;
            exit(0);
        }
    }
}

void read_simulation_file()
{
    ifstream ifs("input_simulation.txt");

    if (!ifs){
        cout << "Input file not found !!" << endl;
        exit(0);
    }

    string param;
    string value;

    while(!ifs.eof())
    {
        ifs >> param >> value;
        if(param.compare("simulation_file_name")==0)
        {
            mesh = value;
        }
        else if (param.compare("max_iterations")==0)
        {
            max_iteration = atof(value.c_str());
        }
        else if (param.compare("prandtl_number")==0)
        {
            Pr = atof(value.c_str());
        }
        else if (param.compare("rayleigh_number")==0)
        {
            Ra = atof(value.c_str());
        }
        else if (param.compare("density")==0)
        {
            rho = atof(value.c_str());
        }
        else if (param.compare("velocity")==0)
        {
            velocity = atof(value.c_str());
        }
        else if (param.compare("T_top")==0)
        {
            T_top = atof(value.c_str());
        }
        else if (param.compare("T_bottom")==0)
        {
            T_bottom = atof(value.c_str());
        }
        else if (param.compare("under_relaxation_P")==0)
        {
            alphaP = atof(value.c_str());
        }
        else if (param.compare("under_relaxation_U")==0)
        {
            alphaU = atof(value.c_str());
        }
        else if (param.compare("under_relaxation_V")==0)
        {
            alphaV = atof(value.c_str());
        }
        else if (param.compare("under_relaxation_T")==0)
        {
            alphaT = atof(value.c_str());
        }
        else if (param.compare("convergence_criteria")==0)
        {
            tol = atof(value.c_str());
        }
        else
        {
            cout << "Invalid parameter " << param << " check input file " << endl;
            exit(0);
        }
    }
}



#endif
