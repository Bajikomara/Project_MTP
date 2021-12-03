/**
    |````````````````````````````````|
    | 2D Incompressible              |
    | Sachin Gupta                   |
    | August 21                      |
    |________________________________|
**/

#ifndef SOL_VECTORLIB_H
#define SOL_VECTORLIB_H

using namespace std;

template<typename T>
vector <double> linspace(T start_in, T end_in, int num_in)
{

    vector<double> linspaced;

    double start_ = static_cast<double>(start_in);
    double end_ = static_cast<double>(end_in);
    double num_ = static_cast<double>(num_in);

    if (num_ == 0) { return linspaced; }
    if (num_ == 1)
    {
        linspaced.push_back(start_);
        return linspaced;
    }

    double delta = (end_ - start_) / (num_ - 1);

    for(int i=0; i < num_-1; ++i)
    {
        linspaced.push_back(start_ + delta * i);
    }

    linspaced.push_back(end_);
    return linspaced;
}


template<typename T>
vector <double> linspace_interval(T start_in, T end_in, double del_in, int num_in)
{

    vector <double> linspaced;

    double start_ = static_cast<double>(start_in);
    double end_ = static_cast<double>(end_in);
    double del_ = static_cast<double>(del_in);
    int num_ = static_cast<int> (num_in);

    if (del_ == 0) { return linspaced; }
    int i = 0;
    for(i=0;i<num_;i++)
    {
        linspaced.push_back(start_ + del_ * i);
    }

    return linspaced;
}

template <typename T>
void print_vector_1d(vector<T> vec)
{
    cout << "Shape: ( " <<vec.size()<< ",1)" << endl;
    for (T d : vec)
        cout << d << "\t";
    cout << endl;
}

template <typename T>
void print_vector_2d(vector<vector<T>> vec)
{
    cout << "Shape:(" <<vec.size()<< "," << vec[0].size() <<")" << endl;
    for(unsigned int i=0 ; i<vec.size(); i++)
    {
        for(unsigned int j=0; j<vec[0].size(); j++)
        {
            cout << vec[i][j] << "\t";
        }
        cout << endl;
    }
}

double vec_dot(vector<double> a, vector<double> b)
{
    double product = 0.0;
    for (int i = 0; i < a.size(); i++)
    {
        product = product + a[i] * b[i];
    }
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

vector<double> vec_cross(vector<double> a, vector<double> b)
{
    vector<double> temp = {0,0,0};
    temp[0] = a[1] * b[2] - a[2] * b[1];
    temp[1] = a[0] * b[2] - a[2] * b[0];
    temp[2] = a[0] * b[1] - a[1] * b[0];
    return temp;
}

double  vec_mag(vector<double> a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

vector<double> unit_vec(vector<double> a)
{
    vector<double> temp = {0,0,0};
    temp[0] = a[0]/vec_mag(a);
    temp[1] = a[1]/vec_mag(a);
    temp[2] = a[2]/vec_mag(a);
    return temp;
}

vector<double> make_vec(double x2, double y2, double z2, double x1, double y1, double z1)
{
    vector<double> temp = {x2-x1, y2-y1, z2-z1};
    return temp;
}

vector<double> vec_scalar_mult(double k, vector<double> a)
{
    vector<double> temp = {a[0]*k, a[1]*k, a[2]*k};
    return temp;
}

vector<double> vec_sum(vector<double> a, vector<double> b)
{
    vector<double> temp = {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
    return temp;
}

vector<double> vec_difference(vector<double> a, vector<double> b)
{
    vector<double> temp = {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
    return temp;
}

double vec_angle(vector<double> a, vector<double> b)
{
    double angle;
    angle = acos(vec_dot(a,b)/(vec_mag(a)*vec_mag(b)));
    return angle;
}
#endif 
