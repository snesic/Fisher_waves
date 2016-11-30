/*
 */

#include "read_write_msg.h"


void initial_message(int lx, int ly, long int t, double dt, double dx, int br_int, double sigma)
{
    cout<< "FISHER WAVES, Splitting step method 2D! " << endl;
    cout<< "System parameters: " << endl;
    cout<< endl;
    cout<< "System size: x coordinate= " << lx << "\t \t" << "y coordinate: " << ly << "\t \t" << "Time=  " << t*dt << endl;
    cout<< "Discretization: " << " dx= " << dx << "\t" << "dt= "<< dt <<endl;
    cout<< "Br interacija: " << br_int<<endl;
    cout << "Sigma= " << sigma << endl;
    cout<<endl;
}

void write1d_data(double *x, double *y, int size, string filename)
{
    ofstream dat;
    dat.open(filename.c_str(), ios::out);
    for(int i=0;i<size;i++)
        dat<< x[i] << "   " << y[i] << endl;
    dat.close();
}


void write2d_data_sq(double **mat, int ly, int no_points, int br_int, string filename)
{   ofstream dat;

    dat.open(filename.c_str(), ios::out);

    for(int iy=0;iy<ly/2;iy++)
    {   dat << iy << "\t";
        for(int i=0;i<no_points;i++)
            dat << mat[iy][i]/br_int << "\t";
        dat << "\n";
    }
    
    dat.close();

}

void write2d_data(double **mat, int ly, int no_points, int br_int, string filename)
{   ofstream dat;
    
    dat.open(filename.c_str(), ios::out);
    
    for(int iy=0;iy<ly;iy++)
    {  dat << iy << "\t";
        for(int i=0;i<no_points;i++)
            dat << mat[i][iy]/br_int << "\t";
        dat << "\n";
    }
    dat.close();
    
}
