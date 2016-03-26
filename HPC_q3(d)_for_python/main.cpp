//
//  main.cpp
//  HPC_q3(d)_for_python
//
//  Created by Ernest on 26/3/2016.
//  Copyright © 2016 Kan Tsz Hang. All rights reserved.
//

#include <iostream>
#include <vector>
#include "TriMatrix.h"
#include <fstream>

using namespace std;

int main(int argc, const char * argv[]){
    
    double L=atof(argv[1]);
    int Nx=atof(argv[2]);
    double T=atof(argv[3]);
    double Nt=atof(argv[4]);
    double alpha=atof(argv[5]);
    double theta=atof(argv[6]);
    double dt=T/Nt;
    double dx = L/Nx;
    double nu = alpha*dt/dx/dx;
    double gamma0 = 0;
    double gamma1 = 0;
    
    
    //Defining all BASIC vectors used in the programme
    
    vector<double> x((Nx+1));
    vector<double> u0((Nx+1));
    vector <double> u1(Nx+1);
    vector<double>  u2(Nx+1);
    
    //Defining x-vector
    for (int i=0;i<Nx+1;i++) {
        x[i]=dx*i; //equallly spaced elements in x-vector
    }
    
    //Defining a vector, which uses the initial condition u0, to hold the solution
    for (int i=0; i<Nx+1; i++) {
        u0[i]=x[i]*(1-x[i]);
    }
    
    //Transfer entries from u0 vector to the u1 vector created earlier
    
    u1[0]=gamma0;
    u1[Nx]=gamma1;
    
    for (int i=1; i<Nx; i++) {
        
        u1[i] = u0[i];
    }
    
    
    
    //create tridiagonal matrix using TriMatrix class constructor
    
    TriMatrix Left(Nx, theta*nu*-1);
    TriMatrix Right(Nx, (1-theta)*nu);
    
    
    
    //forward Euler time integration via matrix multiplication
    
    for (double i=0; i<(T-dt); i+=dt){
        
        u2=Left/Right.multi(u1);
        u1=u2;}
    

    }
    
ofstream fout;
fout.open("results.txt");
fout << (*u1)[(Nx+1)/2] << endl;
fout.close();
cout<<"done"<<endl;

    return 0;

    
}
