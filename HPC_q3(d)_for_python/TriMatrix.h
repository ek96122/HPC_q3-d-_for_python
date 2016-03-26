//
//  TriMatrix.h
//  HPC_q3(d)_for_python
//
//  Created by Ernest on 26/3/2016.
//  Copyright © 2016 Kan Tsz Hang. All rights reserved.
//

#ifndef TriMatrix_h
#define TriMatrix_h

#include <vector>
#include <iostream>
using namespace std;

class TriMatrix{
    
private:
    vector<double> *diag;
    vector<double> *u_diag;
    vector<double> *l_diag;
    
public:
    //CONSTRUCTOR
    
    TriMatrix(vector <double> *p,vector<double> *q, vector<double> *r):
    
    l_diag(p), diag(q), u_diag(r) {};
    
    TriMatrix(double Nx, double nu){
        l_diag = new vector<double> (Nx);
        diag = new vector<double>(Nx+1);
        u_diag = new vector<double>(Nx);
        
        for (int i = 0; i < Nx + 1; i++){
            (*diag)[i] = 1-2*nu;
        }
        
        for (int i = 0; i < Nx; i++){
            (*l_diag)[i] = nu;
            (*u_diag)[i] = nu;
        }
        (*diag)[0] = 1;
        (*diag)[Nx] = 1;
        (*l_diag)[Nx-1] = 0;
        (*u_diag)[0] = 0;
    }
    
    //Create a function which can perform matrix multiplication for the tridiagonal matrix with any Nx1 matrix, where N=number of rows of the tridiagonal matrix
    
    vector<double> multi(vector <double> U){
        int D =diag->size();
        vector <double> U2(D),a(D),b(D),c(D);
        
        
        for (int i=0; i<D; i++){
            a[i]=(*diag)[i]*U[i];
        }
        for (int i=1; i<D; i++){
            b[i]=(*l_diag)[i-1]*U[i-1];
            
        }
        for (int i=0; i<(D-1); i++) {
            c[i]=(*u_diag)[i]*U[i+1];
        }
        for (int i=0; i<(D-1); i++) {
            U2[i] = a[i] + b[i] + c[i];
        }
        U2[0]=a[0]+b[0]; //compatibility check on 1st row
        U2[D-1]=a[D-1]+c[D-1];//compatibility check on last row
        
        return U2;
    }
    
    
    
    vector<double> operator/(vector <double> U){
        int D = diag->size();
        vector <double> newdiag(D);
        vector <double> U2(D);
        vector <double> low_long(D);
        double m;
        
        for (int i = 0; i < D-1; i++){
            low_long[i+1] = (*l_diag)[i];
        }
        
        newdiag[0] = (*diag)[0];
        
        for (int i=1; i<D ; i++){
            m=low_long[i]/newdiag[i-1];
            
            newdiag[i]=(*diag)[i]-m*(*u_diag)[i-1];
            
            U[i]=U[i]-m*U[i-1];
            
        }
        
        U2[D-1]=U[D-1]/newdiag[D-1];
        
        for (int i=D-2; i>=0; i--){
            U2[i]=(U[i]-(*u_diag)[i]*U2[i+1])/newdiag[i];
            
            
        }
        
        
        //cout << endl << "U2: " << endl;
        //for(int i=0; i<D; i++) cout << (U2)[i] << ", ";
        //cout<<endl;
        
        return U2;
        
        
    }
    
    
    
    
    
    void displayM() {
        cout << endl << "upper diagonal: " << endl;
        for(int i=0; i<(*u_diag).size(); i++) cout << (*u_diag)[i] << ", ";
        cout<<endl;
        
        cout << "lower diagonal: " << endl;
        for(int i=0; i<(*l_diag).size(); i++) cout << (*l_diag)[i] << ", ";
        cout<<endl;
        
        cout << endl << "diagonal: " << endl;
        for(int i=0; i<(*diag).size(); i++) cout << (*diag)[i] << ", ";
        cout<<endl;
        
    }
    
    
    
};

#endif /* TriMatrix_h */
