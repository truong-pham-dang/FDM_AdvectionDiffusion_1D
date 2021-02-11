#pragma once
#include <cmath>
#ifndef M_PI
namespace
{
    const double M_PI = acos(-1.0);
}
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

/* Reference: https://www.dam.brown.edu/people/ntrask/CFD/CFD2016.html */

//Data structures
vector<double> x;
vector<double> uold;
vector<double> unew;

//Variables that we'll need to use
int N;
double t, tfinal, dt, dx; 
double a,nu; //advection and diffusion parameters;

//For writing data to disk
void writeFields() {

  ofstream write;

  write.open("num.dat");

  write << "TITLE = \"Numerical solution\"" << endl;
  write << "VARIABLES = \"x\", \"numerical solution\"" << endl;
  write << "ZONE T=\"Only Zone\", I=" << N << ", F=POINT" << endl;

  for(int i = 0; i < N; i++){
    write << x[i] << " " << unew[i] << endl;
  }

  write.close();

  write.open("exact.dat");

  write << "TITLE = \"Exact solution\"" << endl;
  write << "VARIABLES = \"x\", \"exact solution\"" << endl;
  write << "ZONE T=\"Only Zone\", I=" << N << ", F=POINT" << endl;

  for(int j = 0; j < N; j++){
    write << x[j] << " " << exp(-t)*sin(x[j]-t) << endl;
  }

  write.close();

  write.open("numVsExact.dat");

  write << "TITLE = \"Exact solution vs numerical solution\"" << endl;
  write << "VARIABLES = \"x\", \"exact solution\", \"numerical solution\"" << endl;
  write << "ZONE T=\"Only Zone\", I=" << N << ", F=POINT" << endl;

  for(int k = 0; k < N; k++){
    write << x[k] << " " << exp(-t)*sin(x[k]-t) << " " << unew[k] << endl;
  }

  write.close();

}

//Write a function that will update unew and uold at each timestep
void update(){
  
  double dxi = 1.0/dx; //Inverse of dx

  //First update the stencil

  //Interior points
  for(int i = 1; i < N-1; i++){
    unew[i] = uold[i] + dt*(
			    -(uold[i]-uold[i-1])*dxi + (uold[i+1] - 2.0*uold[i] + uold[i-1])*dxi*dxi
			    );
  }

  //End points
    unew[0]   = uold[0]   + dt*(
			    -(uold[0]-uold[N-1])*dxi + (uold[1] - 2.0*uold[0] + uold[N-1])*dxi*dxi
			    );
    unew[N-1] = uold[N-1] + dt*(
			    -(uold[N-1]-uold[N-2])*dxi + (uold[0] - 2.0*uold[N-1] + uold[N-2])*dxi*dxi
			    );
  
  
    //Now kick everything up one timestep
    for(int j = 0; j < N; j++){
      uold[j] = unew[j];
    }

}

int main() {

  cout << "Beginning finite difference code for 1D advection-diffusion equation with periodic BCs." << endl;
  
  cout << "Pre-processing: Initialize parameters and data structures." << endl;
  cout << "Input N:";
  cin >> N;

  dx = 2.0 * M_PI / double(N);
  dt = __min(0.5*dx,0.25*dx*dx);
  tfinal = 0.5;

  //Initialize data structures
  x.resize(N);
  uold.resize(N);
  unew.resize(N);

  for(int i = 0; i < N; i++){
    x[i] = double(i)*dx;
    uold[i] = sin(x[i]);
    unew[i] = 0.0;
  }
    
  cout << "Beginning timestepping." << endl;

  t = 0;

  while(t < tfinal){
	cout << "t = " << t << endl;
    update();
    t += dt;
  }

  cout << "Timestepping complete - outputing final state for post-processing." << endl;
  writeFields();

  cout << "Compute error compared to exact solution." << endl;

  double error = 0.0;

  for(int j = 0; j < N; j++){
    error += abs(unew[j] - exp(-t)*sin(x[j]-t));
  }
  cout << "(dt,dx,error) = (" << dt << "," << dx << "," << error/double(N) << ")"<< endl; 

  cout << "Finished." << endl;

  return 0;
}