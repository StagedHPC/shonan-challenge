#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

using namespace std;

const int NX = 300;
const int NY = 300;
const int NZ = 30;

const double Dt = 0.02;
const double Dx = 0.1;
const double Dy = 0.1;
const double Dz = 0.1;

double Sxx[NZ][NY][NX];
double Syy[NZ][NY][NX];
double Szz[NZ][NY][NX];
double Sxy[NZ][NY][NX];
double Sxz[NZ][NY][NX];
double Syz[NZ][NY][NX];

double Vx[NZ][NY][NX];
double Vy[NZ][NY][NX];
double Vz[NZ][NY][NX];

double dSxx_dx[NZ][NY][NX];
double dSxz_dx[NZ][NY][NX];
double dSxy_dx[NZ][NY][NX];
double dSyy_dy[NZ][NY][NX];
double dSxy_dy[NZ][NY][NX];
double dSyz_dy[NZ][NY][NX];
double dSzz_dz[NZ][NY][NX];
double dSxz_dz[NZ][NY][NX];
double dSyz_dz[NZ][NY][NX];


double dVx_dx[NZ][NY][NX];
double dVy_dx[NZ][NY][NX];
double dVz_dx[NZ][NY][NX];
double dVx_dy[NZ][NY][NX];
double dVy_dy[NZ][NY][NX];
double dVz_dy[NZ][NY][NX];
double dVx_dz[NZ][NY][NX];
double dVy_dz[NZ][NY][NX];
double dVz_dz[NZ][NY][NX];

const double r40 = 9.0/8.0;
const double r41 = 1.0/24.0;



template <class T> T sq(const T &x) { return x*x; }

int init() {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	Sxx[k][j][i] = 0;
	Syy[k][j][i] = 0;
	Szz[k][j][i] = 0;
	Sxy[k][j][i] = 0;
	Sxz[k][j][i] = 0;
	Syz[k][j][i] = 0;

	Vx[k][j][i] = 0;
	Vy[k][j][i] = 0;
	Vz[k][j][i] = 0.1 * exp(- (sq(i-0.5*NX) + sq(j-0.5*NY) + sq(k-0.5*NZ)) 
				 / 100);
      }
    }
  }

}

void diffx3_m4(double f[NZ][NY][NX], double (&df_dx)[NZ][NY][NX] ) {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 2; i < NX-1; ++i) {
	df_dx[k][j][i]
	  = ( (f[k][j][i  ] - f[k][j][i-1] ) * r40 
	    - (f[k][j][i+1] - f[k][j][i-2] ) * r41) / Dx;
      }
    }
  }  
}

void diffx3_p4(double f[NZ][NY][NX], double (&df_dx)[NZ][NY][NX] ) {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 1; i < NX-2; ++i) {
	df_dx[k][j][i]
	  = ( (f[k][j][i+1] - f[k][j][i  ] ) * r40 
	    - (f[k][j][i+2] - f[k][j][i-1] ) * r41) / Dx;
      }
    }
  }  
}

void diffy3_m4(double f[NZ][NY][NX], double (&df_dy)[NZ][NY][NX] ) {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 2; j < NY-1; ++j) {
      for (int i = 0; i < NX; ++i) {
	df_dy[k][j][i]
	  = ( (f[k][j  ][i] - f[k][j-1][i] ) * r40 
	    - (f[k][j+1][i] - f[k][j-2][i] ) * r41) / Dy;              
      }
    }
  }  
}

void diffy3_p4(double f[NZ][NY][NX], double (&df_dy)[NZ][NY][NX] ) {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 1; j < NY-2; ++j) {
      for (int i = 0; i < NX; ++i) {
	df_dy[k][j][i]
	  = ( (f[k][j+1][i] - f[k][j  ][i] ) * r40 
	    - (f[k][j+2][i] - f[k][j-1][i] ) * r41) / Dy;          
      }
    }
  }  
}

void diffz3_m4(double f[NZ][NY][NX], double (&df_dz)[NZ][NY][NX] ) {
  for (int k = 2; k < NZ-1; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	df_dz[k][j][i]
	  = ( (f[k  ][j][i] - f[k-1][j][i] ) * r40 
	    - (f[k+1][j][i] - f[k-2][j][i] ) * r41) / Dz;
      }
    }
  }  
}

void diffz3_p4(double f[NZ][NY][NX], double (&df_dz)[NZ][NY][NX] ) {
  for (int k = 1; k < NZ-2; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	df_dz[k][j][i]
	  = ( (f[k+1][j][i] - f[k  ][j][i] ) * r40 
	    - (f[k+2][j][i] - f[k-1][j][i] ) * r41) / Dz;
      }
    }
  }  
}


void periodic (double (&f)[NZ][NY][NX]) {
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	int i2 = i, j2=j, k2=k;
	if (i<2) i2 = NX-4+i;
	if (i>=NX-2) i2 = i-NX+2;
	if (j<2) j2 = NY-4+j;
	if (j>=NY-2) j2 = j-NY+2;
	if (k<2) k2 = NZ-4+k;
	if (k>=NZ-2) k2 = k-NZ+2;
	f[k][j][i] = f[k2][j2][i2];
      }
    }
  }
}

void diff_V() {
  diffx3_m4( Vx, dVx_dx );
  diffy3_p4( Vx, dVx_dy );
  diffz3_p4( Vx, dVx_dz );
  diffx3_p4( Vy, dVy_dx );
  diffy3_m4( Vy, dVy_dy );
  diffz3_p4( Vy, dVy_dz );
  diffx3_p4( Vz, dVz_dx );
  diffy3_p4( Vz, dVz_dy );
  diffz3_m4( Vz, dVz_dz );
  periodic(Vx);
  periodic(Vy);
  periodic(Vz);
  
}

void diff_S() {
  diffx3_p4( Sxx, dSxx_dx );
  diffy3_p4( Syy, dSyy_dy );
  diffx3_m4( Sxy, dSxy_dx );
  diffx3_m4( Sxz, dSxz_dx );
  diffy3_m4( Sxy, dSxy_dy );
  diffy3_m4( Syz, dSyz_dy );
  diffz3_p4( Szz, dSzz_dz );
  diffz3_m4( Sxz, dSxz_dz );
  diffz3_m4( Syz, dSyz_dz );
  periodic(Sxx);
  periodic(Syy);
  periodic(Szz);
  periodic(Sxy);
  periodic(Sxz);
  periodic(Syz);
}

void update_V() {
  double den = 1.0;
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	Vx[k][j][i] += Dt / den * 
	  (dSxx_dx[k][j][i] + dSxy_dy[k][j][i] + dSxz_dz[k][j][i]);
	Vy[k][j][i] += Dt / den * 
	  (dSxy_dx[k][j][i] + dSyy_dy[k][j][i] + dSyz_dz[k][j][i]);
	Vz[k][j][i] += Dt / den * 
	  (dSxz_dx[k][j][i] + dSyz_dy[k][j][i] + dSzz_dz[k][j][i]);
      }
    }
  }
}

void update_S() {
  double lam = 1.0;
  double mu = 0.35;
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {
	double divV = dVx_dx[k][j][i]
	  + dVy_dy[k][j][i] + dVz_dz[k][j][i];
	Sxx[k][j][i] +=  Dt * 
	  (2*mu * dVx_dx[k][j][i] + lam * divV);
	Syy[k][j][i] +=  Dt * 
	  (2*mu * dVy_dy[k][j][i] + lam * divV);
	Szz[k][j][i] +=  Dt * 
	  (2*mu * dVz_dz[k][j][i] + lam * divV);
	Sxy[k][j][i] +=  Dt * mu * 
	  (dVx_dy[k][j][i] + dVy_dx[k][j][i]);
	Sxz[k][j][i] +=  Dt * mu * 
	  (dVx_dz[k][j][i] + dVz_dx[k][j][i]);
	Syz[k][j][i] +=  Dt * mu * 
	  (dVy_dz[k][j][i] + dVz_dy[k][j][i]);
      }
    }
  }
}


void dump () {
  ofstream ofs("debug.txt");
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {
      ofs << i << " " << j << " " << Vz[NZ/2][j][i] << endl;
    }
    ofs << endl;
  }
}

int main () {
  cout << "hello" << endl;

  init();
  system("mkdir -p frame");


  for (int n = 0; n < 100; ++n) {
    cout << n << " " << Vz[NZ/2][NY/2][NX/2] << endl;
    for (int m = 0; m < 10; ++m) {
      diff_V();
      update_S();
      diff_S();
      update_V();
    }
    
    dump();

    system("gnuplot plot.gnu");
    char cmd[256];
    sprintf(cmd, "mv debug.png frame/%04d.png",n);
    system(cmd);
  }

  return 0;
}
