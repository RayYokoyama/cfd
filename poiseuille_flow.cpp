#include<iostream>
#include<cmath>
#include<fstream>
#include <sstream>
#include <string>
#include<vector>
#define xn 88
#define yn 22
#define Re 100
using namespace std;

inline void update_velocity(
  vector<vector<double>> &u,
  vector<vector<double>> &v,
  vector<vector<double>> &p,
  double dt, 
  double dx,
  double dy,
  double &diff
) {
  for(int i=2; i<xn; i++){
    for(int j=1; j<yn; j++){
      u[j][i]=u[j][i]-dt*((p[j][i]-p[j][i-1])/dx);
    }
  }
  for(int i=1; i<xn; i++){
    for(int j=2; j<yn; j++){
      v[j][i]=v[j][i]-dt*((p[j][i]-p[j-1][i])/dy);
    }
  }

  diff=0.0;
  for(int i=1; i<xn; i++){
    for(int j=1; j<yn; j++){
      diff+=fabs((u[j][i+1]-u[j][i])/dx+(v[j+1][i]-v[j][i])/dy);
    }
  }
}

inline void solve_poisson(
  vector<vector<double>> &u, 
  vector<vector<double>> &v, 
  vector<vector<double>> &p,
  double dx,
  double dy,
  double dt,
  int km,
  double &err
){

  double C1=dy*dy/(2.0*(dx*dx+dy*dy));
	double C2=dx*dx/(2.0*(dx*dx+dy*dy));
	double C3=dx*dx*dy*dy/(2.0*(dx*dx+dy*dy))/dt;

  for(int k=0; k<=km; k++){
    // boundary condition Neumman
    for(int i=0; i<xn+1; i++){
      p[0][i]=p[1][i];
      p[yn][i]=p[yn-1][i];
    }
    for(int i=0; i<yn+1; i++){
      p[i][0]=100.0;
      p[i][xn]=0.0;
    }
    err=0.0;
    for(int i=1; i<xn; i++){
      for(int j=1; j<yn; j++){
        double origin=p[j][i];
        p[j][i]=C1*(p[j][i+1]+p[j][i-1])+C2*(p[j+1][i]+p[j-1][i])-C3*((u[j][i+1]-u[j][i])/dx+(v[j+1][i]-v[j][i])/dy);
        err+=(p[j][i]-origin)*(p[j][i]-origin);
      }
    }
    if(err<=0.001){
      break;
    }
  }
}

inline void predict_velocity(
  vector<vector<double>> &u,
  vector<vector<double>> &v,
  double dx,
  double dy,
  double dt
){
  //Boundary condition for right wall
  for(int i=0; i<yn+1; i++){
    u[i][xn+1]=u[i][xn];
  }
  //Boundary condition for left wall
  for(int i=0; i<yn+1; i++){
    u[i][0]=u[i][1];
  }
  //Boundary condition for bottom wall
  for(int i=0; i<xn+1; i++){
    v[1][i]=0.0;
    v[0][i]=v[2][i];
  }
  for(int i=0; i<xn+2; i++){
    u[0][i]=-u[1][i];
  }
  //Boundary condition for upper wall
  for(int i=0; i<xn+1; i++){
    v[yn][i]=0.0;
    v[yn+1][i]=v[yn-1][i];
  }
  for(int i=0; i<xn+2; i++){
    u[yn][i]=-u[yn-1][i];
  }
  // // predict u
  // for(int i=2; i<xn; i++){
  //   for(int j=1; j<yn; j++){
  //     double udif2 = (u[j][i+1]-2.0*u[j][i]+u[j][i-1])/(dx*dx)+(u[j+1][i]-2.0*u[j][i]+u[j-1][i])/(dy*dy);
  //     double udif1 = u[j][i]*(u[j][i+1]-u[j][i])/dx + v[j][i]*(u[j+1][i]-u[j][i])/dy;
  //     u[j][i] = u[j][i] + dt*(-udif1 + (1.0/Re)*udif2);
  //   }
  // }
  // 
  // // predict v
  // for(int i=2; i<xn; i++){
  //   for(int j=1; j<yn; j++){
  //     double vdif2 = (v[j][i+1]-2.0*v[j][i]+v[j][i-1])/(dx*dx)+(v[j+1][i]-2.0*v[j][i]+v[j-1][i])/(dy*dy);
  //     double vdif1 = u[j][i]*(v[j][i+1]-v[j][i-1])/(2*dx) + v[j][i]*(v[j+1][i]-v[j-1][i])/(2*dy);
  //     v[j][i] = v[j][i] + dt*(-vdif1 + (1.0/Re)*vdif2);
  //   }
  // }
  // predict u
  for(int i=2; i<xn; i++){
    for(int j=1; j<yn; j++){
      double udif2 = (u[j][i+1]-2.0*u[j][i]+u[j][i-1])/(dx*dx)+(u[j+1][i]-2.0*u[j][i]+u[j-1][i])/(dy*dy);
      double udif1 = u[j][i]*(u[j][i+1]-u[j][i-1])/(2*dx) + v[j][i]*(u[j+1][i]-u[j-1][i])/(2*dy);
      u[j][i] = u[j][i] + dt*(-udif1 + (1.0/Re)*udif2);
    }
  }

  // predict v
  for(int i=2; i<xn; i++){
    for(int j=1; j<yn; j++){
      double vdif2 = (v[j][i+1]-2.0*v[j][i]+v[j][i-1])/(dx*dx)+(v[j+1][i]-2.0*v[j][i]+v[j-1][i])/(dy*dy);
      double vdif1 = u[j][i]*(v[j][i+1]-v[j][i-1])/(2*dx) + v[j][i]*(v[j+1][i]-v[j-1][i])/(2*dy);
      v[j][i] = v[j][i] + dt*(-vdif1 + (1.0/Re)*vdif2);
    }
  }
}

int main()
{
  vector<vector<double>> p(yn+1, vector<double>(xn+1));
  vector<vector<double>> u(yn+1, vector<double>(xn+2));
  vector<vector<double>> v(yn+2, vector<double>(xn+1));
  double lx=4.0;
  double ly=1.0;
  double dx=lx/(xn-1);
  double dy=ly/(yn-1);
  int lm=1000;
  double diff;
  double dt=0.001;
  double err;
  int km=600;
  ofstream ff,fc,fu;
	ff.open("pressure2.vtk");

  int point_count=0;
  for(int i=0; i<xn; i++){
    for(int j=0; j<yn; j++){
      if(!(u[j][i]==0 && v[j][i]==0)){
        point_count++;
      }
    }
  }
  for(int l=0; l<=lm; l++){
    predict_velocity(u, v, dx, dy, dt);
    solve_poisson(u, v, p, dx, dy, dt, km, err);
    update_velocity(u, v, p, dt, dx, dy, diff);
    if(l%100==0) cout<<l<<" "<<err<<" "<<diff<<endl;

    if(l%50==0) {
      stringstream filename;
      filename << "velocity" << l/50 << ".vtk";
      string result = filename.str();
      ofstream fk;
      fk.open(result);
      fk << "# vtk DataFile Version 3.0" << endl;
      fk << "vtk output" << endl;
      fk << "ASCII" << endl;
      fk << "DATASET UNSTRUCTURED_GRID" << endl;
      fk << "POINTS " << 1913 << " double" << endl;
      for(int i=0; i<xn; i++){
        for(int j=0; j<yn; j++){
          if(!(u[j][i]==0 && v[j][i]==0)){
            fk<< i*dx << " " << j*dy << " " << 0 << endl;
          }
        }
      }
      fk << "POINT_DATA " << 1913 << endl;
      fk << "VECTORS velocity[m/s] double" << endl; 
      for(int i=0; i<xn; i++){
        for(int j=0; j<yn; j++){
          if(!(u[j][i]==0 && v[j][i]==0)){
            fk << u[j][i] << " " << v[j][i] << " " << 0 << endl;
          }
        }
      }
      fk.close();            
    }
  }

  ff << "# vtk DataFile Version 3.0" << endl;
  ff << "vtk output" << endl;
  ff << "ASCII" << endl;
  ff << "DATASET STRUCTURED_POINTS" << endl;
  ff << "DIMENSIONS " << xn-1 << " " << yn-1 << " " << 1 << endl;
  ff << "ORIGIN " << dx/2 << " " << dy/2 << " " << 0.0 << endl;
  ff << "SPACING " << dx << " " << dy << " " << 1.0 << endl;
  ff << "POINT_DATA " << (xn-1)*(yn-1) << endl;
  ff <<"SCALARS " << "pressure " << "double" << endl;
  ff <<"LOOKUP_TABLE default" << endl;
  for(int i=1; i<yn; i++){
    for(int j=1; j<xn; j++){
      ff << p[i][j] << endl;
    }
  }
  fc.open("y_coordinate.txt");
  fu.open("u.txt");       
  for(int i=1; i<yn; i++){
    fc << i*dy << endl;
    fu << u[i][xn/2] << endl;   
  } 
}
