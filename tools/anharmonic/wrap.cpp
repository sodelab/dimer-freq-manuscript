#include <cmath>
#include <cassert>
#include <cstdlib> 
#include <vector>
#include <string>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "constants.h"
#include "kit.h"

void print_molpro_data(const size_t natoms, const double *ddE, 
                       const size_t ncoord, const double *mass,
                       const double *z,     const double *crd,  
                       const double *freq,  double **Evecs) 
{

  std::vector<std::string> elements(natoms);
  std::vector<std::string> hessians(natoms*3);

  const char ff[] = "molpro.out";
  FILE *fm = fopen(ff, "w");

  fprintf(fm,"                                         ***  PROGRAM SYSTEM MOLPRO  ***\n");
  fprintf(fm,"                                       Copyright, TTI GmbH Stuttgart, 2015\n");
  fprintf(fm,"                                    Version 2015.1 linked Feb  4 2016 13:46:49\n");
  fprintf(fm,"\n");
  fprintf(fm," Point group  C1\n");
  fprintf(fm,"\n");
  fprintf(fm,"\n");
  fprintf(fm," ATOMIC COORDINATES\n");
  fprintf(fm,"\n");
  fprintf(fm," NR  ATOM    CHARGE       X              Y              Z\n\n");
 
  std::string element;
  size_t ccount=0; 
  size_t ocount=0;
  char tmp[3];
  for (size_t i=0; i<natoms; ++i) {
    if (i%3==0) { ccount++; fprintf(fm,"   %zu  C%zu",i+1,ccount); sprintf(tmp,"C%zu",ccount); }
    else { ocount++; fprintf(fm,"   %zu  O%zu",i+1,ocount); sprintf(tmp,"O%zu",ocount); }
    std::string g(tmp);
    elements[i]=g;    
    fprintf(fm,"    %6.2f %14.9f %14.9f %14.9f\n",z[i],crd[i*3+0]/0.529177,crd[i*3+1]/0.529177,crd[i*3+2]/0.529177);
  }
  
  size_t count=0;
  char var[6];

  for (size_t i=0; i<natoms; ++i) {
    count++;
    for (size_t j=0; j<3; ++j) {
      if (i%3==0) {
        if (j==0) sprintf(var,"%s%zu%s%zu","C",count,"X",count);
        else if (j==1) sprintf(var,"%s%zu%s%zu","C",count,"Y",count);
        else if (j==2) sprintf(var,"%s%zu%s%zu","C",count,"Z",count);
      } else {
        if (j==0) sprintf(var,"%s%zu%s%zu","O",count,"X",count);
        else if (j==1) sprintf(var,"%s%zu%s%zu","O",count,"Y",count);
        else if (j==2) sprintf(var,"%s%zu%s%zu","O",count,"Z",count);
      }
      std::string s(var);
      hessians[i*3+j]=s;
//      std::cout << s << std::endl;
//      std::cout << hessians[i*3+j] << std::endl;
      
    }
  }
    
  fprintf(fm,"\n");
  fprintf(fm,"***\n");
  fprintf(fm,"\n");
  fprintf(fm," PROGRAM * HESSIAN\n");
  fprintf(fm,"\n");
  fprintf(fm," Force Constants (Second Derivatives of the Energy) in [a.u.]\n");

  int leni=ncoord; // length of matrix row size
  int lenj=ncoord; // length of matrix column size

  int l=0; int k=0; // row and column index
  int cstop=0; // column stop
  //int rstop; // row stop

  while (l<leni) {
    cstop+=5; // stop length

    int indexi=l;

    fprintf(fm,"    ");

    for (int tmpi=0; tmpi<5; tmpi++) { // maximum five columns
      indexi+=1;
      if (indexi<=cstop && indexi<=leni) {
//        fprintf(fm,"           %3d",indexi);
        fprintf(fm,"        %6s",hessians[indexi-1].c_str());
      }
    }

    int chk=l; // save row counter

    for (k=0; k<lenj; k++) { // all rows
      if (k>=cstop-5) { // only include necessary columns
//        fprintf(fm,"\n %3d",k+1);
        fprintf(fm,"\n %6s",hessians[k].c_str());
 
        l=chk; // restart row counter

        for (int tmpi=0; tmpi<5; tmpi++) { // maximum five columns
          if (l<=k) {
            fprintf(fm,"  %12.4E",ddE[l*ncoord+k]);
          }
          l+=1; //next row
        }
      }
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"\n\n");

  fprintf(fm," Atomic Masses\n");
  fprintf(fm,"   1- %2zu  ",natoms);
  for(size_t i=0; i<natoms; ++i) fprintf(fm,"%10.7f  ",mass[i]);
  fprintf(fm,"\n");
  
  fprintf(fm," Mass weighted Second Derivative Matrix\n");

  leni=ncoord; // length of matrix row size 
  lenj=ncoord; // length of matrix column size   
           
  l=0; k=0; // row and column index
  cstop=0; // column stop
  //rstop=0; // row stop

  while (l<leni) {
    cstop+=5; // stop length
   
    int indexi=l;
 
    fprintf(fm,"    "); 
         
    for (int tmpi=0; tmpi<5; tmpi++) { // maximum five columns
      indexi+=1;
      if (indexi<=cstop && indexi<=leni) {
        fprintf(fm,"        %6s",hessians[indexi-1].c_str());
                   
      }  
    }    
         
    int chk=l; // save row counter
                   
    for (k=0; k<lenj; k++) { // all rows
      if (k>=cstop-5) { // only include necessary columns
//        fprintf(fm,"\n %3d",k+1);
        fprintf(fm,"\n %6s",hessians[k].c_str());

        l=chk; // restart row counter
      
        for (int tmpi=0; tmpi<5; tmpi++) { // maximum five columns
          if (l<=k) {
            fprintf(fm,"  %12.4F",ddE[l*ncoord+k])/sqrt(mass[l/3])/sqrt(mass[k/3]);
          }
          l+=1; //next row
        } 
      }     
    }     
    fprintf(fm,"\n");
  }     
  fprintf(fm,"\n\n");

  fprintf(fm," Mass Weighted 2nd Derivative Matrix Eigenvalues\n");
  fprintf(fm,"   1- %2zu  ",ncoord);
  for(size_t i=0; i<ncoord; ++i) {
    if (freq[i]<0) fprintf(fm,"%10.7f  ",-sqrt(std::abs(freq[i])));
    else fprintf(fm,"%10.7f  ",sqrt(std::abs(freq[i])));
  }
  fprintf(fm,"\n\n");

  cstop=0; // column stop

  fprintf(fm," Mass Weighted 2nd Derivative Matrix Eigenvectors\n");

  l=0; k=0;
  while (l<leni) {
    cstop+=8; // stop length

    int indexi=l;

    fprintf(fm,"    ");
  
    for (int tmpi=0; tmpi<8; tmpi++) { // maximum eight columns
      indexi+=1;
      if (indexi<=cstop && indexi<=leni) {
        fprintf(fm,"           %3d",indexi);
//        fprintf(fm,"        %6s",hessians[indexi-1].c_str());
      }
    }
    int chk=l; // save row counter
                 
    for (k=0; k<lenj; k++) { // all rows
//      if (k>=cstop-8) { // only include necessary columns
        fprintf(fm,"\n %3d",k+1);
//        fprintf(fm,"\n %6s",hessians[k].c_str());
      
        l=chk; // restart row counter
        
        for (int tmpi=0; tmpi<8; tmpi++) { // maximum five columns
          if (l<leni) {
            fprintf(fm,"  %12.4E",Evecs[k][l]);//ddE[l*ncoord+k])/sqrt(mass[l/3])/sqrt(mass[k/3]);
          }
          l+=1; //next row
        }
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"\n\n");

  fprintf(fm,"   Low Vibration      Wavenumber\n");
  fprintf(fm,"        Nr             [1/cm]\n");
  size_t q=0;
  for(size_t i=0; i<ncoord; ++i) {
    if (freq[i]<0.0005) fprintf(fm,"    %5zu             %7.3f\n",++q,sqrt(std::abs(freq[i]))*5140.4862);
  }

  fprintf(fm,"\n\n");
  fprintf(fm,"     Vibration        Wavenumber\n");
  fprintf(fm,"        Nr             [1/cm]\n");
  q=0;
  for(size_t i=0; i<ncoord; ++i) {
    if (freq[i]>=0.0005) fprintf(fm,"    %5zu             %7.3f\n",++q,sqrt(std::abs(freq[i]))*5140.4862);
  }

  fprintf(fm," FREQUENCIES * CALCULATION OF NORMAL MODES FOR HF-SCF\n");

  fprintf(fm,"\n\n");
  fprintf(fm," Atomic Coordinates\n");
  fprintf(fm,"\n");
  fprintf(fm,"  Nr  Atom  Charge       X              Y              Z\n");
  fprintf(fm,"\n");

  for (size_t i=0; i<natoms; ++i){
    fprintf(fm," %3zu %4s  %4.2f  %13.10f  %13.10f  %13.10f\n",i+1,elements[i].c_str(),z[i],
							      crd[i*3+0]/0.529177,
							      crd[i*3+1]/0.529177,
							      crd[i*3+2]/0.529177);
  }
  fprintf(fm,"\n");
  
  fprintf(fm," Frequencies dumped to record   5400.2\n");
  fprintf(fm,"\n");
  fprintf(fm," Gradient norm at reference geometry: 0.16400D-04\n");
  fprintf(fm,"\n");
  fprintf(fm," Normal Modes\n");

//  fprintf(fm,"\n");

  cstop=0; // column stop
  l=0; k=0;
  while (l<leni) {
    cstop+=5; // stop length

    int indexi=l;

    fprintf(fm,"\n                       ");

    for (int tmpi=0; tmpi<5; tmpi++) { // maximum eight columns
      indexi+=1;
      if (indexi<=cstop && indexi<=leni) {
        fprintf(fm,"         %3d",indexi);
//        fprintf(fm,"        %6s",element[indexi-1].c_str());
      }
    }
    int chk=l; // save row counter
    int indexii=l;

    fprintf(fm,"\n");
    fprintf(fm," Wavenumbers [cm-1]    ");
    for(int tmpi=0; tmpi<5; tmpi++){
      indexii+=1;
      if(indexii<=cstop && indexii<=leni){
        fprintf(fm,"     %7.2f",sqrt(std::abs(freq[indexii-1]))*5140.4862);    
      }
    }
    fprintf(fm,"\n");
    fprintf(fm," Intensities [km/mol]  ");
    indexii=l;
    for(int tmpi=0; tmpi<5; tmpi++){
      indexii+=1;
      if(indexii<=cstop && indexii<=leni){
        fprintf(fm,"         0.0");    
      }
    }
    fprintf(fm,"\n");
    fprintf(fm," Intensities [relative]");
    indexii=l;
    for(int tmpi=0; tmpi<5; tmpi++){
      indexii+=1;
      if(indexii<=cstop && indexii<=leni){
        fprintf(fm,"         0.0");                                             
      }
    }

    for (k=0; k<lenj; k++) { // all rows
//        fprintf(fm,"\n %3d",k+1);
        fprintf(fm,"\n         %6s        ",hessians[k].c_str());

        l=chk; // restart row counter

        for (int tmpi=0; tmpi<5; tmpi++) { // maximum five columns
          if (l<leni) {
            fprintf(fm,"  %10.5f",Evecs[k][l]);//ddE[l*ncoord+k])/sqrt(mass[l/3])/sqrt(mass[k/3]);
          }
          l+=1; //next row
        }
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"\n\n");

  fprintf(fm," **********************************************************************************************************************************\n");
  fprintf(fm," Variable memory released\n");

  return;
}

void print_nwchem_data(const size_t natoms, const double *ddE,
                       const size_t ncoord, const double *mass,
                       const double *z,     const double *crd,
                       const double *freq,  double **Evecs)
{

  std::vector<std::string> elements(natoms);
  std::vector<std::string> hessians(natoms*3);

  const char ff[] = "nwchem.out";
  FILE *fm = fopen(ff, "w");

  fprintf(fm,"\n");
  fprintf(fm," ---------------------------- Atom information ----------------------------\n");
  fprintf(fm,"     atom    #        X              Y              Z            mass\n");
  fprintf(fm," --------------------------------------------------------------------------\n");
  
  for (size_t i=0; i<natoms; ++i) {
    fprintf(fm,"  %3s      %3d %14.7E %14.7E %14.7E  %13.7E\n",i%3 == 0 ? "C" : "O",i,
							       crd[i*3+0]/0.529177,crd[i*3+1]/0.529177,crd[i*3+2]/0.529177,
							       mass[i]);
  }
  fprintf(fm," --------------------------------------------------------------------------\n");
  fclose(fm);

  return;

}

void print_nmode(const size_t natoms, double **Evecs,
                 const size_t ncoord, const double *mass,  
                 double *freq,         const size_t nm, 
                 const size_t na,      const int *ac)
{

  const char fm[] = "nmode";
  FILE *fs = fopen(fm, "w");
 
  for (size_t f=0; f<ncoord-nm; ++f) {
    if (freq[f]<0) fprintf(fs,"%15.10f\n",-sqrt(std::abs(freq[f]))*5140.4862);
    else fprintf(fs,"%15.10f\n",sqrt(std::abs(freq[f]))*5140.4862);
  } 
  for (size_t f=0; f<na; ++f) {
    size_t n=ncoord-nm+ac[f];
    fprintf(fs,"%15.10f\n",sqrt(std::abs(freq[n]))*5140.4862);
  }
  for (size_t f=0; f<nm-na; ++f) {
    fprintf(fs,"%15.10f\n",0.00);
  }
 

  for (size_t n=0; n<ncoord-nm; ++n) {
    for (size_t j=0; j<ncoord; ++j) {
      fprintf(fs,"%15.10f\n",Evecs[j][n]);//sqrt(mass[j/3]));
    }
  }
  for (size_t f=0; f<na; ++f) {
    size_t n=ncoord-nm+ac[f];
    for (size_t i=0; i<ncoord; ++i) {
      fprintf(fs,"%15.10f\n",Evecs[i][n]);
    }
  }

  return;
}

void print_sindo(const size_t natoms, const double *ddE,
                 const size_t ncoord, const double *mass,
                 const double *z,     const double *crd,
                 const double *freq,  double **Evecs,
                 const size_t nm,     const size_t na, 
                 const int *ac)
{

  const char fm[] = "sindo.inp";
  FILE *fs = fopen(fm, "w");

//  if (natoms==3) fprintf(fs,"&mol Nat= %3d Nfree= %3d\n",natoms,3*natoms-5);
//  else fprintf(fs,"&mol Nat= %3d Nfree= %3d\n",natoms,3*natoms-6);
  fprintf(fs,"&mol Nat= %3d Nfree= %3d\n",natoms,na);

  fprintf(fs,"mass=\n");  
  for (size_t n=0; n<natoms; ++n) 
    fprintf(fs,"%9.4f\n",mass[n]);

  fprintf(fs,"omega=\n");
//  for (size_t f=0; f<ncoord-nm; ++f) {
//    fprintf(fs,"%9.4f\n",sqrt(std::abs(freq[f]))*5140.4862);
//  } 
  for (size_t f=0; f<na; ++f) {
    size_t n=ncoord-nm+ac[f];
    fprintf(fs,"%9.4f\n",sqrt(std::abs(freq[n]))*5140.4862);   
  }
//  if (natoms==3) {
//    for (size_t f=5; f<natoms*3; ++f) {
//      fprintf(fs,"%9.4f\n",sqrt(std::abs(freq[f]))*5140.4862);
//    }
//  } else {
//    for (size_t f=6; f<natoms*3; ++f) {
//      fprintf(fs,"%9.4f\n",sqrt(std::abs(freq[f]))*5140.4862);
//    }
//  }

  fprintf(fs,"x=\n");
  for (size_t x=0; x<natoms; ++x)
    fprintf(fs,"%9.4f  %9.4f  %9.4f\n",crd[x*3+0],crd[x*3+1],crd[x*3+2]);

  fprintf(fs,"L=\n");
  for (size_t f=0; f<ncoord-nm; ++f) {
    for (size_t i=0; i<ncoord; ++i) {
      fprintf(fs,"%9.4f",Evecs[i][f]);
      if (i%3==2) fprintf(fs,"\n");
    }
  }
  for (size_t f=0; f<na; ++f) {
    size_t n=ncoord-nm+ac[f];
    for (size_t i=0; i<ncoord; ++i) {
      fprintf(fs,"%9.4f",Evecs[i][n]);
      if (i%3==2) fprintf(fs,"\n");
    }
  }
  
//  if (natoms==3) {
//    for (size_t l=5; l<natoms*3; ++l) {
//      for (size_t i=0; i<ncoord; ++i) {
//        fprintf(fs,"%9.4f",Evecs[i][l]);
//        if (i%3==2) fprintf(fs,"\n");
//      }
//    }
//  } else {
//    for (size_t l=6; l<natoms*3; ++l) {
//      for (size_t i=0; i<ncoord; ++i) {
//        fprintf(fs,"%9.4f",Evecs[i][l]);
//        if (i%3==2) fprintf(fs,"\n");
//      }
//    }
//  }
  return;
}

void distance(const double *a, const double *b, double &r) 
{
  
  double x,y,z;
  
  x=b[0]-a[0];
  y=b[1]-a[1];
  z=b[2]-a[2];

  r=sqrt(x*x+y*y+z*z);
  return;  

}

void angle(const double *a, const double *b, const double *c, double &s)
{

  double x1,y1,z1;
  double x2,y2,z2;

  x1=a[0]-b[0];y1=a[1]-b[1];z1=a[2]-b[2];
  x2=c[0]-b[0];y2=c[1]-b[1];z2=c[2]-b[2];

  double r1=sqrt(x1*x1+y1*y1+z1*z1);
  double r2=sqrt(x2*x2+y2*y2+z2*z2);  

  double v1norm[3];
  double v2norm[3];

  v1norm[0]=x1/r1;v1norm[1]=y1/r1;v1norm[2]=z1/r1;
  v2norm[0]=x2/r2;v2norm[1]=y2/r2;v2norm[2]=z2/r2;

  double res=v1norm[0]*v2norm[0] + v1norm[1]*v2norm[1] + v1norm[2]*v2norm[2];
  
  if (res<-1.0) res=-1.0;
  else if (res>1.0) res=1.0;

  s=acos(res)*180.0/M_PI;
  return;

}

void dihedral(const double *a, const double *b, const double *c, const double *d, double &t)
{

  double x1,y1,z1;
  double x2,y2,z2;
  double x3,y3,z3;

  x1=a[0]-b[0];y1=a[1]-b[1];z1=a[2]-b[2];
  x2=c[0]-b[0];y2=c[1]-b[1];z2=c[2]-b[2];
  x3=d[0]-c[0];y3=d[1]-c[1];z3=d[2]-c[2];

  double r1=sqrt(x1*x1+y1*y1+z1*z1);
  double r2=sqrt(x2*x2+y2*y2+z2*z2);  
  double r3=sqrt(x3*x3+y3*y3+z3*z3);  
  
  double v1norm[3],v1[3];
  double v2norm[3],v2[3];
  double v3norm[3],v3[3];

  double w[3],v[3];
  
  v1[0]=x1;v1[1]=y1;v1[2]=z1;
  v2[0]=x2;v2[1]=y2;v2[2]=z2;
  v3[0]=x3;v3[1]=y3;v3[2]=z3;

  v1norm[0]=x1/r1;v1norm[1]=y1/r1;v1norm[2]=z1/r1;
  v2norm[0]=x2/r2;v2norm[1]=y2/r2;v2norm[2]=z2/r2;
  v3norm[0]=x3/r3;v3norm[1]=y3/r3;v3norm[2]=z3/r3;

  double vdot=v1[0]*v2norm[0]+v1[1]*v2norm[1]+v1[2]*v2norm[2];
  double wdot=v3[0]*v2norm[0]+v3[1]*v2norm[1]+v3[2]*v2norm[2];

  for (size_t i=0; i<3; i++) {
    v[i]=vdot*v2norm[i];
    w[i]=wdot*v2norm[i];
   
    v[i]=v1[i]-v[i];
    w[i]=v3[i]-w[i];
  }

  double yyy[3];

  double xdot=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
  yyy[0]=v2norm[1]*v[2]-v2norm[2]*v[1];
  yyy[1]=v2norm[2]*v[0]-v2norm[0]*v[2];
  yyy[2]=v2norm[0]*v[1]-v2norm[1]*v[0];

  double ydot=yyy[0]*w[0]+yyy[1]*w[1]+yyy[2]*w[2];

  t=atan2(ydot,xdot)*180.0/M_PI;

  return; 

}

void print_nitrogen(const size_t natoms, const double *coord)
{

  const char fm[] = "nitrogen.dat";
  FILE *fs = fopen(fm, "w");

  //fprintf(fs,"&mol Nat= %3d Nfree= %3d\n",natoms,na);

  /*
  COORD = ZMAT
  %ZMAT
  C1  
  O2 1 R12
  O3 1 R13 2 A213
  C4 1 R14 2 A214 3 D3214
  O5 4 R45 1 A145 2 D2145 // 
  O6 4 R46 5 A546 3 D3546 //dihedreal needs debugginb
  % 
  */

  double r12,r13,r14,r45,r46; // distances
  double a213,a214,a145,a546; // angles
  double d3214,d2145,d3546; // dihedrals

  double c1[3],o2[3],o3[3];
  double c4[3],o5[3],o6[3];

  std::copy(coord + 0*3, coord + 0*3 + 3, c1);
  std::copy(coord + 1*3, coord + 1*3 + 3, o2);
  std::copy(coord + 2*3, coord + 2*3 + 3, o3);
  std::copy(coord + 3*3, coord + 3*3 + 3, c4);
  std::copy(coord + 4*3, coord + 4*3 + 3, o5);
  std::copy(coord + 5*3, coord + 5*3 + 3, o6);

  distance(o2,c1,r12);
  distance(o3,c1,r13);
  distance(c4,c1,r14);
  distance(o5,c4,r45);
  distance(o6,c4,r46);

  angle(o2,c1,o3,a213);
  angle(o2,c1,c4,a214);
  angle(c1,c4,o5,a145);
  angle(o5,c4,o6,a546);

  dihedral(o3,o2,c1,c4,d3214);
  dihedral(o2,c1,c4,o5,d2145);
  dihedral(o3,o5,c4,o6,d3546);

  fprintf(fs,"%20.17f %20.17f D%-20.15f %20.17f D%-20.15f D%-20.15f"
             "%20.17f D%-20.15f D%-20.15f %20.17f D%-20.15f D%-20.15f\n",
             r12,r13,a213,r14,a214,d3214,r45,a145,d2145,r46,a546,d3546);

  fclose(fs);
  return;  
  
}
