#include <cmath>
#include <cassert>
#include <cstdlib> 

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "constants.h"
#include "mbCO2.h"
#include "io-xyz.h"
#include "wrap.h"

#include "dsyev.h"

//#define TRAINING 0

#ifdef _OPENMP
#include "omp.h"
#endif

#define STEPSIZE 0.005

double get_energy(const double *crd, const size_t ncoord) 
{
  size_t natoms = ncoord/3;
  size_t nfrags = natoms/3;

  std::cout << std::scientific << std::setprecision(9);

  double pos[ncoord];

#ifdef TRAINING
  const char ff[] = "ts-vibration.xyz";
  FILE *fo = fopen(ff, "a");
#endif

  for (size_t n = 0; n < ncoord; ++n)
    pos[n] = crd[n];

  double E = 0.0;

#pragma omp parallel for schedule(dynamic) reduction(+:E) 
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    E += p1b(pos + 9*ifrag);
  }

#pragma omp parallel for schedule(dynamic) reduction(+:E) 
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    for(size_t jfrag=ifrag+1; jfrag<nfrags; jfrag++){

      double dim[18];
      std::copy(pos + ifrag*9, pos + ifrag*9 + 9, dim + 0);
      std::copy(pos + jfrag*9, pos + jfrag*9 + 9, dim + 9);

      double nnn = p2b(dim);
      E += nnn;

#ifdef TRAINING
        double x = dim[9+0] - dim[0+0];
        double y = dim[9+1] - dim[0+1];
        double z = dim[9+2] - dim[0+2];
        double d = sqrt ( x*x + y*y + z*z ); // C-C distance

 	fprintf(fo,"6\n");
   	fprintf(fo,"%2d %2d %lf %lf\n",ifrag,jfrag,d,nnn); 
 	for (size_t n = 0; n < 6; ++n)      
	    fprintf(fo,"%-2s %20.10f %20.10f %20.10f\n",n%3 == 0 ? "C" : "O",
							dim[3*n + 0],
							dim[3*n + 1],
							dim[3*n + 2]);
#endif

    }
  }

#ifdef TRAINING
  fclose(fo);  
#endif

  return E;
}

double get_derivative(const double *crd, const size_t ncoord,
		      const size_t n)
{

  double dE;

    double e_num[2];
    for (size_t i=0; i<2; ++i) { // two energy evals 

      double tmp_coord[ncoord];

      for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
        if (ndiff == n) {
          tmp_coord[ndiff] = (i ? crd[ndiff] + STEPSIZE : crd[ndiff] - STEPSIZE);
        } else {
          tmp_coord[ndiff] = crd[ndiff];
        }
      }
      e_num[i]=get_energy(tmp_coord, ncoord);

    }
    // do numerical differentiaions with energies
    dE=(e_num[1]-e_num[0])/(2*STEPSIZE);
    return dE;
    
}

void second_derivative(const double *crd, const size_t ncoord,
                       double *ddE)
{

  // do numerical derivative for each coordinate
  for (size_t n0=0; n0<ncoord; ++n0) {
    for (size_t n1=0; n1<ncoord; ++n1) {

      double d[2];
      for (size_t i=0; i<2; ++i) {

        double tmp_coord[ncoord];

        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
          if (ndiff == n0) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + STEPSIZE : crd[ndiff] - STEPSIZE);
          } else {
            tmp_coord[ndiff] = crd[ndiff];
          }
        }

        d[i]=get_derivative(tmp_coord,ncoord,n1);

      }
        
      ddE[n0*ncoord+n1]=(d[1]-d[0])/(2*STEPSIZE);
      ddE[n0*ncoord+n1]=ddE[n0*ncoord+n1]*constants::Bohr_A*constants::Bohr_A/constants::Eh_kcalmol;
    }
  }
  return;

}

double get_grad(const double *crd, const size_t ncoord,
                double const* const* Evecs, const size_t n2)
{

    double e[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n2]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n2]*STEPSIZE);
        }
        
        e[i]=get_energy(tmp_coord, ncoord);
        
    }
    
    return (e[1]-e[0])/(2*STEPSIZE);
    
}

double get1MR(const double *crd, const size_t ncoord,
              double **Evecs, const size_t n2)
{
    
    double e[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n2]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n2]*STEPSIZE);
        }
        
        e[i]=get_energy(tmp_coord, ncoord);
        
    }

    return (e[1]-e[0])/(2*STEPSIZE);
    
}

double get_hess(const double *crd, const size_t ncoord,
                double const* const* Evecs,
                const size_t n1, const size_t n2)
{
    
    double d[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n1]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n1]*STEPSIZE);
        }
        
        d[i]=get_grad(tmp_coord,ncoord,Evecs,n2);
        
    }
    
    return (d[1]-d[0])/(2*STEPSIZE);
    
}

double get2MR(const double *crd, const size_t ncoord,
              double **Evecs,
              const size_t n1, const size_t n2)
{
    
    double d[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n1]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n1]*STEPSIZE);
        }
        
        d[i]=get1MR(tmp_coord,ncoord,Evecs,n2);
        
    }
    
    return (d[1]-d[0])/(2*STEPSIZE);
    
}

double get_cub(const double *crd, const size_t ncoord,
               double const* const* Evecs,
               const size_t n1, const size_t n2, const size_t n3)
{
    double dd[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n1]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n1]*STEPSIZE);
        }
        
        dd[i]=get_hess(tmp_coord,ncoord,Evecs,n2,n3);
        
    }
    
    return (dd[1]-dd[0])/(2*STEPSIZE);
    
}

double get3MR(const double *crd, const size_t ncoord,
              double **Evecs,
              const size_t n1, const size_t n2, const size_t n3)
{
    
    double dd[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n1]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n1]*STEPSIZE);
        }
        
        dd[i]=get2MR(tmp_coord,ncoord,Evecs,n2,n3);
        
    }
    
    return (dd[1]-dd[0])/(2*STEPSIZE);
    
}

double get_quar(const double *crd, const size_t ncoord,
                double const* const* Evecs,
                const size_t n1, const size_t n2, const size_t n3, const size_t n4)
{
    double dd[2];
    for (size_t i=0; i<2; ++i) {
        
        double tmp_coord[ncoord];
        
        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n1]*STEPSIZE :
                                    crd[ndiff] - Evecs[ndiff][n1]*STEPSIZE);
        }
        
        dd[i]=get_cub(tmp_coord,ncoord,Evecs,n2,n3,n4);
        
    }
    
    return (dd[1]-dd[0])/(2*STEPSIZE);
    
}

void third_derivative(const double *crd, const size_t ncoord,
                      double **Evecs, double *dddE)
{
    // do numerical derivative for each normal coordinate
    for (size_t n0=0; n0<ncoord; ++n0) {
        for (size_t n1=0; n1<ncoord; ++n1) {
            for (size_t n2=0; n2<ncoord; ++n2) {
            
                double dd[2];
                for (size_t i=0; i<2; ++i) {
                    
                    double tmp_coord[ncoord];
                    
                    for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
                            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n0]*STEPSIZE :
                                                    crd[ndiff] - Evecs[ndiff][n0]*STEPSIZE);

                    }
                    
                    dd[i]=get2MR(tmp_coord,ncoord,Evecs,n1,n2);
                    
                }
                dddE[n0*ncoord+n1*ncoord+n2]=(dd[1]-dd[0])/(2*STEPSIZE);
                dddE[n0*ncoord+n1*ncoord+n2]=dddE[n0*ncoord+n1*ncoord+n2]*constants::Bohr_A*constants::Bohr_A/constants::Eh_kcalmol;
            }
        }
    }
    return;
    
}

void fourth_derivative(const double *crd, const size_t ncoord,
                       double **Evecs, double *ddddE)
{
    
    // do numerical derivative for each normal coordinate
    for (size_t n0=0; n0<ncoord; ++n0) {
        for (size_t n1=0; n1<ncoord; ++n1) {
            for (size_t n2=0; n2<ncoord; ++n2) {
                for (size_t n3=0; n3<ncoord; ++n3) {
                
                    double ddd[2];
                    for (size_t i=0; i<2; ++i) {
                        
                        double tmp_coord[ncoord];
                        
                        for (size_t ndiff = 0; ndiff < ncoord; ++ndiff) {
                            tmp_coord[ndiff] = (i ? crd[ndiff] + Evecs[ndiff][n0]*STEPSIZE :
                                                    crd[ndiff] - Evecs[ndiff][n0]*STEPSIZE);
                            
                        }
                        
                        ddd[i]=get3MR(tmp_coord,ncoord,Evecs,n1,n2,n3);
                        
                    }
                    
                    ddddE[n0*ncoord+n1*ncoord+n2*ncoord+n3]=(ddd[1]-ddd[0])/(2*STEPSIZE);
                    ddddE[n0*ncoord+n1*ncoord+n2*ncoord+n3]=ddddE[n0*ncoord+n1*ncoord+n2*ncoord+n3]*constants::Bohr_A*constants::Bohr_A/constants::Eh_kcalmol;
                }
            }
        }
    }
    return;
    
}

void calculate_freq(const double *ddE, const size_t ncoord,
                    const double *mass, double *freq, 
                    double **Evecs) 
{
  double **H;
  H = new double*[ncoord];
  for (size_t i=0; i<ncoord; ++i) {
    H[i] = new double[ncoord];
  }  

  for (size_t n0=0; n0<ncoord; ++n0) {
    for (size_t n1=0; n1<ncoord; ++n1) {
      
      H[n0][n1]=ddE[n0*ncoord+n1]/sqrt(mass[n0/3])/sqrt(mass[n1/3]);

    }
  }

  dsyev(H,ncoord,freq,Evecs);

  const char fm[] = "frequencies.txt";
  FILE *fs = fopen(fm, "w");
  
  for (size_t i=0; i<ncoord; ++i) {
    if (freq[i]<0) fprintf(fs,"%15.10f cm-1\n",-sqrt(std::abs(freq[i]))*5140.4862);
    else fprintf(fs,"%15.10f cm-1\n",sqrt(std::abs(freq[i]))*5140.4862);
  }

  const char fn[] = "modes.txt";
  FILE *fr = fopen(fn, "w");

  for (size_t i=0; i<ncoord; ++i) { 
    if (freq[i]<0) fprintf(fr,"%7.2f cm-1\n",-sqrt(std::abs(freq[i]))*5140.4862);
    else fprintf(fr,"%7.2f cm-1\n",sqrt(std::abs(freq[i]))*5140.4862);

    for (size_t j=0; j<ncoord; ++j) {
      fprintf(fr,"%15.10f\n",Evecs[j][i]);
    }
    fprintf(fr,"\n");
  }

  return;

}

void mavi_header(const size_t natoms, const size_t ncoord, const double *crd, 
                 const double *mass, const double nrg, const size_t nm)
{
    const char ff[] = "001.hs";
    FILE *fm = fopen(ff, "w");
 
//    size_t nm;
//    if (linear) nm = natoms*3-5
//    else f
//    else nm = natoms*3-6;
    
    fprintf(fm,"# Energy / hartree \n");
    fprintf(fm,"%20.10f\n",nrg);
   
    size_t c = nm/3;
    if (nm%3>0) c+=1;

    fprintf(fm," #Geometry / Angs amu1/2\n");
    for (size_t i=0; i<(c); i++)
	fprintf(fm," %15.10f  %15.10f  %15.10f\n",0.0,0.0,0.0);

    fclose(fm);
    return;
}

void oneMR(const size_t natoms, const size_t ncoord, const double *crd, const double *mass, 
           double **Evecs, double *Gi, double *Hii, double *Ciii, double *Qiiii, const size_t nm,
           const size_t na, int *ac)
{

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    for (size_t i=0; i<ncoord; i++) {
      for (size_t j=0; j<ncoord; j++) {

        Evecs[j][i]=Evecs[j][i]*constants::Bohr_A/sqrt(mass[j/3]);
    
      }
    }

    for (size_t i=0; i<na; i++) {
	size_t n=ncoord-nm+ac[i];
	//size_t n=ac[i];
	//std::cerr << n << " - " << ncoord << " - " << nm << " - " << ac[i] << " : n,ncoord,nm,ac" << std::endl;

        Gi[i]=get_grad(crd,ncoord,Evecs,n);
        Hii[i]=get_hess(crd,ncoord,Evecs,n,n);
        Ciii[i]=get_cub(crd,ncoord,Evecs,n,n,n);
        Qiiii[i]=get_quar(crd,ncoord,Evecs,n,n,n,n);
    }
    return;
}

void twoMR(const size_t natoms, const size_t ncoord, const double *crd, const double *mass,
           double **Evecs, double *Hij, double *Ciij, double *Qiijj, double *Qiiij, const size_t nm,
           const size_t na, int *ac)
{

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    for (size_t i=0; i<ncoord; i++) {
      for (size_t j=0; j<ncoord; j++) {
        Evecs[j][i]=Evecs[j][i];//*constants::Bohr_A/sqrt(mass[j/3]);
      }
    }

    double a=constants::Bohr_A;
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        size_t ni=ncoord-nm+ac[i];
        size_t nj=ncoord-nm+ac[j];
        //size_t ni=ncoord-nm+i;
        //size_t nj=ncoord-nm+j;

        Hij[i*na+j]  =get_hess(crd,ncoord,Evecs,ni,nj);
	Qiijj[i*na+j]=get_quar(crd,ncoord,Evecs,ni,ni,nj,nj);
        Ciij[i*na+j] =get_cub(crd,ncoord,Evecs,ni,ni,nj);
        Ciij[j*na+i] =get_cub(crd,ncoord,Evecs,nj,nj,ni);
        Qiiij[i*na+j]=get_quar(crd,ncoord,Evecs,ni,ni,ni,nj);
        Qiiij[j*na+i]=get_quar(crd,ncoord,Evecs,nj,nj,nj,ni);

      }
    }
    return;
}

void threeMR(const size_t natoms, const size_t ncoord, const double *crd, const double *mass,
             double **Evecs, double *Cijk, double *Qiijk, const size_t nm, const size_t na, int *ac)
{

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    for (size_t i=0; i<ncoord; i++) {
      for (size_t j=0; j<ncoord; j++) {
        Evecs[j][i]=Evecs[j][i];//*constants::Bohr_A/sqrt(mass[j/3]);
      }
    }

    double a=constants::Bohr_A;
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
	for (size_t k=j+1; k<na; k++) {
	  size_t ni=ncoord-nm+ac[i];
          size_t nj=ncoord-nm+ac[j];
          size_t nk=ncoord-nm+ac[k];
          //size_t ni=ncoord-nm+i;
          //size_t nj=ncoord-nm+j;
	  //size_t nk=ncoord-nm+k;

          Cijk[i*na*na+j*na+k] =get_cub(crd,ncoord,Evecs,ni,nj,nk);
          Qiijk[i*na*na+j*na+k]=get_quar(crd,ncoord,Evecs,ni,ni,nj,nk);
          Qiijk[j*na*na+k*na+i]=get_quar(crd,ncoord,Evecs,nj,nj,nk,ni);
	  Qiijk[k*na*na+i*na+j]=get_quar(crd,ncoord,Evecs,nk,nk,ni,nj);
	}
      }
    }
    return;
}

void fourMR(const size_t natoms, const size_t ncoord, const double *crd, const double *mass,
            double **Evecs, double *Qijkl, const size_t nm, const size_t na, int *ac)
{

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    for (size_t i=0; i<ncoord; i++) {
      for (size_t j=0; j<ncoord; j++) {
        Evecs[j][i]=Evecs[j][i];//*constants::Bohr_A/sqrt(mass[j/3]);
      }
    }

    double a=constants::Bohr_A;
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        for (size_t k=j+1; k<na; k++) {
	  for (size_t l=k+1; l<na; l++) {
            size_t ni=ncoord-nm+ac[i];
	    size_t nj=ncoord-nm+ac[j];
	    size_t nk=ncoord-nm+ac[k];
            size_t nl=ncoord-nm+ac[l];
            //size_t ni=ncoord-nm+i;
            //size_t nj=ncoord-nm+j;
            //size_t nk=ncoord-nm+k;
            //size_t nl=ncoord-nm+l;

            Qijkl[i*na*na*na+j*na*na+k*na+l]=get_quar(crd,ncoord,Evecs,ni,nj,nk,nl);
	  }
        }
      }
    }
    return;
}
void mavi_1MR(const size_t natoms,const size_t ncoord,
              const double *tmp, const double *mass,
              double **Evecs, const size_t nm, const size_t na, int *ac)
{
    const char ff[] = "001.hs";
    FILE *fm = fopen(ff, "a");
    
//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    double Gi[na];
    double Hii[na];
    double Ciii[na];
    double Qiiii[na];
    
    oneMR(natoms,ncoord,tmp,mass,Evecs,Gi,Hii,Ciii,Qiiii,nm,na,ac);
    
    fprintf(fm,"# 1MR\n");
    fprintf(fm,"# Gradient / hartree Angs^-1 amu^-1/2 \n");
    for (size_t i=0; i<na; i++) {
      double a=constants::Bohr_A;
      fprintf(fm,"%4zu  %15.10f\n",i+1,Gi[i]/constants::Eh_kcalmol/a);
    }
    fprintf(fm,"# Hessian(i,i) / hartree Angs^-2 amu^-1 \n");
    for (size_t i=0; i<na; i++) {
      double a=constants::Bohr_A;
      fprintf(fm,"%4zu  %15.10f\n",i+1,Hii[i]/constants::Eh_kcalmol/a/a);
    }
    fprintf(fm,"# Cubic(i,i,i) / hartree Angs^-3 amu^-3/2 \n");
    for (size_t i=0; i<na; i++) {
      double a=constants::Bohr_A;
      fprintf(fm,"%4zu  %15.10f\n",i+1,Ciii[i]/constants::Eh_kcalmol/a/a/a);
    }
    fprintf(fm,"# Quartic(i,i,i,i) / hartree Angs^-4 amu^-2 \n");
    for (size_t i=0; i<na; i++) {
      double a=constants::Bohr_A;
      fprintf(fm,"%4zu  %15.10f\n",i+1,Qiiii[i]/constants::Eh_kcalmol/a/a/a/a);
    }
    fclose(fm);
    return;
}

void mavi_2MR(const size_t natoms,const size_t ncoord,
              const double *tmp, const double *mass,
              double **Evecs, const size_t nm, const size_t na, int *ac)
{

    const char ff[] = "001.hs";
    FILE *fm = fopen(ff, "a");

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    double Hij[na*na];
    double Qiijj[na*na];
    double Ciij[na*na];
    double Qiiij[na*na];

    for (size_t i=0; i<na; i++) {
      for (size_t j=0; j<na; j++) {
        Hij[i*na+j]=0.0;
        Qiijj[i*na+j]=0.0;
        Ciij[i*na+j]=0.0;
        Qiiij[i*na+j]=0.0;
      }
    }

    twoMR(natoms,ncoord,tmp,mass,Evecs,Hij,Ciij,Qiijj,Qiiij,nm,na,ac);

    fprintf(fm,"# 2MR\n");
    fprintf(fm,"# Hessian(i,j) / hartree Angs^-2 amu^-1 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        double a=constants::Bohr_A;
        fprintf(fm," %3zu %3zu %15.10f\n",j+1,i+1,Hij[i*na+j]/constants::Eh_kcalmol/a/a);
      }
    }
    fprintf(fm,"# Quartic(i,i,j,j) / hartree Angs^-4 amu^-2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
	double a=constants::Bohr_A;
        fprintf(fm," %3zu %3zu %15.10f\n",j+1,i+1,Qiijj[i*na+j]/constants::Eh_kcalmol/a/a/a/a);
      }
    }
    fprintf(fm,"# Cubic(i,i,j) / hartree Angs^-3 amu^-3/2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        double a=constants::Bohr_A;
        fprintf(fm," %3zu %3zu %15.10f\n",j+1,i+1,Ciij[j*na+i]/constants::Eh_kcalmol/a/a/a);
        fprintf(fm," %3zu %3zu %15.10f\n",i+1,j+1,Ciij[i*na+j]/constants::Eh_kcalmol/a/a/a);
      }
    }
    fprintf(fm,"# Quartic(i,i,i,j) / hartree Angs^-4 amu^-2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        double a=constants::Bohr_A;
        fprintf(fm," %3zu %3zu %15.10f\n",j+1,i+1,Qiiij[j*na+i]/constants::Eh_kcalmol/a/a/a/a);
        fprintf(fm," %3zu %3zu %15.10f\n",i+1,j+1,Qiiij[i*na+j]/constants::Eh_kcalmol/a/a/a/a);
      }
    }
    fclose(fm);
    return;
}

void mavi_3MR(const size_t natoms,const size_t ncoord,
              const double *tmp, const double *mass,
              double **Evecs, const size_t nm, const size_t na, int *ac)
{

    const char ff[] = "001.hs";
    FILE *fm = fopen(ff, "a");

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    double Cijk[na*na*na];
    double Qiijk[na*na*na];

    for (size_t i=0; i<na; i++) {
      for (size_t j=0; j<na; j++) {
	for (size_t k=0; k<na; k++) {
          Cijk[i*na*na+j*na+k]=0.0;
          Qiijk[i*na*na+j*na+k]=0.0;
	}
      }
    }

    threeMR(natoms,ncoord,tmp,mass,Evecs,Cijk,Qiijk,nm,na,ac);

    fprintf(fm,"# 3MR\n");
    fprintf(fm,"# Cubic(i,j,k) / hartree Angs^-3 amu^-3/2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
	for (size_t k=j+1; k<na; k++) {
          double a=constants::Bohr_A;
          fprintf(fm," %3zu %3zu %3zu %15.10f\n",k+1,j+1,i+1,Cijk[i*na*na+j*na+k]/constants::Eh_kcalmol/a/a/a);
        }
      }
    }
    fprintf(fm,"# Quartic(i,i,j,k) / hartree Angs^-4 amu^-2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
	for (size_t k=j+1; k<na; k++) {
          double a=constants::Bohr_A;
          fprintf(fm," %3zu %3zu %3zu %15.10f\n",k+1,j+1,i+1,Qiijk[i*na*na+j*na+k]/constants::Eh_kcalmol/a/a/a/a);
          fprintf(fm," %3zu %3zu %3zu %15.10f\n",j+1,i+1,k+1,Qiijk[k*na*na+i*na+j]/constants::Eh_kcalmol/a/a/a/a);
          fprintf(fm," %3zu %3zu %3zu %15.10f\n",i+1,k+1,j+1,Qiijk[j*na*na+k*na+i]/constants::Eh_kcalmol/a/a/a/a);
	}
      }
    }
    fclose(fm);
    return;
}

void mavi_4MR(const size_t natoms,const size_t ncoord,
              const double *tmp, const double *mass,
              double **Evecs, const size_t nm, const size_t na, int *ac)
{

    const char ff[] = "001.hs";
    FILE *fm = fopen(ff, "a");

//    size_t nm;
//    if (natoms==3) nm = 4;
//    else nm = natoms*3-6;

    double Qijkl[na*na*na*na];

    for (size_t i=0; i<na; i++) {
      for (size_t j=0; j<na; j++) {
        for (size_t k=0; k<na; k++) {
	  for (size_t l=0; l<na; l++) {
            Qijkl[i*na*na*na+j*na*na+k*na+l]=0.0;
          }
        }
      }
    }

    fourMR(natoms,ncoord,tmp,mass,Evecs,Qijkl,nm,na,ac);

    fprintf(fm,"# 4MR\n");
    fprintf(fm,"# Quartic(i,j,k,l) / hartree Angs^-4 amu^-2 \n");
    for (size_t i=0; i<na; i++) {
      for (size_t j=i+1; j<na; j++) {
        for (size_t k=j+1; k<na; k++) {
	  for (size_t l=k+1; l<na; l++) {
            double a=constants::Bohr_A;
            fprintf(fm," %3zu %3zu %3zu %3zu %15.10f\n",l+1,k+1,j+1,i+1,Qijkl[i*na*na*na+j*na*na+k*na+l]/constants::Eh_kcalmol/a/a/a/a);
	  }
        }
      }
    }
    fclose(fm);
    return;
}

void read_mass(const size_t ncoord, double *mass)
{

  const char* filename = "mass";

  std::ifstream ifs(filename);
  if (!ifs) {
    std::ostringstream oss;
    oss << "could not open 'mass' for reading";
    throw std::runtime_error(oss.str());
  }

  std::string line;

  if (ifs.eof())
      return;

  for (size_t i=0; i<ncoord; i++) {
    std::getline(ifs, line);
    if (ifs.eof())
      return;

    std::istringstream iss(line);
    iss >> mass[i];
  }

  fprintf(stderr,"mass file read.\n");
  return;

}

void read_atomic(const size_t ncoord, double *z)
{

  const char* filename = "atomic";

  std::ifstream ifs(filename);
  if (!ifs) {
    std::ostringstream oss;
    oss << "could not open 'atomic' for reading";
    throw std::runtime_error(oss.str());
  }

  std::string line;

  if (ifs.eof())
      return;

  for (size_t i=0; i<ncoord; i++) {
    std::getline(ifs, line);
    if (ifs.eof())
      return;

    std::istringstream iss(line);
    iss >> z[i];
  }
  fprintf(stderr,"atomic charge read.\n");
  return;

}


void read_second_derivative(const size_t ncoord, double *ddE) 
{

  const char* filename = "hess";
    
  std::ifstream ifs(filename);
  if (!ifs) {
    std::ostringstream oss;
    oss << "could not open 'hess' for reading";
    throw std::runtime_error(oss.str());
  }

  std::string line;

  if (ifs.eof())
      return;

  for (size_t i=0; i<ncoord*ncoord; i++) {
    std::getline(ifs, line);
    if (ifs.eof())
      return;

    std::istringstream iss(line);
    iss >> ddE[i];
  }
  fprintf(stderr,"Hessians read.\n");
  return;

}

int main(int argc, char** argv)
{

  const char* filename;
  size_t linear = 0;
  //size_t reduced = 0;
  size_t fixed = 0;
  size_t nmodes = 0;
  int *temp;
  int *active = 0;
  size_t nactive = 0;

  for(int i = 0; i < argc; i++) {
    if (std::string(argv[i]) == "-i") {
      filename = argv[i + 1];
    } else if (std::string(argv[i]) == "-l") {
      linear = 1;
    } else if (std::string(argv[i]) == "-n") {
      nmodes = atoi(argv[i+1]);
      
    } else if (std::string(argv[i]) == "-a") {
      //reduced = 1;
      temp = new int[nmodes];
      for (int j = 0; j < nmodes; j++) { 
        temp[j] = atoi(argv[i+j+1]);
        if (temp[j]==1) nactive++;
      }
      active = new int[nactive];
      nactive = 0;
      for (int j = 0; j < nmodes; j++) { 
        if (temp[j]==1) { active[nactive] = j; nactive++; }
       
      //fixed = atoi(argv[i+1]);
      }
    } else if (std::string(argv[i]) == "-h") {
      std::cerr << std::endl; 
      std::cerr << "Usage: anharmonic.x -i input-file [-l] [-m int] [-a (int int int ...)] [-h]" << std::endl;
      std::cerr << std::endl; 
      std::cerr << "Options/Arguments:" << std::endl;
      std::cerr << "    -i char             input file to read [required]." << std::endl;
      std::cerr << "    -l                  linear molecule." << std::endl;
      std::cerr << "    -n int              number of degrees of freedom. (required for specifying active modes)." << std::endl; 
      std::cerr << "    -a int int ...      activity of all mode vibrations." << std::endl;
      std::cerr << "    -h                  display help." << std::endl;
      std::cerr << std::endl;
      exit(0);
    }
  }

#ifdef _OPENMP   
  int threads = omp_get_num_threads();
//  int threads = 1;
  omp_set_num_threads(threads);
  fprintf(stderr,"Parallel execution on %d allocated threads.\n", threads);
#endif

  fprintf(stderr,"\t***********************************************************\n");
  fprintf(stderr,"\t**             Numerical Cluster Frequencies             **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t**                  Principal Author:                    **\n");
  fprintf(stderr,"\t**                  Olaseni Sode                         **\n");
  fprintf(stderr,"\t**                  email: osode@calstatela.edu          **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t***********************************************************\n");

  std::ifstream ifs(filename);
  if (!ifs) {
    std::ostringstream oss;
    oss << "could not open '" << *argv << "' for reading";
    throw std::runtime_error(oss.str());
  }
 
  size_t natoms=0;
  std::string comment;
  std::vector<std::string> elements;
  std::vector<double> xyz;
 
  /* Load starting structure */
  io::load_xyz(ifs, natoms, comment, elements, xyz);
    
//  ++argv;
//  std::ofstream ofs(*argv);

  /* Initialize coordinates */
  size_t ncoord = natoms*3;
  double tmp[ncoord];
  double ddE[ncoord*ncoord];
  double dddE[ncoord*ncoord*ncoord];
  double ddddE[ncoord*ncoord*ncoord*ncoord];
  double mass[natoms];
  double z[natoms];

  for (size_t n0=0; n0<ncoord; ++n0) {
    tmp[n0]=xyz[n0];
    for (size_t n1=0; n1<ncoord; ++n1) {
      ddE[n0*ncoord+n1]=0.0;
      for (size_t n2=0; n2<ncoord; ++n2) {
        dddE[n0*ncoord+n1*ncoord+n2]=0.0;
        for (size_t n3=0; n3<ncoord; ++n3) {
          dddE[n0*ncoord+n1*ncoord+n2*ncoord+n3]=0.0;
        }
      }
    }
  }

  read_mass(natoms,mass);
  read_atomic(natoms,z);
  double nrg=get_energy(tmp,ncoord)/627.509;

  second_derivative(tmp,ncoord,ddE);

  double freq[ncoord];

  double **Evecs;
  Evecs = new double*[ncoord];
  for (size_t i=0; i<ncoord; ++i) {
    freq[i]=0.0;
    Evecs[i] = new double[ncoord];
    for (size_t j=0; j<ncoord; ++j) {
      Evecs[i][j]=0.0;
    }
  }
  
  calculate_freq(ddE,ncoord,mass,freq,Evecs);

  if (nactive>0) ;
  else if (linear) {
    nactive = natoms*3-5;
    nmodes = nactive;
    active = new int[nactive];
    for (size_t i=0; i<nactive; i++) active[i] = i;
  }
  else {
     nactive = natoms*3-6;
     nmodes = nactive;
     active = new int[nactive];
     for (size_t i=0; i<nactive; i++) active[i] = i;
  }

  print_molpro_data(natoms,ddE,ncoord,mass,z,tmp,freq,Evecs);
  fprintf(stderr,"Molpro output file 'molpro.out' printed.\n");
  print_nwchem_data(natoms,ddE,ncoord,mass,z,tmp,freq,Evecs);
  fprintf(stderr,"NWChem output file 'nwchem.out' printed.\n");
  print_nmode(natoms,Evecs,ncoord,mass,freq,nmodes,nactive,active);
  fprintf(stderr,"Normal mode file 'nmode' printed.\n");
  print_sindo(natoms,ddE,ncoord,mass,z,tmp,freq,Evecs,nmodes,nactive,active);
  fprintf(stderr,"SINDO input file 'sindo.inp' printed.\n");
  print_nitrogen(natoms,tmp);
  fprintf(stderr,"NITROGEN ref coordinates 'nitrogen.dat' printed.\n");
  
  mavi_header(natoms,ncoord,tmp,mass,nrg,nactive);
  mavi_1MR(natoms,ncoord,tmp,mass,Evecs,nmodes,nactive,active);
  mavi_2MR(natoms,ncoord,tmp,mass,Evecs,nmodes,nactive,active);
  mavi_3MR(natoms,ncoord,tmp,mass,Evecs,nmodes,nactive,active);
  mavi_4MR(natoms,ncoord,tmp,mass,Evecs,nmodes,nactive,active);
  fprintf(stderr,"MaVi input file '001.hs' printed.\n");

  return 0;
}

