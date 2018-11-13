#include <cmath>
#include <cassert>
#include <cstdlib> 

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"
#include "mbCO2.h"

#include <gsl/gsl_multimin.h>

#ifdef _OPENMP
#include "omp.h"
#endif

#define STEPSIZE 0.001

void save_frame(size_t t, const std::vector<double>& pos, int natoms, char* comment)
{

  std::cout << std::setw(2) << (natoms) << '\n'
            << std::setw(4) << t << ' ' << comment << '\n';

  for (int n = 0; n < natoms; ++n)
    std::cout << (n%3 == 0 ? "C " : "O ")
              << std::setw(20) << pos[3*n + 0]
              << std::setw(20) << pos[3*n + 1]
              << std::setw(20) << pos[3*n + 2]
              << '\n';
}

void save_frame(size_t t, double* pos, int natoms, char* comment)
{

  std::cout << std::setw(2) << (natoms) << '\n'
            << std::setw(4) << t << ' ' << comment << '\n';

  for (int n = 0; n < natoms; ++n)
    std::cout << (n%3 == 0 ? "C " : "O ")
              << std::setw(20) << pos[3*n + 0]
              << std::setw(20) << pos[3*n + 1]
              << std::setw(20) << pos[3*n + 2]
              << '\n';
}

double get_energy(const double *crd, const size_t ncoord) 
{
  size_t natoms = ncoord/3;
  size_t nfrags = natoms/3;

  std::cout << std::scientific << std::setprecision(9);

  double pos[ncoord];

  for (size_t n = 0; n < ncoord; ++n)
    pos[n] = crd[n];

  double E = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:E) 
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    E += p1b(pos + 9*ifrag);
  }
#else
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    E += p1b(pos + 9*ifrag);
  }
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:E) 
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    for(size_t jfrag=ifrag+1; jfrag<nfrags; jfrag++){

      double dim[18];
      std::copy(pos + ifrag*9, pos + ifrag*9 + 9, dim + 0);
      std::copy(pos + jfrag*9, pos + jfrag*9 + 9, dim + 9);

      double nnn = p2b(dim);
      E += nnn;

      if (nnn < -2.0) {
        double x = dim[9+0] - dim[0+0];
        double y = dim[9+1] - dim[0+1];
        double z = dim[9+2] - dim[0+2];
        double d = sqrt ( x*x + y*y + z*z ); // C-C distance

   	std::cerr << ifrag << " " << jfrag << " " << d << " " << nnn << std::endl;
        for (size_t n = 0; n < 6; ++n)
          std::cout << (n%3 == 0 ? "C " : "O ")
                    << std::setw(20) << dim[3*n + 0]
                    << std::setw(20) << dim[3*n + 1]
                    << std::setw(20) << dim[3*n + 2]
                    << '\n';
	exit(1);
      }
    }
  }
#else
  for(size_t ifrag=0; ifrag<nfrags; ifrag++){
    for(size_t jfrag=ifrag+1; jfrag<nfrags; jfrag++){

      double dim[18];
      std::copy(pos + ifrag*9, pos + ifrag*9 + 9, dim + 0);
      std::copy(pos + jfrag*9, pos + jfrag*9 + 9, dim + 9);

      double nnn = p2b(dim);
      E += nnn;

      if (nnn < -2.0) {
        double x = dim[9+0] - dim[0+0];
        double y = dim[9+1] - dim[0+1];
        double z = dim[9+2] - dim[0+2];
        double d = sqrt ( x*x + y*y + z*z ); // C-C distance

        std::cerr << ifrag << " " << jfrag << " " << d << " " << nnn << std::endl;
        for (size_t n = 0; n < 6; ++n)
          std::cout << (n%3 == 0 ? "C " : "O ")
                    << std::setw(20) << dim[3*n + 0]
                    << std::setw(20) << dim[3*n + 1]
                    << std::setw(20) << dim[3*n + 2]
                    << '\n';
        exit(1);
      }
    }
  }
#endif

  return E;
}

double cluster_f (const gsl_vector *v, void *params)
{
    size_t ncoord = v->size; // number of coordinates in optimization space
    
    double crd[ncoord];
    
    for (size_t n = 0; n < ncoord; ++n)
        crd[n] = gsl_vector_get(v, n);
	
    double E = get_energy(crd,ncoord);
    std::cerr << "E: " << std::setw(20) << E << std::endl;

    return E;
}

void cluster_df (const gsl_vector *v, void *params,
             gsl_vector *df)
{
  size_t ncoord = v->size; // number of coordinates in optimization space

  double dE[ncoord];
  double crd[ncoord];

  for (size_t n = 0; n < ncoord; ++n) {
    crd[n] = gsl_vector_get(v, n);
  }

  // do numerical derivative for each coordinate
  for (size_t n=0; n<ncoord; ++n) {

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
    dE[n]=(e_num[1]-e_num[0])/(2*STEPSIZE); 
  }

  for (size_t n=0; n<ncoord; ++n) {
    gsl_vector_set (df, n, dE[n]);
  }
  return;
}

void cluster_fdf (const gsl_vector *v, void *params,
              double *f, gsl_vector *df)
{

  *f = cluster_f(v, params);
  cluster_df(v, params, df);

  return;
}

int main(int argc, char** argv)
{

  if (argc < 2) {
    std::cerr << "usage: opt-cluster.x initial_file.xyz  > output_file.xyz"
              << std::endl;
    return 0;
  }

#ifdef _OPENMP   
//  int threads = omp_get_num_threads();
  int threads = 1;
  omp_set_num_threads(threads);
  fprintf(stderr,"Parallel execution on %d allocated threads.\n", threads);
#endif

  fprintf(stderr,"\t***********************************************************\n");
  fprintf(stderr,"\t**             Numerical Cluster Optimization            **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t**                  Principal Author:                    **\n");
  fprintf(stderr,"\t**                  Olaseni Sode                         **\n");
  fprintf(stderr,"\t**                  email: osode@ut.edu                  **\n");
  fprintf(stderr,"\t**                                                       **\n");
  fprintf(stderr,"\t***********************************************************\n");

  ++argv; --argc; 
  const char* filename = *argv;

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
  kit::io::load_xyz(ifs, natoms, comment, elements, xyz);
    
  ++argv;
  std::ofstream ofs(*argv);

  /* Initialize coordinates */
  size_t ncoord = natoms*3;
  double tmp[ncoord];

  /* Initialize function minimizer */
  size_t iter = 0;
  int status;

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  double par[1] = { 0.0 };

  gsl_vector *x;
  gsl_multimin_function_fdf minex_func;

  minex_func.n      = ncoord;
  minex_func.f      = cluster_f;
  minex_func.df     = cluster_df;
  minex_func.fdf    = cluster_fdf;
  minex_func.params = par;

  x = gsl_vector_alloc (ncoord);
  for (size_t n=0; n<ncoord; ++n) {
    gsl_vector_set (x, n, xyz[n]);
  }

//  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  T = gsl_multimin_fdfminimizer_vector_bfgs;
  s = gsl_multimin_fdfminimizer_alloc (T, ncoord);

  gsl_multimin_fdfminimizer_set (s, &minex_func, x, 0.002, 1e-2);

  double E = cluster_f(s->x, par);

  char comm[64];
  sprintf(comm,"%15.10f",E);
  save_frame(iter,xyz,natoms,comm);

  do
    {

      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status) 
        break;

      status = gsl_multimin_test_gradient (s->gradient, 2e-3);

      if (status == GSL_SUCCESS)
      {
        std::cerr<< "*** Converged ***" << std::endl;
      }
 
      for (size_t n = 0; n < ncoord; ++n)
        tmp[n] = gsl_vector_get(s->x, n);
      
      {
        char comment[128];
        sprintf(comment,"%15.10f %15.10f",s->f, s->gradient);   
   
        std::string comm(comment);
        save_frame(iter,tmp,natoms,comment);       
      } 
    }
  while (status == GSL_CONTINUE && iter < 5000);

  gsl_vector_free (x);
  gsl_multimin_fdfminimizer_free (s);

  return 0;
}

