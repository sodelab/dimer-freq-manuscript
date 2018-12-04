#include <sstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////

void print_molpro_data(const size_t natoms, const double *ddE,
                       const size_t ncoord, const double *mass,
                       const double *z,     const double *crd,
                       const double *freq,  double **Evecs);

void print_nwchem_data(const size_t natoms, const double *ddE,
                       const size_t ncoord, const double *mass,
                       const double *z,     const double *crd,
                       const double *freq,  double **Evecs);

void print_nmode(const size_t natoms, double **Evecs,
                 const size_t ncoord, const double *mass,
                 double *freq,         const size_t nmodes,
                 const size_t nac,     const int *ac);

void print_sindo(const size_t natoms, const double *ddE,
                 const size_t ncoord, const double *mass,
                 const double *z,     const double *crd,
                 const double *freq,  double **Evecs,
                 const size_t nmodes, const size_t nac,
                 const int *ac);

void print_nitrogen(const size_t natoms, const double *crd);
////////////////////////////////////////////////////////////////////////////////
