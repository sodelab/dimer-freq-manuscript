#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "x1b.h"
#include "x2b.h"
//#include "kit.h"

const double Eh_J = 4.35974434e-18; // CODATA 2010
const double Na = 6.02214129e+23; // CODATA 2010

const double kcal_J = 4184.0;
const double Eh_kcalmol = Eh_J*Na/kcal_J;

#define Eh_J 4.35974434e-18
#define Na 6.02214129e+23
#define kcal_J 4184.0
#define Eh_kcalmol Eh_J*Na/kcal_J
#define Nq 6                            //Number of internal coordinates
#define N 6				//Number of atoms

//
// getN() is a required exported function
// Returns N, the number of atoms
#ifdef _WIN32
__declspec (dllexport)
#endif
extern "C" int getN(){
        return N;
}


// calcSurface(...) is a required exported function
// Returns the value of the PES given the nuclear
// cartesian coordinates (X). The expected atom ordering
// is O C O O C O
#ifdef _WIN32
__declspec (dllexport)
#endif
extern "C" double calcSurface(double *X, int *PERM){

	double mon1[N/2*3];
 	std::copy(X + 0, X + 3, mon1 + 0);
	std::copy(X + 3, X + 6, mon1 + 3);
	std::copy(X + 6, X + 9, mon1 + 6);

        double mon2[N/2*3];
        std::copy(X +  9, X + 12, mon2 + 0);
        std::copy(X + 12, X + 15, mon2 + 3);
        std::copy(X + 15, X + 18, mon2 + 6);

	double dim[N*3];
        std::copy(mon1 + 0, mon1 +  9, dim +  0);
        std::copy(mon2 + 9, mon2 + 18, dim +  9);

	double value = 0.0;
        value += x1b::value(mon1);
        value += x1b::value(mon2);

	value += x2b::value(X);

        // Return value of PES in wavenumbers
	return value*349.75;
}

double p2b(double *X){

        double mon1[3*3];
        std::copy(X + 0, X + 3, mon1 + 0);
        std::copy(X + 3, X + 6, mon1 + 3);
        std::copy(X + 6, X + 9, mon1 + 6);
        
        double mon2[3*3];
        std::copy(X +  9, X + 12, mon2 + 0);
        std::copy(X + 12, X + 15, mon2 + 3);
        std::copy(X + 15, X + 18, mon2 + 6);

        double dim[6*3];
        std::copy(mon1 + 0, mon1 +  9, dim +  0);
        std::copy(mon2 + 9, mon2 + 18, dim +  9);

        double value = 0.0;
//        value += x1b::value(mon1);
//        value += x1b::value(mon2);

        value += x2b::value(X);

        // Return value of PES in kcal/mol
        return value;
}

double p1b(double *X){

        double mon[3*3];
        std::copy(X + 0, X + 3, mon + 0);
        std::copy(X + 3, X + 6, mon + 3);
        std::copy(X + 6, X + 9, mon + 6);

        double value = 0.0;
        value += x1b::value(mon);
  
        // Return value of PES in kcal/mol
        return value;
}


