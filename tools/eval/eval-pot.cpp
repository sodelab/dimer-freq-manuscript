#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"
#include "mbCO2.h"

int main(int argc, char** argv)
{
    if (argc != 2) {
	std::cerr << "usage: eval monomer.xyz" << std::endl;
	return 0;
    }

    size_t natoms=0;
    std::string comment;
    std::vector<std::string> elements;
    std::vector<double> xyz;
    const char* filename;

    ++argv; --argc;

    std::ifstream ifs(*argv);
    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << *argv << "' for reading";
        throw std::runtime_error(oss.str());
    }

    /* Load starting structure */
    io::load_xyz(ifs, natoms, comment, elements, xyz);

   double tmp[natoms*3];
   for (size_t n=0; n<natoms*3; ++n) tmp[n]=xyz[n];

    std::cout << std::scientific << std::setprecision(9);
    double E = 0.0;
    double nfrags=natoms/3;

#pragma omp parallel for schedule(dynamic) reduction(+:E) 
    for(size_t ifrag=0; ifrag<nfrags; ifrag++){
        E += p1b(tmp + 9*ifrag);
    }

#pragma omp parallel for schedule(dynamic) reduction(+:E) 
    for(size_t ifrag=0; ifrag<nfrags; ifrag++){
        for(size_t jfrag=ifrag+1; jfrag<nfrags; jfrag++){

            double dim[18];
            std::copy(tmp + ifrag*9, tmp + ifrag*9 + 9, dim + 0);
            std::copy(tmp + jfrag*9, tmp + jfrag*9 + 9, dim + 9);

            double nnn = p2b(dim);
            E += nnn;

	}
    }

    std::cout << E/627.509 << std::endl;

    return 0;
}


