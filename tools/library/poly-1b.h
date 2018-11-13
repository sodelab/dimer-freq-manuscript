#ifndef POLY_1B_H
#define POLY_1B_H

//
// reduced (no bonds breaking) permutational symmetry
// including only 1B terms
//
//    x[0] = @VAR@-|x-intra-CO|(Ca, Oa1);
//    x[1] = @VAR@-|x-intra-CO|(Ca, Oa2);
//    x[2] = @VAR@-|x-intra-OO|(Oa1, Oa2);
//
//
namespace x1b {

void eval(const double x[3], double*);

}

#endif // POLY_1B_H
