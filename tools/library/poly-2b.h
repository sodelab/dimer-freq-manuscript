#ifndef POLY_2B_H
#define POLY_2B_H

//
// reduced (no bonds breaking) permutational symmetry
// including only 2B terms
//
//    x[0] = @VAR@-|x-intra-CO|(Ca, Oa1);
//    x[1] = @VAR@-|x-intra-CO|(Ca, Oa2);
//    x[2] = @VAR@-|x-intra-OO|(Oa1, Oa2);
//    x[3] = @VAR@-|x-intra-CO|(Cb, Ob1);
//    x[4] = @VAR@-|x-intra-CO|(Cb, Ob2);
//    x[5] = @VAR@-|x-intra-OO|(Ob1, Ob2);
//    x[6] = @VAR@-|x-CC|(Ca, Cb);
//    x[7] = @VAR@-|x-CO|(Ca, Ob1);
//    x[8] = @VAR@-|x-CO|(Ca, Ob2);
//    x[9] = @VAR@-|x-CO|(Cb, Oa1);
//    x[10] = @VAR@-|x-CO|(Cb, Oa2);
//    x[11] = @VAR@-|x-OO|(Oa1, Ob1);
//    x[12] = @VAR@-|x-OO|(Oa1, Ob2);
//    x[13] = @VAR@-|x-OO|(Oa2, Ob1);
//    x[14] = @VAR@-|x-OO|(Oa2, Ob2);
//
namespace x2b {

void eval(const double x[15], double*);

}

#endif // POLY_2B_H
