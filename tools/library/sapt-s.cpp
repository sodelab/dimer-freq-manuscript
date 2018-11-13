//#ifndef X2B_H
//#define X2B_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>

#include "x2b.h"
#include "kit.h"

////////////////////////////////////////////////////////////////////////////////

namespace x2b {

inline void site(const double* Ca, const double* Oa, double *Ba)
{
    const double dx = Oa[0] - Ca[0];
    const double dy = Oa[1] - Ca[1];
    const double dz = Oa[2] - Ca[2];
    
    double *vector = new double [3];
    
    vector[0]=dx;
    vector[1]=dy;
    vector[2]=dz;
    
    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);
    
    vector[0] = vector[0]/r;
    vector[1] = vector[1]/r;
    vector[2] = vector[2]/r;
    
    const double msite = 0.8456;

    Ba[0] = vector[0]*msite+Ca[0];
    Ba[1] = vector[1]*msite+Ca[1];
    Ba[2] = vector[2]*msite+Ca[2];
    
    delete vector;
    
    return;
}

inline double exponential(const double* A, const double* B,
                          const double &alpha, const double &beta)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    return std::exp(alpha - beta*r);

}
   
inline double qq(const double* A, const double* B, const double &delta,
                 const double &qA, const double &qB)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    
    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);
    const double tt = kit::tang_toennies(1, delta*r);   

    return tt*qA*qB/r;
    
}

inline double xC6(const double* A, const double* B,
                  const double &delta, const double &C6)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    
    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);
    const double tt = kit::tang_toennies(6, delta*r);

    return tt*C6/std::pow(r,6);
}
    
inline double xC8(const double* A, const double* B,
                  const double &delta, const double &C8)
{
    const double dx = A[0] - B[0];
    const double dy = A[1] - B[1];
    const double dz = A[2] - B[2];
    
    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);
    const double tt = kit::tang_toennies(8, delta*r);

    return tt*C8/std::pow(r,8);
}
    
double sapt_s(const double* crd) 
{
//Bukowski, R., Sadlej, J., Jeziorski, B., Jankowski, P., 
//Szalewicz, K., Kucharski, S. A., et al. (1999).
//Intermolecular potential of carbon dioxide dimer from
//symmetry-adapted perturbation theory. 
//The Journal of Chemical Physics, 110(8), 3785–3803. 
//http://doi.org/10.1063/1.479108
    
    const double alpha_OO = 1.1210441e1;
    const double alpha_OC = 1.1333682e1;
    const double alpha_CC = 1.1399839e1;
    
    const double beta_OO  = 4.0202795;
    const double beta_OC  = 4.5074401;
    const double beta_CC  = 5.0932632;
    
    const double qO =  2.3786535e-1*18.22262373;
    const double qC =  1.6316722*18.22262373;
    const double qM = -1.0537015*18.22262373;
    
    const double delta_1_OO = 1.4968485;
    const double delta_1_OC = 1.8797629;
    const double delta_1_CC = 2.1958809;
    const double delta_1_OM = 1.9648279;
    const double delta_1_CM = 2.6032461;
    const double delta_1_MM = 5.2350982;
    
    const double delta_6_OO = 2.5924278;
    const double delta_6_OC = 1.8139783;
    const double delta_6_CC = 1.7584847;
    
    const double delta_8_OO = 1.0769328;
    const double delta_8_OC = 6.7472291e-1;
    const double delta_8_CC = 3.0176726;
    
    const double C6_OO =  1.0426642e3;
    const double C6_OC = -1.3834479e3;
    const double C6_CC =  3.4808543e3;
    
    const double C8_OO = -1.3516797e4;
    const double C8_OC =  2.0217414e4;
    const double C8_CC = -2.6899552e4;
    
    const double* Ca  = crd + 0;
    const double* Oa1 = crd + 3;
    const double* Oa2 = crd + 6;

    const double* Cb  = crd + 9;
    const double* Ob1 = crd + 12;
    const double* Ob2 = crd + 15;

    double *Ba1 = new double [3];
    double *Ba2 = new double [3];
    double *Bb1 = new double [3];
    double *Bb2 = new double [3];
    
    site(Ca,Oa1,Ba1);
    site(Ca,Oa2,Ba2);
    site(Cb,Ob1,Bb1);
    site(Cb,Ob2,Bb2);

    // exp
    const double OO_exp =
            exponential(Oa1,Ob1,alpha_OO,beta_OO)
          + exponential(Oa1,Ob2,alpha_OO,beta_OO)
          + exponential(Oa2,Ob1,alpha_OO,beta_OO)
          + exponential(Oa2,Ob2,alpha_OO,beta_OO);
    
    const double OC_exp =
            exponential(Oa1,Cb,alpha_OC,beta_OC)
          + exponential(Oa2,Cb,alpha_OC,beta_OC)
          + exponential(Ob1,Ca,alpha_OC,beta_OC)
          + exponential(Ob2,Ca,alpha_OC,beta_OC);

    const double CC_exp =
            exponential(Ca,Cb,alpha_CC,beta_CC);
    
    // charge-charge
    const double OO_QQ =
            qq(Oa1,Ob1,delta_1_OO,qO,qO)
          + qq(Oa1,Ob2,delta_1_OO,qO,qO)
          + qq(Oa2,Ob1,delta_1_OO,qO,qO)
          + qq(Oa2,Ob2,delta_1_OO,qO,qO);
    
    const double OC_QQ =
            qq(Oa1,Cb,delta_1_OC,qO,qC)
          + qq(Oa2,Cb,delta_1_OC,qO,qC)
          + qq(Ob1,Ca,delta_1_OC,qO,qC)
          + qq(Ob2,Ca,delta_1_OC,qO,qC);
    
    const double CC_QQ =
            qq(Ca,Cb,delta_1_CC,qC,qC);
    
    const double OM_QQ =
            qq(Oa1,Bb1,delta_1_OM,qO,qM)
          + qq(Oa1,Bb2,delta_1_OM,qO,qM)
          + qq(Oa2,Bb1,delta_1_OM,qO,qM)
          + qq(Oa2,Bb2,delta_1_OM,qO,qM)
          + qq(Ob1,Ba1,delta_1_OM,qO,qM)
          + qq(Ob1,Ba2,delta_1_OM,qO,qM)
          + qq(Ob2,Ba1,delta_1_OM,qO,qM)
          + qq(Ob2,Ba2,delta_1_OM,qO,qM);

    const double CM_QQ =
            qq(Ca,Bb1,delta_1_CM,qC,qM)
          + qq(Ca,Bb2,delta_1_CM,qC,qM)
          + qq(Cb,Ba1,delta_1_CM,qC,qM)
          + qq(Cb,Ba2,delta_1_CM,qC,qM);

    const double MM_QQ =
            qq(Ba1,Bb1,delta_1_MM,qM,qM)
          + qq(Ba1,Bb2,delta_1_MM,qM,qM)
          + qq(Ba2,Bb1,delta_1_MM,qM,qM)
          + qq(Ba2,Bb2,delta_1_MM,qM,qM);
    
    // C6
    const double OO_C6 =
            xC6(Oa1,Ob1,delta_6_OO,C6_OO)
          + xC6(Oa1,Ob2,delta_6_OO,C6_OO)
          + xC6(Oa2,Ob1,delta_6_OO,C6_OO)
          + xC6(Oa2,Ob2,delta_6_OO,C6_OO);
    
    const double OC_C6 =
            xC6(Oa1,Cb,delta_6_OC,C6_OC)
          + xC6(Oa2,Cb,delta_6_OC,C6_OC)
          + xC6(Ob1,Ca,delta_6_OC,C6_OC)
          + xC6(Ob2,Ca,delta_6_OC,C6_OC);
    
    const double CC_C6 =
            xC6(Ca,Cb,delta_6_CC,C6_CC);
    
    // C8
    const double OO_C8 =
            xC8(Oa1,Ob1,delta_8_OO,C8_OO)
          + xC8(Oa1,Ob2,delta_8_OO,C8_OO)
          + xC8(Oa2,Ob1,delta_8_OO,C8_OO)
          + xC8(Oa2,Ob2,delta_8_OO,C8_OO);
    
    const double OC_C8 =
            xC8(Oa1,Cb,delta_8_OC,C8_OC)
          + xC8(Oa2,Cb,delta_8_OC,C8_OC)
          + xC8(Ob1,Ca,delta_8_OC,C8_OC)
          + xC8(Ob2,Ca,delta_8_OC,C8_OC);
    
    const double CC_C8 =
            xC8(Ca,Cb,delta_8_CC,C8_CC);

    return (OO_exp + OC_exp + CC_exp
          + OO_QQ  + OC_QQ  + CC_QQ
          + OM_QQ  + CM_QQ  + MM_QQ
          - OO_C6  - OC_C6  - CC_C6
          - OO_C8  - OC_C8  - CC_C8);

}

} // namespace x2b

////////////////////////////////////////////////////////////////////////////////
