#ifndef X2B_H
#define X2B_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include "kit.h" 
#include "poly-2b.h"
#include "x2b.h"

////////////////////////////////////////////////////////////////////////////////
namespace x2b {

double f_switch(const double* xyz)
{
    //find smallest intermolecular distance
    double r = kit::distance_short(xyz);

    if (r > m_r2f_long) {
        return 0.0;
    } else if (r > m_r2i_long) {
        const double x = (r - m_r2i_long)/(m_r2f_long - m_r2i_long);
        return (1.0 + std::cos(M_PI*x))/2;
    } else if (r > m_r2f_small) {
        return 1.0;
    } else if (r > m_r2i_small) {
        const double x = (r - m_r2i_small)/(m_r2f_small - m_r2i_small);
        return std::pow(x,1/5);
    } else {
        return 0.0;
    }
}

void cart_to_vars(const double* xyz, double* v, double& s)
{
    const double* Ca  = xyz;
    const double* Oa1 = xyz + 3;
    const double* Oa2 = xyz + 6;

    const double* Cb  = xyz + 9;
    const double* Ob1 = xyz + 12;
    const double* Ob2 = xyz + 15;

    using kit::distance;

    v[0]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Ca, Oa1));
    v[1]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Ca, Oa2));
    v[2]  = var_intra(d0_intra_OO, m_k_OO_intra, distance( Oa1, Oa2));

    v[3]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Cb, Ob1));
    v[4]  = var_intra(d0_intra_CO, m_k_CO_intra, distance(  Cb, Ob2));
    v[5]  = var_intra(d0_intra_OO, m_k_OO_intra, distance( Ob1, Ob2));

    v[6]  = var_inter(d0_inter_CC, m_k_CC_inter, distance( Ca,  Cb));

    v[7]  = var_inter(d0_inter_CO, m_k_CO_inter, distance( Ca, Ob1));
    v[8]  = var_inter(d0_inter_CO, m_k_CO_inter, distance( Ca, Ob2));
    v[9]  = var_inter(d0_inter_CO, m_k_CO_inter, distance( Cb, Oa1));
    v[10] = var_inter(d0_inter_CO, m_k_CO_inter, distance( Cb, Oa2));

    v[11] = var_inter(d0_inter_OO, m_k_OO_inter, distance( Oa1, Ob1));
    v[12] = var_inter(d0_inter_OO, m_k_OO_inter, distance( Oa1, Ob2));
    v[13] = var_inter(d0_inter_OO, m_k_OO_inter, distance( Oa2, Ob1));
    v[14] = var_inter(d0_inter_OO, m_k_OO_inter, distance( Oa2, Ob2));

    s = f_switch(xyz);

} 

double value(const double xyz[18])
{
    double v[num_vars], s;
    cart_to_vars(xyz, v, s);

    double E_poly(0);

    {   
        double mono[num_linear_params];
        eval(v, mono);
		
        for (size_t n = 0; n < num_linear_params; ++n) 
            E_poly += poly[n]*mono[n];
    }

    return s*E_poly + sapt_s(xyz);
}

}
#endif // X2B_H

////////////////////////////////////////////////////////////////////////////////
