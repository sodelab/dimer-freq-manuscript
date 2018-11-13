#ifndef X1B_H
#define X1B_H

#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include "kit.h" 
#include "poly-1b.h"
#include "x1b.h"

////////////////////////////////////////////////////////////////////////////////
namespace x1b {

void cart_to_vars(const double* xyz, double* v) 
{
    const double* C  = xyz;
    const double* O1 = xyz + 3;
    const double* O2 = xyz + 6;

    using kit::distance;

    v[0]  = var_intra(d0_intra_CO, m_k_CO_intra, distance( C, O1));
    v[1]  = var_intra(d0_intra_CO, m_k_CO_intra, distance( C, O2));
    v[2]  = var_intra(d0_intra_OO, m_k_OO_intra, distance(O1, O2));
} 

//----------------------------------------------------------------------------//

double value(const double xyz[9])
{
    double v[num_vars];
    cart_to_vars(xyz, v);

    double E_poly(0);

    {
        double mono[num_linear_params];
        eval(v, mono);

        for (size_t n = 0; n < num_linear_params; ++n)
            E_poly += poly[n]*mono[n];
    }

    return E_poly;
}

} // namespace x1b

////////////////////////////////////////////////////////////////////////////////

#endif // X1B_H
