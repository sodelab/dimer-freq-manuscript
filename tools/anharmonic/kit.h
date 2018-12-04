#ifndef KIT_H
#define KIT_H


#include <cmath>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <iomanip>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////

namespace kit {

inline double distance(const double* p1, const double* p2)
{
    double result(0);

    for (size_t k = 0; k < 3; ++k) {
        const double delta = p1[k] - p2[k];
        result += delta*delta;
    }

    return std::sqrt(result);
}

//----------------------------------------------------------------------------//

inline double distance_short(const double *xyz)
{
    const double* Ca  = xyz;
    const double* Oa1 = xyz + 3;
    const double* Oa2 = xyz + 6;

    const double* Cb  = xyz + 9;
    const double* Ob1 = xyz + 12;
    const double* Ob2 = xyz + 15;

    double r_small = distance( Ca, Cb);
    if (distance( Ca,Ob1) < r_small) r_small = distance( Ca,Ob1);
    if (distance( Ca,Ob2) < r_small) r_small = distance( Ca,Ob2);
    if (distance(Oa1,Ob1) < r_small) r_small = distance(Oa1,Ob1);
    if (distance(Oa1,Ob2) < r_small) r_small = distance(Oa1,Ob2);
    if (distance(Oa1, Cb) < r_small) r_small = distance(Oa1, Cb);
    if (distance(Oa2,Ob1) < r_small) r_small = distance(Oa2,Ob1);
    if (distance(Oa2,Ob2) < r_small) r_small = distance(Oa2,Ob2);
    if (distance(Oa2, Cb) < r_small) r_small = distance(Oa2, Cb);

    return r_small;

} 

//----------------------------------------------------------------------------//

template <int N>
inline int factorial()
{
    return N*factorial<N-1>();
}

template<>
inline int factorial<0>()
{
    return 1;
}

//----------------------------------------------------------------------------//

inline double tang_toennies(int n, const double& x)
{
//
//  T(n, x) = 1 - exp(-x)*(1 + x/1! + x^2/2! + x^3/3! + ... + x^n/n!)
//  diff(T(n, x), x) = T(n - 1, x) - T(n, x) = exp(-x)*x^n/n!
//
    assert(n >= 0);

    int nn = n;

    double sum = 1.0 + x/nn;
    while (--nn != 0)
        sum = 1.0 + sum*x/nn;

    double tt = 1.0 - sum*std::exp(-x);

    if (std::fabs(tt) < 1.0e-8) {

        double term(1);
        for (nn = n; nn != 0; --nn)
            term *= x/nn;

        sum = 0.0;
        for (nn = n + 1; nn < 1000; ++nn) {
            term *= x/nn;
            sum += term;

            if (std::fabs(term/sum) < 1.0e-8)
                break;
        }

        tt = sum*std::exp(-x);
    }

    return tt;
}
    

void inline to_cartesian(const double *v, double *pos, size_t ncoord)
{

    size_t nfrags = ncoord/6;
    size_t natoms = nfrags*3;
    //  double crd[ncoord];
    //  for (size_t n = 0; n < ncoord; ++n) {
    //    crd[n] = gsl_vector_get(v->x, n);
    //  }
    //    for (size_t n = 0; n < nfrags; ++n)
    //        std::cout
    //        << std::setw(20) << crd[6*n + 0]
    //        << std::setw(20) << crd[6*n + 1]
    //        << std::setw(20) << crd[6*n + 2]
    //        << std::setw(20) << crd[6*n + 3]
    //        << std::setw(20) << crd[6*n + 4]
    //        << std::setw(20) << crd[6*n + 5]
    //        << '\n';
    
    // Initialize fragment position
    // carbon atom won't move
    //  std::cout<<natoms<<std::endl;
    for (size_t i=0; i<nfrags; ++i) {
        //    std::cout << natoms/nfrags*3<<std::endl;
        for (size_t n=0; n<natoms/nfrags*3; ++n) {
            // std::cout <<6*i+n%3<< " " <<  v[6*i+n%3]<< std::endl;
            pos[9*i+n] = v[6*i+n%3]; // carbon atom positions
        }
    }
    
    // Determine position of oxygen one
    // Use the values: p,phi,theta
    // Convert spherical coordinates to cartesian
    for (size_t i=0; i<nfrags; ++i) {
        //    double r     =  v[6*i+3];
        double r     =  1.162047;
        double phi   =  v[6*i+4];
        double theta =  v[6*i+5];
        pos[9*i+3+0] += r*std::sin(theta)*std::cos(phi);
        pos[9*i+3+1] += r*std::sin(theta)*std::sin(phi);
        pos[9*i+3+2] += r*std::cos(theta);
        
        //    pos[9*i+3+0] += r*std::sin(phi*M_PI/180)*std::cos(theta*M_PI/180);
        //    pos[9*i+3+1] += r*std::sin(phi*M_PI/180)*std::sin(theta*M_PI/180);
        //    pos[9*i+3+2] += r*std::cos(phi*M_PI/180);
    }
    
    // Determine position of oxygen two
    // Use the values: r,phi,theta
    // Convert spherical coordinates to cartesian
    for (size_t i=0; i<nfrags; ++i) {
        //    double r     = -v[6*i+3];
        double r     = -1.162047;
        double phi   =  v[6*i+4];
        double theta =  v[6*i+5];
        pos[9*i+6+0] += r*std::sin(theta)*std::cos(phi);
        pos[9*i+6+1] += r*std::sin(theta)*std::sin(phi);
        pos[9*i+6+2] += r*std::cos(theta);
        
        //    pos[9*i+6+0] += r*std::sin(phi*M_PI/180)*std::cos(theta*M_PI/180);
        //    pos[9*i+6+1] += r*std::sin(phi*M_PI/180)*std::sin(theta*M_PI/180);
        //    pos[9*i+6+2] += r*std::cos(phi*M_PI/180);
    }
    
       // for (size_t n = 0; n < natoms; ++n)
       //     std::cout << (n%3 == 0 ? "C " : "O ")
       //     << std::setw(20) << pos[3*n + 0]
       //     << std::setw(20) << pos[3*n + 1]
       //     << std::setw(20) << pos[3*n + 2]
       //     << '\n';
//     exit(0); 
    
}

void inline to_spherical(const double *pos, double *var, size_t ncoord)
{
    
    size_t nfrags = ncoord/6;
    size_t natoms = nfrags*3;
    
    
        //for (size_t n = 0; n < natoms; ++n)
        //    std::cout << (n%3 == 0 ? "C " : "O ")
        //    << std::setw(20) << pos[3*n + 0]
        //    << std::setw(20) << pos[3*n + 1]
        //    << std::setw(20) << pos[3*n + 2]
        //    << '\n';
    
    for (size_t i=0; i<nfrags; ++i) {
        for (size_t n=0; n<3; ++n) {
            //      std::cout << 9*i+n%3 << " " << pos[9*i+n%3] << std::endl;
            var[6*i+n] = pos[9*i+n%3]; // carbon atom positions
        }
        
        double x1 = pos[9*i+1*3+0]-pos[9*i+0*3+0];
        double y1 = pos[9*i+1*3+1]-pos[9*i+0*3+1];
        double z1 = pos[9*i+1*3+2]-pos[9*i+0*3+2];
        
        double r1      = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        //double phi1  = ( y1/x1 != y1/x1 ? atan2(y1,x1) : atan2(y1,x1));
        //    double theta1 = atan2(std::sqrt(x1*x1 + y1*y1)/z1);
        //double theta1 = std::acos(z1/std::sqrt(x1*x1 + y1*y1 + z1*z1));
        
        double theta1    = std::acos(z1/r1);
        double phi1      = atan2(y1,x1);
        
        //    std::cerr << theta1 << " " << phi1<< std::endl;
        //    double theta1  = ( y1/x1 != y1/x1 ? 0 : atan2(y1/x1)*180/M_PI);
        //    double phi1 = atan2(std::sqrt(x1*x1 + y1*y1)/z1)*180/M_PI;
        
        double x2 = pos[9*i+2*3+0]-pos[9*i+0*3+0];
        double y2 = pos[9*i+2*3+1]-pos[9*i+0*3+1];
        double z2 = pos[9*i+2*3+2]-pos[9*i+0*3+2];
        
        double r2     = std::sqrt(x2*x2 + y2*y2 + z2*z2);
        double phi2 = atan2(y2,x2);
        double theta2 = atan2(std::sqrt(x2*x2 + y2*y2),z2);
        
        //    double theta2 = atan2(y2/x2)*180/M_PI;
        //    double phi2 = atan2(std::sqrt(x2*x2 + y2*y2)/z2)*180/M_PI;

//	std::cout << theta1 << " " << theta2 << " THETA" << std::endl;
//        std::cout << phi1 << " " << phi2 << " PHI" << std::endl;        

        var[6*i+3] = (r1+r2)/2.0;
        var[6*i+4] = phi1;
        var[6*i+5] = theta1;
        
    }
    //    std::cout <<nfrags<<std::endl;
    //    for (size_t n = 0; n < nfrags; ++n)
    //        std::cout
    //        << std::setw(20) << var[6*n + 0]
    //        << std::setw(20) << var[6*n + 1]
    //        << std::setw(20) << var[6*n + 2]
    //        << std::setw(20) << var[6*n + 3]
    //        << std::setw(20) << var[6*n + 4]
    //        << std::setw(20) << var[6*n + 5]
    //        << " the end" << "\n";
    
}
    

} // namespace kit

////////////////////////////////////////////////////////////////////////////////

#endif // KIT_H
