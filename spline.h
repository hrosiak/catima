/*
 * 
 * This is modification of Tino Kluge tk spline
 * calculation is optimized for tridiagonal matrices
 * 
 *  Copyright(C) 2017
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CATIMA_SPLINE_H
#define CATIMA_SPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <array>
#include "catima/constants.h"

#ifdef GSL_INTERPOLATION
#include <gsl/gsl_spline.h>
#endif

namespace catima
{

enum interpolation_t {cspline, linear};

/**
 * Tridiagonal matrix solver
 */

template<int N>
class tridiagonal_matrix
{
private:
    std::array<double,N> a;
    std::array<double,N> d;
    std::array<double,N> c;
public:
    tridiagonal_matrix() {}

    // access operator
    double & operator () (unsigned int i, unsigned int j);            // write
    double   operator () (unsigned int i, unsigned int j) const;      // read
    std::array<double, N> trig_solve(const std::array<double, N>& b) const;
};

template<int N>
double & tridiagonal_matrix<N>::operator () (unsigned int i, unsigned int j)
{
    int k=j-i;
    if(k == -1)return c[i];
    else if(k==0) return d[i];
    else return a[i];
}

template<int N>
double tridiagonal_matrix<N>::operator () (unsigned int i, unsigned int j) const
{
    int k=j-i;
    if(k==-1)return c[i];
    else if(k==0) return d[i];
    else if(k==1)return a[i];
    else return 0.0;
}

template<int N>
std::array<double, N> tridiagonal_matrix<N>::trig_solve(const std::array<double, N>& b) const
{
    std::array<double, N> x;
    if(d[0] == 0.0){return x;}
    std::array<double, N> g;
    x[0] = b[0]/d[0];
    double bet = d[0];
    for(std::size_t j=1, max=N;j<max;j++){
        g[j] = c[j-1]/bet;
        bet = d[j] - (a[j]*g[j]);
        if(bet == 0.0){
            x.fill(0.0);
            return x;
        }
        x[j] = (b[j]-a[j]*x[j-1])/bet;
    }
    for(int j=N-2;j>=0;j--){
        x[j] -= g[j+1]*x[j+1];
    }
    return x;
}


/**
 * Cubic Spline class, accepting EnergyTable type as x-variable
 */
template<typename T>
struct cspline_special{
    constexpr static int N = T::size();
    cspline_special(const T& x,
                    const std::vector<double>& y,
                    bool boundary_second_deriv = true);
    cspline_special() = default;

    const T *table;
    const double *m_x;
    const double *m_y;
    std::array<double,N> m_a,m_b,m_c;
    double  m_b0, m_c0;


    double operator()(double x)const{return evaluate(x);}
    double evaluate(double x) const
    {
        int idx=std::max( table->index(x), 0);
        double h=x-m_x[idx];
        double interpol;
        if(x<m_x[0]) {
            // extrapolation to the left
            interpol=(m_b0*h + m_c0)*h + m_y[0];
        } else if(x>m_x[N-1]) {
            // extrapolation to the right
            interpol=(m_b[N-1]*h + m_c[N-1])*h + m_y[N-1];
        } else {
            // interpolation
            interpol=((m_a[idx]*h + m_b[idx])*h + m_c[idx])*h + m_y[idx];
        }
        return interpol;
    }


    double deriv(double x) const
    {
        int idx=std::max( table->index(x), 0);

        double h=x-m_x[idx];
        double interpol;
        if(x<m_x[0]) {
            // extrapolation to the left
            interpol=2.0*m_b0*h + m_c0;
        } else if(x>m_x[N-1]) {
            // extrapolation to the right
            interpol=2.0*m_b[N-1]*h + m_c[N-1];
        } else {
            // interpolation
            interpol=(3.0*m_a[idx]*h + 2.0*m_b[idx])*h + m_c[idx];
        }
        return interpol;
    }
    static_assert (T::size()>2, "N must be > 2");
};

template<typename T>
cspline_special<T>::cspline_special(const T &x,
                      const std::vector<double>& y,
                      bool boundary_second_deriv
                      ):table(&x),m_y(y.data()),m_x(x.values)
{
    static_assert (N>2, "N must be > 2");
    tridiagonal_matrix<N> A{};
    std::array<double, N> rhs;
    for(std::size_t i=1; i<N-1; i++) {
        A(i,i-1)=1.0/3.0*(x[i]-x[i-1]);
        A(i,i)=2.0/3.0*(x[i+1]-x[i-1]);
        A(i,i+1)=1.0/3.0*(x[i+1]-x[i]);
        rhs[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    }
        // boundary conditions
    if(boundary_second_deriv) {
        // 2*b[0] = f''
        A(0,0)=2.0;
        A(0,1)=0.0;
        rhs[0]=0.0; // 0.0 is value of derivative
        A(N-1,N-1)=2.0;
        A(N-1,N-2)=0.0;
        rhs[N-1]=0.0; // 0.0 is value of derivative
       } else {
        // c[0] = f', needs to be re-expressed in terms of b:
        // (2b[0]+b[1])(x[1]-x[0]) = 3 ((y[1]-y[0])/(x[1]-x[0]) - f')
        A(0,0)=2.0*(x[1]-x[0]);
        A(0,1)=1.0*(x[1]-x[0]);
        rhs[0]=3.0*((y[1]-y[0])/(x[1]-x[0])-0.0); // 0.0 is deriv value

        // c[n-1] = f', needs to be re-expressed in terms of b:
        // (b[n-2]+2b[n-1])(x[n-1]-x[n-2])
        // = 3 (f' - (y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
        A(N-1,N-1)=2.0*(x[N-1]-x[N-2]);
        A(N-1,N-2)=1.0*(x[N-1]-x[N-2]);
        rhs[N-1]=3.0*(0.0-(y[N-1]-y[N-2])/(x[N-1]-x[N-2]));
    }


        m_b=A.trig_solve(rhs);

        // calculate parameters a[] and c[] based on b[]
        for(int i=0; i<N-1; i++) {
            m_a[i]=1.0/3.0*(m_b[i+1]-m_b[i])/(x[i+1]-x[i]);
            m_c[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
                   - 1.0/3.0*(2.0*m_b[i]+m_b[i+1])*(x[i+1]-x[i]);
        }


    // for left extrapolation coefficients
    //s.m_b0 = (m_force_linear_extrapolation==false) ? s.m_b[0] : 0.0;
    m_b0 =  0.0;
    m_c0 = m_c[0];

    double h=x[N-1]-x[N-2];
    m_a[N-1]=0.0;
    m_c[N-1]=3.0*m_a[N-2]*h*h+2.0*m_b[N-2]*h+m_c[N-2];   // = f'_{n-2}(x_{n-1})
    m_b[N-1]=0.0;
}

} // namespace end

#endif
