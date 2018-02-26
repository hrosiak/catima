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

#ifndef TK_SPLINE_H
#define TK_SPLINE_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>

namespace catima
{

// band matrix solver
class band_matrix
{
private:
    std::vector<double> a;
    std::vector<double> d;
    std::vector<double> c;
    std::vector<double> save;
public:
    band_matrix() {};                             // constructor
    band_matrix(int dim);       // constructor
    ~band_matrix() {};                            // destructor
    void resize(int dim);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension

    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diogonal (in m_lower)
    std::vector<double> trig_solve(const std::vector<double>& b) const;
};


// spline interpolation
class spline
{
public:
    enum class bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    const double *m_x, *m_y;            // x,y coordinates of points
    size_t n=0;
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    std::vector<double> m_a,m_b,m_c;        // spline coefficients
    double  m_b0, m_c0;                     // for left extrapol
    bd_type m_left = bd_type::second_deriv;
    bd_type m_right = bd_type::second_deriv;
    double  m_left_value = 0.0;
    double  m_right_value = 0.0;
    bool    m_force_linear_extrapolation = false;

public:
    // set default boundary condition to be zero curvature at both ends
    spline(){}

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const double *x,
                    const double *y,
                    const size_t num);
    double operator() (double x) const;
    double deriv(int order, double x) const;
};

} // namespace tk

#endif /* TK_SPLINE_H */
