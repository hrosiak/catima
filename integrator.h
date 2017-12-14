/*
 *  Author: Andrej Prochazka
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

#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "catima/build_config.h"
#include "gsl/gsl_integration.h"
#include <functional>
#include <array>
#ifdef USE_THREADS
#include <mutex>
#endif

namespace catima{
    
/// helper class to integrate functions using the GSL library
class IntegratorGSL{
	public:
	IntegratorGSL(bool adapt=true);
	~IntegratorGSL();

	double integrate(std::function<double(double)> f, double min, double  max, double precision=0.001);

	private:
    gsl_integration_workspace *w;
    bool adaptive;
	double error;
	double result;
	double min;
	double max;
	#ifdef USE_THREADS
	std::mutex integration_mutex;
	#endif
    };

// built in integrator 
template<int order>
class GaussLegendreIntegration{
public:
    template<typename F>
    double integrate(F f, double a, double b) const;
    template<typename F>
    double operator()(F f, double a, double b) const {return integrate(f, a, b);}
    double get_w(int i) const {return w[i];}
    double get_x(int i) const {return x[i];}
    int n() const {return order;}
    std::array<double,order> get_points(double a = -1.0, double b = 1.0)const;

public:
    static std::array<double,order/2> w;
    static std::array<double,order/2> x;
};

template<int order>
template<typename F>
double GaussLegendreIntegration<order>::integrate(F f, double a, double b) const{
    double res=0.0;
    double p = 0.5*(b-a);
    double q = 0.5*(b+a);

    for(int i=0;i<order/2;i++){
        res += w[i] * (f(p*x[i] + q) + f(-p*x[i] + q));
    }
    return p*res;
}

template<int order>
std::array<double,order> GaussLegendreIntegration<order>::get_points(double a,  double b)const{
    std::array<double,order> points;
    double p = 0.5*(b-a);
    double q = 0.5*(b+a);
    
    int num = (order/2);
    for(int i=0;i< num;i++){
        points[num-i-1] = -p*x[i] + q;
        points[num+i] = p*x[i] + q;
    }
    return points;
}

#ifdef GSL_INTEGRATION
using integrator_type = IntegratorGSL;
#else
using integrator_type = GaussLegendreIntegration<8>;
#endif

extern integrator_type integrator;
extern IntegratorGSL integratorGSL;
}

#endif
