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

template<int order>
struct GL_data{
};

template<int order>
class GaussLegendreIntegration{
public:
    template<typename F>
    double integrate(F f, double a, double b) const;
    template<typename F>
    double operator()(F f, double a, double b) const {return integrate(f, a, b);}
    double w(int i) const {return GL_data<order>::w()[i];}
    double x(int i) const {return GL_data<order>::x()[i];}
    int n() const {return order;}
    std::array<double,order> get_points(double a = -1.0, double b = 1.0)const;
};

template<int order>
template<typename F>
double GaussLegendreIntegration<order>::integrate(F f, double a, double b) const{
    double res=0.0;
    double p = 0.5*(b-a);
    double q = 0.5*(b+a);

    if(order%2){res+= w(0) * (f(p*x(0) + q));} // in case odd-order
    for(int i=order%2;i<order/2 + order%2;i++){
        res += w(i) * (f(p*x(i) + q) + f(-p*x(i) + q));
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
        points[num-i-1] = -p*x(i) + q;
        points[num+i] = p*x(i) + q;
    }
    return points;
}

// order = 8
template<>
struct GL_data<8>{
    static const std::array<double,4>& x(){
        static const std::array<double,4> _x = 
            {0.1834346424956498049394761,0.5255324099163289858177390,0.7966664774136267395915539,0.9602898564975362316835609};
        return _x;
    }
    static const std::array<double,4>& w(){
        static const std::array<double,4> _w = 
            {0.3626837833783619829651504,0.3137066458778872873379622,0.2223810344533744705443560,0.1012285362903762591525314};
        return _w;
    }
};

// order = 10
template<>
struct GL_data<10>{
    static const std::array<double,5>& x(){
        static const std::array<double,5> _x = 
            {0.1488743389816312108848260,0.4333953941292471907992659,0.6794095682990244062343274,0.8650633666889845107320967,0.9739065285171717200779640};
        return _x;
    }
    static const std::array<double,5>& w(){
        static const std::array<double,5> _w = 
            {0.2955242247147528701738930,0.2692667193099963550912269,0.2190863625159820439955349,0.1494513491505805931457763,0.0666713443086881375935688};
        return _w;
    }
};

// order  = 12
template<>
struct GL_data<12>{
    static const std::array<double,6>& x(){
        static const std::array<double,6> _x = 
            {0.1252334085114689154724414,0.3678314989981801937526915,0.5873179542866174472967024,0.7699026741943046870368938,0.9041172563704748566784659,0.9815606342467192506905491};
        return _x;
    }
    
    static const std::array<double,6>& w(){
        static const std::array<double,6> _w = 
            {0.2491470458134027850005624,0.2334925365383548087608499,0.2031674267230659217490645,0.1600783285433462263346525,0.1069393259953184309602547,0.0471753363865118271946160};
        return _w;
    }
};

// order  = 14
template<>
struct GL_data<14>{
    static const std::array<double,7>& x(){
        static const std::array<double,7> _x = 
            {0.1080549487073436620662447,0.3191123689278897604356718,0.5152486363581540919652907,0.6872929048116854701480198,0.8272013150697649931897947,0.9284348836635735173363911,0.9862838086968123388415973};
        return _x;
    }
    
    static const std::array<double,7>& w(){
        static const std::array<double,7> _w = 
            {0.2152638534631577901958764,0.2051984637212956039659241,0.1855383974779378137417166,0.1572031671581935345696019,0.1215185706879031846894148,0.0801580871597602098056333,0.0351194603317518630318329};
        return _w;
    }
};

// order = 16
template<>
struct GL_data<16>{
    static const std::array<double,8>& x(){
        static const std::array<double,8> _x = 
           {0.0950125098376374401853193,0.2816035507792589132304605,0.4580167776572273863424194,0.6178762444026437484466718,0.7554044083550030338951012,0.8656312023878317438804679,0.9445750230732325760779884,0.9894009349916499325961542};
        return _x;
    }
    
    static const std::array<double,8>& w(){
        static const std::array<double,8> _w = 
           {0.1894506104550684962853967,0.1826034150449235888667637,0.1691565193950025381893121,0.1495959888165767320815017,0.1246289712555338720524763,0.0951585116824927848099251,0.0622535239386478928628438,0.0271524594117540948517806};
        return _w;
    }
};

#ifdef GSL_INTEGRATION
using integrator_type = IntegratorGSL;
#else
using integrator_type = GaussLegendreIntegration<8>;
#endif

extern integrator_type integrator;
extern IntegratorGSL integratorGSL;
}

#endif
