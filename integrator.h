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
#include <functional>
#include <array>

//#ifdef USE_THREADS
//#include <mutex>
//#endif
#include "catima/glq_integrator.h"
#include "catima/gkq_integrator.h"
#ifdef GSL_INTEGRATION
#include "gsl/gsl_integration.h"
#endif


using integrators::GaussLegendreIntegration;
using integrators::GaussKronrodIntegration;

namespace catima{
    
#ifdef GSL_INTEGRATION
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
#endif

#ifdef GSL_INTEGRATION
using integrator_type = IntegratorGSL;
#else
using integrator_type = GaussLegendreIntegration<8>;
#endif
using integrator_adaptive_type = GaussKronrodIntegration<21>;

extern integrator_type integrator;
extern integrator_adaptive_type integrator_adaptive;
}

#endif
