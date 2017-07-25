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
#include "gsl/gsl_integration.h"
#include <functional>
#ifdef USE_THREADS
#include <mutex>
#endif

namespace catima{
    
/// helper class to integrate functions using the GSL library
class IntegratorGSL{
	public:
	IntegratorGSL();
	~IntegratorGSL();

	double Integrate(std::function<double(double)> f, double min, double  max, double precision, bool adaptive=true);

	private:
	gsl_integration_workspace *w;
	double error;
	double precision;
	double result;
	double min;
	double max;
	#ifdef USE_THREADS
	std::mutex integration_mutex;
	#endif
    };

extern IntegratorGSL integratorGSL;
    
}

#endif
