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

#ifndef REACTIONS_H
#define REACTIONS_H
#include "catima/build_config.h"
#ifdef NUREX
#include "nurex/nurex.h"
#include "catima/structures.h"
#include "catima/config.h"
#include "catima/integrator.h"
#include <cmath>

namespace catima{
    
    /**
     * return reaction probability 
     * @param sigma - cross section in mb
     * @param t - number of targets per cm2 in 10^23 unit 
     */
    inline double reaction_rate(double sigma, double t){
        return 1.0 - std::exp(-sigma*t*0.0001);
    }

    /**
     * return nonreaction rate
     * @param sigma - cross section in mb
     * @param t - number of targets per cm2 in 10^23 unit 
     */
    inline double nonreaction_rate(double sigma, double t){
        return std::exp(-sigma*t*0.0001);
    }
    
    template<typename F>
    double reaction_rate(F& f, double t){
        GaussLegendreIntegration<8> ii;        
        double i = ii.integrate(f,0,t);
        return 1.0 - std::exp(-i*0.0001);
    }

    double reaction_rate1(Projectile &projectile, const Material &target, const Config &c=default_config);
    
}

#endif //NUREX
#endif
