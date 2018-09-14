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
#include "catima/structures.h"
#include "catima/config.h"
#include "catima/integrator.h"
#include <cmath>
#endif

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
    double nonreaction_rate(Projectile &projectile, const Material &target, const Config &c=default_config);
    double production_rate(double cs, double rcs_projectile, double rcs_product, const Material &target, const Config &c=default_config);
    
#ifndef NUREX
double SigmaR_Kox(int Ap, int Zp, double E, int At, int Zt);
inline double p_from_T(double T, double M=1.0){
    return M*sqrt(T*T + 2*T*atomic_mass_unit);
}
inline double Ecm_from_T_relativistic(double T, double Ap, double At){
    double mp = Ap*atomic_mass_unit;
    double mt = At*atomic_mass_unit;
    double plab= p_from_T(T,Ap);
    double elab = sqrt(plab*plab + mp*mp);
    double ecm = sqrt(mp*mp + mt*mt + 2*elab*mt);
    double pcm = plab * mt / ecm;
    return sqrt(pcm*pcm+mp*mp)-mp;
}

#endif //NUREX
#endif
} // end of catime namespace
