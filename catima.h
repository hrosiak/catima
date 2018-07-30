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

#ifndef CPPATIMA_H
#define CPPATIMA_H

#include <utility>
#include <vector>

// #define NDEBUG
#include "catima/config.h"
#include "catima/constants.h"
#include "catima/structures.h"
#include "catima/calculations.h"
#include "catima/material_database.h"
#include "catima/storage.h"

namespace catima{
    
    /**
      * calculate dEdx for projectile-Material combination
      * @param p - Projectile
      * @param mat - Material
      * @return dEdx
      */
    double dedx(Projectile &p, double T, const Material &mat, const Config &c=default_config);

    /**
      * calculate energy loss straggling variance for projectile-Material combination
      * @param p - Projectile
      * @param mat - Material
      * @return dOmega^2/dx
      */
    double domega2dx(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * calculates variance of angular scattering of Projectile p on Material m
      */
    double da2dx(Projectile &p, double T, const Material &m, const Config &c=default_config);

    /**
      * returns the range of the Projectile in Material calculated from range spline
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return range
      */
    double range(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * returns the dEdx calculated from range spline as derivative
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return range
      */
    double dedx_from_range(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * returns the dEdx calculated from range spline as derivative
      * @param p - Projectile
      * @param T - energy vector
      * @param mat - Material
      * @return range
      */
    std::vector<double> dedx_from_range(Projectile &p, const std::vector<double> &T, const Material &t, const Config &c=default_config);

    /**
      * returns the  range straggling of the Projectile in Material from spline
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return range straggling
      */
    double range_straggling(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * returns the  range variance of the Projectile in Material from spline
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return range straggling
      */
    double range_variance(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * returns the  range variance per dE, calculated as derivative of range variance spline
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return range variance / dE
      */
    double domega2de(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * returns the  angular variance per dE, calculated as derivative of angular variance spline
      * @param p - Projectile
      * @param T - energy in MeV/u
      * @param mat - Material
      * @return angular variance / dE
      */
    double da2de(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * calculates angular scattering in the material from difference of incoming a nd outgoing energies
      * @param p - Projectile
      * @param T - incoming energy
      * @param Tout - outcoming energy
      * @param mat - Material
      * @return angular straggling
      */
    double angular_straggling_from_E(Projectile &p, double T, double Tout,const Material &t, const Config &c=default_config);

    /**
      * calculates Energy straggling in the material from difference of incoming a nd outgoing energies
      * @param p - Projectile
      * @param T - incoming energy
      * @param Tout - outcoming energy
      * @param mat - Material
      * @return angular straggling
      */
    double energy_straggling_from_E(Projectile &p, double T, double Tout,const Material &t, const Config &c=default_config);

    /**
      * calculates outcoming energy from range spline
      * @param T - incoming energy
      * @thickness - thicnkess of the target in g/cm2
      * @range_spline - precaclulated range spline for material 
      * @return outcoming energy after the thickness in Mev/u
      */
    double energy_out(double T, double thickness, Interpolator &range_spline);

    /**
      * calculates outcoming energy 
      * @p - Projectile
      * @t - Material
      * @param T - incoming energy
      * @return outcoming energy after the material in Mev/u
      */
    double energy_out(Projectile &p, double T, const Material &t, const Config &c=default_config);

    /**
      * calculates outcoming energy 
      * @p - Projectile
      * @t - Material
      * @param T - incoming energy vector
      * @return outcoming energy after the material in Mev/u
      */
    std::vector<double> energy_out(Projectile &p, const std::vector<double> &T, const Material &t, const Config &c=default_config);

    /**
      * calculates all observables for projectile passing material
      * @param p - Projectile
      * @param mat - Material
      * @return structure of Result
      */
    Result calculate(Projectile &p, const Material &t, const Config &c=default_config);
    inline Result calculate(Projectile &p, const Material &t, double T, const Config &c=default_config){
        p.T = T;
        return calculate(p, t, c);
    }

      /**
      * wrapper to other calculate function with simplified arguments
      * @param p - Projectile
      * @param mat - Material
      * @return structure of Result
      */
    Result calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density);


    /**
      * calculate observables for multiple layers of material defined by Layers
      * @return results stored in MultiResult structure
      *
      */
    MultiResult calculate(Projectile &p, const Layers &layers, const Config &c=default_config);
    inline MultiResult calculate(Projectile &p, double T, const Layers &layers, const Config &c=default_config){
        p.T = T;
        return calculate(p, layers, c);
    }

    /// the following functions are used to calculates array of data points for whole range of energies
    /// usually used to construct splines
    std::vector<double> calculate_range(Projectile p, const Material &t, const Config &c=default_config);
    std::vector<double> calculate_range_straggling(Projectile p, const Material &t, const Config &c=default_config);
    std::vector<double> calculate_angular_variance(Projectile p, const Material &t, const Config &c=default_config);
    std::vector<double> calculate_tof(Projectile p, const Material &t, const Config &c=default_config);
    
    /**
      * calculates TOF of the Projectile in Material
      * this is used instead of precalculated TOF spline
      * @return TOF in ns
      */
    double calculate_tof_from_E(Projectile p, double Eout, const Material &t, const Config &c=default_config);
    
    /**
     * returns energy magnification after passing material t
     */
    std::pair<double,double> w_magnification(Projectile p, double Ein, const Material &t, const Config &c=default_config);

    class DataPoint;
    /**
      * calculates DataPoint for Projectile Material combinatino
      * it substitute series of calls to calculate_* functions
      * they are all combined here in 1 single function
      * it has a perfomance gain to call this function if all splines are to be caclulated
      */
    DataPoint calculate_DataPoint(Projectile p, const Material &t, const Config &c=default_config);

    bool operator==(const Config &a, const Config&b);
}
#endif
