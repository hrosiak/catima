#include "catima/reactions.h"

#ifdef NUREX
#include "nurex/Parametrization.h"
#include "catima/catima.h"
#include "catima/abundance_database.h"
#include "catima/storage.h"
#include <cmath>
#include <iostream>
namespace catima{
    
double nonreaction_rate(Projectile &projectile, const Material &target, const Config &c){

    if(projectile.T<emin_reaction)return -1.0;
    if(target.thickness()<=0.0)return 1.0;

    int ap = lround(projectile.A);
    int zp = lround(projectile.Z);
    int zt = target.get_element(0).Z;
    int at = abundance::get_isotope_a(zt,0); // most abundand natural isotope mass

    auto data = _storage.Get(projectile,target,c);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    if(energy_out(projectile.T, target.thickness(), range_spline) < emin_reaction)return -1.0;
    auto sigma_r = [&](double th){
        double stn_sum=0.0, sum=0.0;
        double e = energy_out(projectile.T, th, range_spline);
        for(unsigned int i = 0;i<target.ncomponents();i++){
            stn_sum += target.molar_fraction(i);
            sum += target.molar_fraction(i)*nurex::SigmaR_Kox(ap,zp,e,at,zt); 
        }
        return sum/stn_sum;
    };

    //nurex::Nucleus nurex_projectile = nurex::get_default_nucleus(ap,zp);
    //nurex::Nucleus nurex_target = nurex::get_default_nucleus(at,zt);
    //nurex::GlauberModelOLA_ZeroRange gm(nurex_projectile, nurex_target);
    //double cs = nurex::SigmaR_Kox(ap,zp,projectile.T,);

    double cs0 = sigma_r(0);
    double cs1 = sigma_r(target.thickness());
    double cs;
    if(std::abs(cs0-cs1)/cs0 < 0.05){
        cs = target.number_density_cm2()*(cs0 + cs1)/2.0;
    }
    else{
        cs = catima::integrator.integrate(sigma_r,0,target.number_density_cm2());
    }
    
    return exp(-cs*0.0001);
    }

}

#endif
