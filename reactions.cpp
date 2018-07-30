#include "catima/reactions.h"

#ifdef NUREX
#include "nurex/Parametrization.h"
#include "catima/catima.h"
#include "catima/abundance_database.h"
#include "catima/storage.h"
#include <cmath>
#include <iostream>
namespace catima{
    
double reaction_rate1(Projectile &projectile, const Material &target, const Config &c){
    int ap = lround(projectile.A);
    int zp = lround(projectile.Z);
    int zt = target.get_element(0).Z;
    int at = abundance::get_isotope_a(zt,0); // most abundand natural isotope mass

    auto data = _storage.Get(projectile,target,c);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);

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
    double cs = catima::integrator.integrate(sigma_r,0,target.thickness());

    double rr = reaction_rate(cs,target.number_density_cm2(0));
    std::cout<<rr<<"\n";
    
    return 1.0;
    }

}

#endif
