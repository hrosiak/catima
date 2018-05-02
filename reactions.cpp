#include "catima/reactions.h"
#ifdef NUREX
#include "catima/catima.h"
#include "catima/abundance_database.h"
#include <cmath>
#include <iostream>
namespace catima{
    
double reaction_rate1(Projectile &projectile, const Material &target, const Config &c){

    int num_elements = target.ncomponents();
    int ap = lround(projectile.A);
    int zp = lround(projectile.Z);
    nurex::Nucleus nurex_projectile = nurex::get_default_nucleus(ap,zp);
    
    int zt = target.get_element(0).Z;
    int at = abundance::get_isotope_a(zt,0);
    nurex::Nucleus nurex_target = nurex::get_default_nucleus(at,zt);

    double eout = energy_out(projectile,projectile.T, target,c);
    std::cout<<eout<<"\n";
    nurex::GlauberModelOLA_ZeroRange gm(nurex_projectile, nurex_target);
    double cs = nurex::SigmaR(gm, projectile.T);

    double rr = reaction_rate(cs,target.number_density_cm2(0));
    std::cout<<rr<<"\n";
    
    return 1.0;
    }

}

#endif
