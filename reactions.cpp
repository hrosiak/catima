#include "catima/reactions.h"
#ifdef NUREX
#include <cmath>
namespace catima{
    
double reaction_rate(Projectile &projectile, const Material &target, const Config &c){
    int num_elements = target.ncomponents();
    int ap = lround(projectile.A);
    int zp = lround(projectile.Z);
    nurex::Nucleus nurex_projectile = get_default_nucleus(ap,zp);
    
    int zt = lround(target.get_element(0).Z;
    int at = abundance::get_isotope_a(zt,0);
    nurex::Nucleus nurex_target = get_default_nucleus(at,zt);
    
    nurex::GlauberModelOLA_ZeroRange gm(nurex_projectile, nurex_target);
    double cs = nurex::SigmaR(gm, projectile.T);
    
    
    return 1.0;
    }

}

#endif
