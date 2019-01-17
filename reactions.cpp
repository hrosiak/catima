#include "catima/reactions.h"
#include "catima/catima.h"
#include "catima/abundance_database.h"
#include "catima/storage.h"
#include <cmath>
#include <iostream>

#ifdef NUREX
#include "nurex/Parametrization.h"
using nurex::SigmaR_Kox;
#else
using catima::SigmaR_Kox;
#endif


namespace catima{
    
double nonreaction_rate(Projectile &projectile, const Material &target, const Config &c){

    if(projectile.T<emin_reaction)return -1.0;
    if(target.thickness()<=0.0)return 1.0;

    int ap = lround(projectile.A);
    int zp = lround(projectile.Z);

    auto& data = _storage.Get(projectile,target,c);
    spline_type range_spline = get_range_spline(data);
    if(energy_out(projectile.T, target.thickness(), range_spline) < emin_reaction)return -1.0;
    
    auto sigma_r = [&](double th){
        double stn_sum=0.0, sum=0.0;
        double e = energy_out(projectile.T, th, range_spline);
        for(unsigned int i = 0;i<target.ncomponents();i++){
            int zt = target.get_element(i).Z;
            int at = abundance::get_isotope_a(zt,0); // most abundand natural isotope mass
            stn_sum += target.molar_fraction(i);
            sum += target.molar_fraction(i)*SigmaR_Kox(ap,zp,e,at,zt); 
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
        cs = Avogadro*catima::integrator.integrate(sigma_r,0,target.thickness())/target.M();
    }
    return exp(-cs*0.0001);
    }
    
double production_rate(double cs, double rcs_projectile, double rcs_product, const Material &target, const Config &c){
    double t = target.number_density_cm2();
    double res = 0.0;
    if( std::abs(rcs_product - rcs_projectile)<5.0){
        res = 0.0001*cs*t*(exp(-0.0001*rcs_product*t));
    }
    else{
        res = cs*(exp(-0.0001*rcs_projectile*t) - exp(-0.0001*rcs_product*t))/(rcs_product-rcs_projectile);
    }
    return res;
}

#ifndef NUREX
double SigmaR_Kox(int Ap, int Zp, double E, int At, int Zt){
    constexpr double rc = 1.3;
    constexpr double r0 = 1.1;
    constexpr double a = 1.85;
    constexpr double c1 = 2-(10.0/(1.5*1.5*1.5*1.5*1.5));
    double Ap13 = pow(Ap,1.0/3.0);
    double At13 = pow(At,1.0/3.0);
    double D = 5.0*(At-2*Zt)*Zp/(Ap*At);
    double Bc = Zp*Zt/(rc*(Ap13+At13));
    double logE = std::log10(E);
    double c = 0;
    if(logE < 1.5){
        c = c1*std::pow(logE/1.5,3);
    }
    else{
        c = (-10.0/std::pow(logE,5)) + 2;
    }
    double Rs = r0 * ((a*Ap13*At13)/(Ap13+At13)-c)+D;
    double Rv = r0 * (Ap13 + At13);
    double Ri = Rv + Rs;
    double Ecm = Ecm_from_T(E,Ap,At);
    return 10.0*PI*Ri*Ri*(1-(Bc/Ecm));
    }
#endif

} //end of catima namespace
