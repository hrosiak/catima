#include "catima/cwrapper.h"
#include "catima/catima.h"
#include "catima/material_database.h"
#include <cstring>

extern "C" {

    CatimaResult catima_calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density){
        catima::default_config.z_effective = catima_defaults.z_effective;
        catima::Material mat;
        catima::Projectile p(pa,pz);
        if(tz>200){
            mat = catima::get_material(tz);
        }
        else{
            mat.add_element(ta,tz,1.0);
        }
        mat.density(density).thickness(thickness);
        catima::Result r =  catima::calculate(p(T),mat);    
        CatimaResult res;
        std::memcpy(&res,&r,sizeof(res));
        return res;
    }

    double catima_angular_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz){
        catima::Projectile p(pa,pz);
        
        catima::Material mat;
        if(tz>200){
            mat = catima::get_material(tz);
        }
        else{
            mat.add_element(ta,tz,1.0);
        }

        return catima::angular_straggling_from_E(p,Tin,Tout,mat);
    }

    double catima_energy_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz){
        catima::Projectile p(pa,pz);
        catima::Material mat;
        if(tz>200){
            mat = catima::get_material(tz);
        }
        else{
            mat.add_element(ta,tz,1.0);
        }
        return catima::energy_straggling_from_E(p,Tin,Tout,mat);
    }

}