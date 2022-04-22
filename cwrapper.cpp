#include "catima/cwrapper.h"
#include "catima/catima.h"
#include "catima/material_database.h"
#include "catima/nucdata.h"
#include "catima/reactions.h"
#include <cstring>

extern "C" {
    struct CatimaConfig catima_defaults = {1};

    catima::Material make_material(double ta, double tz, double thickness, double density){
        catima::Material mat;        
        if(tz>200){
            mat = catima::get_material(tz);
        }
        else{            
            mat.add_element(ta,tz,1.0);
        }
        mat.density(density).thickness(thickness);
        if(density<=0.0){
            catima::Material m0 = catima::get_material(tz);
            mat.density(m0.density());
        }    
        return mat;
    }

    CatimaResult catima_calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density){
        catima::default_config.z_effective = catima_defaults.z_effective;
        catima::Projectile p(pa,pz);
        catima::Material mat = make_material(ta,tz, thickness, density);
        catima::Result r =  catima::calculate(p(T),mat);    
        CatimaResult res;
        res.Ein = r.Ein;
        res.Eout = r.Eout;
        res.Eloss = r.Eloss;
        res.range = r.range;
        res.dEdxi = r.dEdxi;
        res.dEdxo = r.dEdxo;
        res.sigma_E = r.sigma_E;
        res.sigma_a = r.sigma_a;
        res.sigma_r = r.sigma_r;
        res.tof = r.tof;
        // printf("%d\n",catima::_storage.get_index());
        return res;
    }

    double catima_Eout(double pa, int pz, double T, double ta, double tz, double thickness, double density){
        catima::default_config.z_effective = catima_defaults.z_effective;
        catima::Projectile p(pa,pz);
        catima::Material mat = make_material(ta,tz, thickness, density);          
        return energy_out(p(T), mat);
    }

    double catima_range(double pa, int pz, double T, double ta, double tz){
        catima::default_config.z_effective = catima_defaults.z_effective;
        catima::Projectile p(pa,pz);
        catima::Material mat = make_material(ta,tz, 0, -1);
        return range(p, mat);
    }

    double catima_range_straggling(double pa, int pz, double T, double ta, double tz){
        catima::default_config.z_effective = catima_defaults.z_effective;
        catima::Projectile p(pa,pz);
        catima::Material mat = make_material(ta,tz, 0, -1);
        return range_straggling(p, T, mat);
    }

    double catima_angular_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz){
        catima::Projectile p(pa,pz);        
        catima::Material mat = make_material(ta,tz, 0, -1);

        return catima::angular_straggling_from_E(p(Tin),Tout,mat);
    }

    double catima_energy_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz){
        catima::Projectile p(pa,pz);
        catima::Material mat = make_material(ta,tz, 0, -1);
        return catima::energy_straggling_from_E(p,Tin,Tout,mat);
    }

    double atomic_weight(int i){return catima::element_atomic_weight(i);}

    double catima_nonreaction_rate(double pa, int pz, double T, double ta, double tz, double thickness){
        catima::Projectile p(pa,pz);
        p.T = T;
        catima::Material mat = make_material(ta,tz, thickness, -1);
        return catima::nonreaction_rate(p,mat);
    }

}