#include <iostream>
#include <math.h>
#include <algorithm>
#include "catima/catima.h"
#include "catima/constants.h"
#include "catima/data_ionisation_potential.h"
#include "catima/data_atima.h"
#include "catima/integrator.h"
#include "catima/storage.h"
#include "catima/nucdata.h"
#include "catima/calculations.h"

namespace catima{

Config default_config;


bool operator==(const Config &a, const Config&b){
    if(std::memcmp(&a,&b,sizeof(Config)) == 0){
        return true;
        }
    else
        return false;
    }


double dedx(Projectile &p, double T, const Material &mat, const  Config &c){
    double sum = 0;
    Target t;
    double w=0;
    if(T<=0)return 0.0;
    for(int i=0;i<mat.ncomponents();i++){
        auto res = mat.get_element(i);
        t = res.first;  //struct of target 
        w = res.second; //number of atoms of the element
        p.T = T;
        sum += t.A*w*dedx(p,t);
    }
    return sum/mat.M();
}

double domega2dx(Projectile &p, double T, const Material &mat, const Config &c){
    double sum = 0;
    Target t;
    double w=0;
    
    for(int i=0;i<mat.ncomponents();i++){
        auto res = mat.get_element(i);
        t = res.first;  //struct of target 
        w = res.second; //number of atoms of the element
        p.T = T;
        sum += t.A*w*dedx_rms(p,t);
    }
    return sum/mat.M();
}

double da2dx(Projectile &p, double T, const Material &mat, const Config &c){
    double sum = 0;
    Target t;
    double w=0;
    
    for(int i=0;i<mat.ncomponents();i++){
        auto res = mat.get_element(i);
        t = res.first;  //struct of target 
        w = res.second; //number of atoms of the element
        p.T = T;
        sum += t.A*w*angular_scattering_variance(p,t);
    }
    return sum/mat.M();
}


double range(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    return range_spline(T);
}

double dedx_from_range(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    return p.A/range_spline.derivative(T);
}

double range_straggling(Projectile &p, double T, const Material &t, const Config &c){
    double r=0;
    auto data = _storage.Get(p,t,c);
    Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    return sqrt(range_straggling_spline(T));
}

double range_variance(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    return range_straggling_spline(T);
}

double domega2de(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    return range_straggling_spline.derivative(T);
}

double da2de(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator angular_variance_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
    return angular_variance_spline.derivative(T);
}

double angular_straggling_from_E(Projectile &p, double T, double Tout, const Material &t, const Config &c){
    double r=0;
    auto data = _storage.Get(p,t,c);
    Interpolator angular_straggling_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
    return sqrt(angular_straggling_spline(T) - angular_straggling_spline(Tout));
}

double energy_straggling_from_E(Projectile &p, double T, double Tout,const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);

    Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    double dEdxo = p.A/range_spline.derivative(Tout);
    return dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(Tout))/p.A;
}

double energy_out(double T, double thickness, Interpolator &range_spline){
    constexpr double epsilon = 1E-5;
    int counter = 0;
    double lo=0,hi=T;
    double range;
    double dedx;
    double e,r;
    double step;
    
    range = range_spline(T);
    dedx = 1.0/range_spline.derivative(T);
    if(range<= thickness) return 0.0;
    
    e = T - (thickness*dedx);
    while(1){
        r = range - range_spline(e) - thickness;
        if(fabs(r)<epsilon)return e;
        step = -r*dedx;
        e = e-step;
        if(e<Ezero)return 0.0;
        dedx = 1.0/range_spline.derivative(T);
        counter++;
        if(counter>100){printf("too many iterations finding Eout");return -1;}
    }
    return -1;
}

double energy_out(Projectile &p, double T, const Material &t, const Config &c){
    auto data = _storage.Get(p,t,c);
    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    return energy_out(T,t.thickness(),range_spline);
    }

Result calculate(Projectile &p, const Material &t, const Config &c){
    Result res;
    double T = p.T;
    auto data = _storage.Get(p,t,c);

    Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    res.Ein = T;
    res.range = range_spline(T);
    res.dEdxi = p.A/range_spline.derivative(T);
    res.Eout = energy_out(T,t.thickness(),range_spline);

    Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);

    if(res.Eout<Ezero){
        res.dEdxo = 0.0;
        res.sigma_a = 0.0;
        res.tof = 0.0;
        res.sigma_E = 0.0;
    }
    else{
        res.dEdxo = p.A/range_spline.derivative(res.Eout);
        
        #ifdef THIN_TARGET_APPROXIMATION
        if(thin_target_limit*res.Ein<res.Eout){
            double edif = (res.Ein-res.Eout);
            double s1 = range_straggling_spline.derivative(T);
            double s2 = range_straggling_spline.derivative(res.Eout);
            res.sigma_E = res.dEdxo*sqrt(edif*0.5*(s1+s2))/p.A;
        
            Interpolator angular_variance_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
            s1 = angular_variance_spline.derivative(T);
            s2 = angular_variance_spline.derivative(res.Eout);
            res.sigma_a = sqrt(0.5*(s1+s2)*edif);
        }
        else{
            res.sigma_E = res.dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(res.Eout))/p.A;
            Interpolator angular_variance_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
            res.sigma_a = sqrt(angular_variance_spline(T) - angular_variance_spline(res.Eout));    
        }
        
        #else
        res.sigma_E = res.dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(res.Eout))/p.A;
        Interpolator angular_variance_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
        res.sigma_a = sqrt(angular_variance_spline(T) - angular_variance_spline(res.Eout));
        #endif
        if( !(c.skip&skip_tof) && t.thickness()>0){
            //auto tofdata = calculate_tof(p,t,c);
            //Interpolator tof_spline(energy_table.values, tofdata.data(), energy_table.num,interpolation_t::linear);
            //res.tof = tof_spline(res.Ein) - tof_spline(res.Eout);
            res.tof = calculate_tof_from_E(p,res.Eout,t);
            }
    }
    res.sigma_r = sqrt(range_straggling_spline(T));
    res.Eloss = (res.Ein - res.Eout)*p.A;
    return res;
}

MultiResult calculate(Projectile &p, const Layers &layers, const Config &c){
    MultiResult res;
    double e = p.T;
    res.total_result.Ein = e;
    res.results.reserve(layers.num());

    for(auto&m:layers.get_materials()){
        Result r = calculate(p,m,e,c);
        e = r.Eout;
        res.total_result.sigma_a += r.sigma_a*r.sigma_a;
        res.total_result.Eloss += r.Eloss;
        res.total_result.sigma_E += r.sigma_E*r.sigma_E; 
        res.total_result.tof += r.tof;
        res.total_result.Eout = r.Eout;
        res.results.push_back(r);
    }
    if(e>Ezero){
        res.total_result.sigma_a = sqrt(res.total_result.sigma_a);
        res.total_result.sigma_E = sqrt(res.total_result.sigma_E);
    }
    else{
        res.total_result.sigma_a = 0.0;
        res.total_result.sigma_E = 0.0;
        }
    return res;
}

Result calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density){
    Projectile p(pa,pz);
    Material m(ta,tz,density,thickness);
    return calculate(p(T),m);
}

std::vector<double> calculate_range(Projectile p, const Material &t, const Config &c){
    double res;
    std::vector<double>values;
    values.reserve(max_datapoints);
    auto fdedx = [&](double x)->double{return 1.0/dedx(p,x,t);};
    
    //calculate 1st point to have i-1 element ready for loop
    res = integratorGSL.Integrate(fdedx,Ezero,energy_table(0),int_eps_range,false);
    res = p.A*res;
    values.push_back(res);
    
    for(int i=1;i<max_datapoints;i++){
        res = integratorGSL.Integrate(fdedx,energy_table(i-1),energy_table(i),int_eps_range,false);
        res = p.A*res;
        res += values[i-1];
        values.push_back(res);
    }
    
    return values;
}


std::vector<double> calculate_range_straggling(Projectile p, const Material &t, const Config &c){
    double res;
    std::vector<double>values;
    values.reserve(max_datapoints);
    auto function = [&](double x)->double{return 1.0*domega2dx(p,x,t)/pow(dedx(p,x,t),3);};
    //auto function = [&](double x)->double{
                //double de = dedx(p,x,t);
                //return 1.0*domega2dx(p,x,t)/(de*de*de);
                //};
    //calculate 1st point to have i-1 element ready for loop
    res = integratorGSL.Integrate(function,Ezero,energy_table(0),int_eps_range_str,false);
    res = p.A*res;
    values.push_back(res);
    for(int i=1;i<max_datapoints;i++){
        res = integratorGSL.Integrate(function,energy_table(i-1),energy_table(i),int_eps_range_str,false);
        res = p.A*res;
        res += values[i-1];
        values.push_back(res);
    }
    
    return values;
}

std::vector<double> calculate_da2dx(Projectile p, const Material &t, const Config &c){
    double res;
    std::vector<double>values;
    values.reserve(max_datapoints);
    //auto function = [&](double x)->double{return p.A*da2dx(p,x,t)/dedx(p,x,t);};
    auto function = [&](double x)->double{return 1.0/dedx(p,x,t);};
    res = integratorGSL.Integrate(function,Ezero,energy_table(0),int_eps_ang_str,false);
    res = p.A*da2dx(p,energy_table(0),t)*res;
    values.push_back(res);
    for(int i=1;i<max_datapoints;i++){
        res = integratorGSL.Integrate(function,energy_table(i-1),energy_table(i),int_eps_ang_str,false);
        res = p.A*da2dx(p,energy_table(i),t)*res;
        res += values[i-1];
        values.push_back(res);
    }
    return values;
}

std::vector<double> calculate_tof(Projectile p, const Material &t, const Config &c){
    double res;
    std::vector<double> values;
    values.reserve(max_datapoints);
    auto function = [&](double x)->double{return 1.0/(dedx(p,x,t)*beta_from_T(x));};
    res = integratorGSL.Integrate(function,Ezero,energy_table(0),int_eps_tof,false);
    res = res*10.0*p.A/(c_light*t.density());
    values.push_back(res);
    for(int i=1;i<max_datapoints;i++){
        res = integratorGSL.Integrate(function,energy_table(i-1),energy_table(i),int_eps_tof,false);
        res = res*10.0*p.A/(c_light*t.density());
        res += values[i-1];
        values.push_back(res);
    }
    return values;
}

DataPoint calculate_DataPoint(Projectile p, const Material &t, const Config &c){
    DataPoint dp(p,t,c);
    std::vector<double>range_values;
    range_values.reserve(max_datapoints);
    std::vector<double>range_straggling_values;
    range_straggling_values.reserve(max_datapoints);
    std::vector<double>angular_straggling_values;
    angular_straggling_values.reserve(max_datapoints);

    double dedxval;
    auto fdedx = [&](double x)->double{
            return 1.0/dedx(p,x,t);    
            };
    auto fomega = [&](double x)->double{
            //return 1.0*domega2dx(p,x,t)/pow(dedx(p,x,t),3);
            return domega2dx(p,x,t)/catima::power(dedx(p,x,t),3);
            };

    double res;
    //calculate 1st point to have i-1 element ready for loop
    res = integratorGSL.Integrate(fdedx,Ezero,energy_table(0),int_eps_range,false);
    res = p.A*res;
    range_values.push_back(res);
    res = da2dx(p,energy_table(0),t)*res;
    angular_straggling_values.push_back(res);
    
    res = integratorGSL.Integrate(fomega,Ezero,energy_table(0),int_eps_range_str,false);
    res = p.A*res;
    range_straggling_values.push_back(res);

    for(int i=1;i<max_datapoints;i++){
        res = p.A*integratorGSL.Integrate(fdedx,energy_table(i-1),energy_table(i),int_eps_range,false);
        range_values.push_back(res + range_values[i-1]);

        res = da2dx(p,energy_table(i),t)*res;
        angular_straggling_values.push_back(res + angular_straggling_values[i-1]);
    
        res = integratorGSL.Integrate(fomega,energy_table(i-1),energy_table(i),int_eps_range_str,false);
        res = p.A*res;
        res += range_straggling_values[i-1];
        range_straggling_values.push_back(res);
    }

    dp.range = range_values;
    dp.range_straggling = range_straggling_values;
    dp.angular_variance = angular_straggling_values;
    return dp;
}

double calculate_tof_from_E(Projectile p, double Eout, const Material &t, const Config &c){
    double res;
    //double beta_in = beta_from_T(p.T);
    //double beta_out = beta_from_T(Eout);
    auto function = [&](double x)->double{return 1.0/(dedx(p,x,t)*beta_from_T(x));};
    res = integratorGSL.Integrate(function,Eout,p.T,int_eps_tof,false);
    res = res*10.0*p.A/(c_light*t.density());
    return res;
}

} // end of atima namespace
