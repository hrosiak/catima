#include <iostream>
#include <cmath>
#include <algorithm>
#include "catima/catima.h"
#include "catima/constants.h"
#include "catima/data_ionisation_potential.h"
#include "catima/data_atima.h"
#include "catima/integrator.h"
#include "catima/storage.h"
#include "catima/nucdata.h"
#include "catima/calculations.h"
#ifdef REACTIONS
#include "catima/reactions.h"
#endif

namespace catima{

Config default_config;

bool operator==(const Config &a, const Config&b){
    if(a.scattering != b.scattering)return false;
    if(a.z_effective != b.z_effective)return false;
    if(a.corrections != b.corrections)return false;
    if(a.calculation != b.calculation)return false;
    if(a.low_energy != b.low_energy)return false;
    if(fabs(a.scattering_factor-b.scattering_factor)>1e-3) return false;
    return true;
    }


double dedx(const Projectile &p, const Material &mat, const  Config &c){
    double sum = 0;
    if(p.T<=0)return 0.0;
    sum += dedx_n(p,mat);
    double se=0;
    if(p.T<=10){
        se = sezi_dedx_e(p,mat,c );
    }
    else if(p.T>10 && p.T<30){
        double factor = 0.05 * ( p.T - 10.0 );
        se = (1-factor)*sezi_dedx_e(p,mat,c) + factor*bethek_dedx_e(p,mat,c);
    }
    else {
        se = bethek_dedx_e(p,mat,c);
    }
    sum+=se;

    return sum;
}

double domega2dx(const Projectile &p, const Material &mat, const Config &c){
    double sum = 0;
    
    for(int i=0;i<mat.ncomponents();i++){
        auto t= mat.get_element(i);
        double w = mat.weight_fraction(i);
        sum += w*dedx_variance(p,t,c);
    }
    return sum;
}

double range(const Projectile &p, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);
    return range_spline(p.T);
}

double dedx_from_range(const Projectile &p, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);
    return p.A/range_spline.derivative(p.T);
}

std::vector<double> dedx_from_range(const Projectile &p, const std::vector<double> &T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);
    std::vector<double> dedx;
    dedx.reserve(T.size());
    for(auto e:T){
        if(e<catima::Ezero){
            dedx.push_back(0.0);
        }
        else{
            dedx.push_back(p.A/range_spline.derivative(e));
        }
    }
    return dedx;
}

double range_straggling(const Projectile &p, double T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    spline_type range_straggling_spline = get_range_straggling_spline(data);
    return sqrt(range_straggling_spline(T));
}

double range_variance(const Projectile &p, double T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    spline_type range_straggling_spline = get_range_straggling_spline(data);
    return range_straggling_spline(T);
}

double domega2de(const Projectile &p, double T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_straggling_spline(energy_table.values,data.range_straggling.data(),energy_table.num);
    spline_type range_straggling_spline = get_range_straggling_spline(data);
    return range_straggling_spline.derivative(T);
}
double da2dx(const Projectile &p, const Material &mat, const Config &c){  
    double Es2 = 198.81;
    //double f = 1.0;    
    if(c.scattering == scattering_types::dhighland)Es2 = Es_squared;
    if(c.scattering == scattering_types::fermi_rossi)Es2 = Es_squared;
    return c.scattering_factor*angular_scattering_power(p,mat, Es2);
}

/*
double da2de(const Projectile &p, double T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator angular_variance_spline(energy_table.values,data.angular_variance.data(),energy_table.num);
    spline_type angular_variance_spline = get_angular_variance_spline(data);
    return angular_variance_spline.derivative(T);
}
*/

double da2de(const Projectile &p, const Material &mat, const Config &c){
    return da2dx(p,mat,c)/dedx(p,mat,c);    
}

double Tfr(const Projectile &p, double X, double Es2){
    if(p.T<=0)return 0.0;
    double _p = p_from_T(p.T,p.A);
    double beta = _p/((p.T+atomic_mass_unit)*p.A);
    return Es2 /(X*ipow(_p*beta,2));
}

double angular_variance(Projectile p, const Material &t, const Config &c, int order){
    const double T = p.T;    
    const double p1 = p_from_T(T,p.A);
    const double beta1 = p1/((T+atomic_mass_unit)*p.A);    
    assert(T>0.0);
    assert(t.density()>0.0);
    assert(t.thickness()>0.0);    
    auto& data = _storage.Get(p,t,c);    
    spline_type range_spline = get_range_spline(data);    
    double range = range_spline(T);    
    double rrange = std::min(range/t.density(), t.thickness_cm()); // residual range, in case of stopping inside material
    double X0 = radiation_length(t);
    double Es2 = 198.81;
    if(c.scattering == scattering_types::dhighland)Es2 = Es_squared;
    if(c.scattering == scattering_types::fermi_rossi)Es2 = Es_squared;

    auto fx0p = [&](double x)->double{         
        double e =energy_out(T,x*t.density(),range_spline);
        double d = ipow((rrange-x),order);
        double ff = 1;        
        if(c.scattering == scattering_types::dhighland){
            double l = x*t.density()/X0;
            double lnl = log(l);
            ff = 0.97*(1+lnl/20.7)*(1+lnl/22.7);
        }
	    //return d*ff*da2dx(p(e), t, c);
        return d*ff*Tfr(p(e),1, 1.0);
            };

    auto fx0p_2 = [&](double x)->double{         
        double e =energy_out(T,x*t.density(),range_spline);                
        double d = ipow((rrange-x),order);
        return d*angular_scattering_power_xs(p(e),t,p1,beta1);
            };
    
    // corrections      
    if(c.scattering == scattering_types::gottschalk){
        return integrator.integrate(fx0p_2,0, rrange)*t.density();
    }
    return integrator.integrate(fx0p,0, rrange)*t.density()*ipow(p.Z,2)*Es2/X0;
}

double angular_straggling(Projectile p, const Material &t, const Config &c){
    return sqrt(angular_variance(p,t,c));
}

double angular_straggling_from_E(const Projectile &p, double Tout, Material t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    spline_type range_spline = get_range_spline(data);    
    double th = range_spline(p.T)-range_spline(Tout);    
    t.thickness(th);
    return angular_straggling(p,t,c);
}

double energy_straggling_from_E(const Projectile &p, double T, double Tout,const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);    
    spline_type range_spline = get_range_spline(data);
    spline_type range_straggling_spline = get_range_straggling_spline(data);
    double dEdxo = p.A/range_spline.derivative(Tout);
    return dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(Tout))/p.A;
}

double energy_out(double T, double thickness, const Interpolator &range_spline){
    int counter = 0;
    double range;
    double dedx;
    double e,r;
    
    range = range_spline(T);
    dedx = 1.0/range_spline.derivative(T);
    if(range<= thickness) return 0.0;
    
    e = T - (thickness*dedx);
    while(1){
        r = range - range_spline(e) - thickness;
        if(fabs(r)<Eout_th_epsilon)return e;
        double step = -r*dedx;
        e = e-step;
        if(e<Ezero)return 0.0;
        dedx = 1.0/range_spline.derivative(e);
        counter++;
        assert(counter<=100);
        if(counter>100)return -1;
    }
    return -1;
}

double energy_out(const Projectile &p, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);
    return energy_out(p.T,t.thickness(),range_spline);
    }

std::vector<double> energy_out(const Projectile &p, const std::vector<double> &T, const Material &t, const Config &c){
    auto& data = _storage.Get(p,t,c);
    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);

    std::vector<double> eout;
    eout.reserve(T.size());
    for(auto e:T){
        if(e<catima::Ezero){
            eout.push_back(0.0);
        }
        else{
            eout.push_back(energy_out(e,t.thickness(),range_spline));
        }
    }
    
    return eout;
    }

std::vector<double> calculate_tof(Projectile p, const Material &t, const Config &c){
    double res;
    std::vector<double> values;
    values.reserve(max_datapoints);
    auto function = [&](double x)->double{return 1.0/(dedx(p(x),t,c)*beta_from_T(x));};
    res = integrator.integrate(function,Ezero,energy_table(0));
    res = res*10.0*p.A/(c_light*t.density());
    values.push_back(res);
    for(int i=1;i<max_datapoints;i++){
        res = integrator.integrate(function,energy_table(i-1),energy_table(i));
        res = res*10.0*p.A/(c_light*t.density());
        res += values[i-1];
        values.push_back(res);
    }
    return values;
}

Result calculate(Projectile p, const Material &t, const Config &c){
    Result res;
    double T = p.T;
    if(T<catima::Ezero && T<catima::Ezero-catima::numeric_epsilon){return res;}
    auto& data = _storage.Get(p,t,c);

    bool use_angular_spline = false;
    if(c.scattering == scattering_types::atima_scattering){
        use_angular_spline = true;
    }

    //Interpolator range_spline(energy_table.values,data.range.data(),energy_table.num);
    spline_type range_spline = get_range_spline(data);
    spline_type range_straggling_spline = get_range_straggling_spline(data);    

    res.Ein = T;
    res.range = range_spline(T);
    res.dEdxi = p.A/range_spline.derivative(T);   
    res.sigma_r = sqrt(range_straggling_spline(T));        

    if(t.thickness()==0){
        res.dEdxo = res.dEdxi;
	res.Eout = res.Ein;
        return res;
    }
    
    res.Eout = energy_out(T,t.thickness(),range_spline);
    res.Eloss = (res.Ein - res.Eout)*p.A;
        
    if(res.Eout<Ezero){
        res.dEdxo = 0.0;
        res.sigma_a = 0.0;
        res.tof = 0.0;
        res.sigma_E = 0.0;
    }    
    else{        
        spline_type angular_variance_spline = get_angular_variance_spline(data);
        res.dEdxo = p.A/range_spline.derivative(res.Eout);        
        #ifdef THIN_TARGET_APPROXIMATION
        if(thin_target_limit*res.Ein<res.Eout){
            double edif = (res.Ein-res.Eout);            
            double s1 = range_straggling_spline.derivative(T);
            double s2 = range_straggling_spline.derivative(res.Eout);
            res.sigma_E = res.dEdxo*sqrt(edif*0.5*(s1+s2))/p.A;
            if(use_angular_spline){
                s1 = angular_variance_spline.derivative(T);
                s2 = angular_variance_spline.derivative(res.Eout);
                res.sigma_a = sqrt(0.5*(s1+s2)*edif);
            }        
        }
        else{
            res.sigma_E = res.dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(res.Eout))/p.A;            
            if(use_angular_spline){
                res.sigma_a = sqrt(angular_variance_spline(T) - angular_variance_spline(res.Eout));
            }
            
        }        
        #else
        res.sigma_E = res.dEdxo*sqrt(range_straggling_spline(T) - range_straggling_spline(res.Eout))/p.A;
        if(use_angular_spline){
            res.sigma_a = sqrt(angular_variance_spline(T) - angular_variance_spline(res.Eout));
        }        
        #endif
        if( (!use_angular_spline) && res.range>t.thickness()){ // do not calculate angle scattering when stopped inside material        
            res.sigma_a = angular_straggling(p(T),t,c);
        }
        //Interpolator tof_spline(energy_table.values, tofdata.data(), energy_table.num,interpolation_t::linear);
        //res.tof = tof_spline(res.Ein) - tof_spline(res.Eout);
        res.tof = calculate_tof_from_E(p,res.Eout,t);
            
    } //end of else for non stopped case
    
    // position straggling in material    
    double rrange = std::min(res.range/t.density(), t.thickness_cm());   
    res.sigma_x = angular_variance(p(T),t,c,2);
    res.sigma_x = sqrt(res.sigma_x);

    rrange = std::min(res.range/t.density(), t.thickness_cm());
    // position vs angle covariance, needed later for final position straggling    
    res.cov = angular_variance(p(T),t,c,1);

    #ifdef REACTIONS
    res.sp = nonreaction_rate(p,t,c);
    #endif
    return res;
}

MultiResult calculate(const Projectile &p, const Phasespace &ps, const Layers &layers, const Config &c){
    MultiResult res;
    double e = p.T;
    res.total_result.Ein = e;
    res.total_result.sigma_a = ps.sigma_a*ps.sigma_a;
    res.total_result.sigma_x = ps.sigma_x*ps.sigma_x;
    res.total_result.cov = ps.cov_x;
    res.results.reserve(layers.num());
    for(auto&m:layers.get_materials()){
        Result r = calculate(p,m,e,c);
        e = r.Eout;        
        res.total_result.Eloss += r.Eloss;
        res.total_result.sigma_E += r.sigma_E*r.sigma_E; 
        res.total_result.tof += r.tof;
        res.total_result.Eout = r.Eout;        
        double a2 = res.total_result.sigma_a;
        res.total_result.sigma_x += (2*m.thickness_cm()*res.total_result.cov) 
                                 + (a2*m.thickness_cm()*m.thickness_cm()) 
                                 + r.sigma_x*r.sigma_x;
        res.total_result.cov += a2*m.thickness_cm() + r.cov;
        res.total_result.sigma_a += r.sigma_a*r.sigma_a;
        #ifdef REACTIONS
        res.total_result.sp = (r.sp>=0.0)?res.total_result.sp*r.sp:-1;
        #endif
        res.results.push_back(r);
    }
    if(e>Ezero){
        res.total_result.sigma_a = sqrt(res.total_result.sigma_a);
        res.total_result.sigma_E = sqrt(res.total_result.sigma_E);
        res.total_result.sigma_x = sqrt(std::abs(res.total_result.sigma_x));
        
    }
    else{
        res.total_result.sigma_a = 0.0;
        res.total_result.sigma_E = 0.0;
        res.total_result.sigma_x = sqrt(std::abs(res.total_result.sigma_x));
        }
    return res;
}

Result calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density){
    Projectile p(pa,pz);
    Material m(ta,tz,density,thickness);
    return calculate(p(T),m);
}

DataPoint calculate_DataPoint(Projectile p, const Material &t, const Config &c){
    DataPoint dp(p,t,c);
    dp.range.resize(max_datapoints);
    dp.range_straggling.resize(max_datapoints);
    dp.angular_variance.resize(max_datapoints);
    auto fdedx = [&](double x)->double{
            return 1.0/dedx(p(x),t,c);
            };
    auto fomega = [&](double x)->double{
            return domega2dx(p(x),t,c)/catima::power(dedx(p(x),t,c),3);
            };
    auto ftheta = [&](double x)->double{
          return da2de(p(x),t,c);
          };

    //double res=0.0;
    //calculate 1st point to have i-1 element ready for loop
    //res = integrator.integrate(fdedx,Ezero,energy_table(0));
    //res = p.A*res;
    //dp.range[0] = res;
    
    dp.range[0] = 0.0;
    dp.angular_variance[0] = 0.0;

    //res = integrator.integrate(fomega,Ezero,energy_table(0));
    //res = p.A*res;
    dp.range_straggling[0]=0.0;
    //p.T = energy_table(0);    
    for(int i=1;i<max_datapoints;i++){
        double res = p.A*integrator.integrate(fdedx,energy_table(i-1),energy_table(i));
        dp.range[i] = res + dp.range[i-1];
        //res = da2dx(p(energy_table(i)),t)*res;
        //dp.angular_variance[i] = res + dp.angular_variance[i-1];        
        dp.angular_variance[i] = p.A*integrator.integrate(ftheta,energy_table(i-1),energy_table(i))
                                + dp.angular_variance[i-1];

        res = integrator.integrate(fomega,energy_table(i-1),energy_table(i));
        res = p.A*res;
        dp.range_straggling[i] = res + dp.range_straggling[i-1];
    }
    return dp;
}

double calculate_tof_from_E(Projectile p, double Eout, const Material &t, const Config &c){
    double res;
    auto function = [&](double x)->double{return 1.0/(dedx(p(x),t,c)*beta_from_T(x));};
    res = integrator.integrate(function,Eout,p.T);
    res = res*10.0*p.A/(c_light*t.density());
    return res;
}

std::pair<double,double> w_magnification(const Projectile &p, double Ein, const Material &t, const Config &c){
    std::pair<double, double> res{1.0,1.0};
    if(t.density()<= 0.0 || t.thickness()<=0){
        return res;
    }
    std::vector<double> energies{0.99*Ein, Ein, 1.01*Ein};
    auto eres = energy_out(p,energies,t,c);
    if(eres[0]>0.0 && eres[1]>0.0 && eres[2]>0.0){
        res.first = energies[1]*(eres[2]-eres[0])/(eres[1]*(energies[2]-energies[0]));
        res.second = p_from_T(energies[1],p.A)*(p_from_T(eres[2],p.A)-p_from_T(eres[0],p.A))/( p_from_T(eres[1],p.A)*( p_from_T(energies[2],p.A)-p_from_T(energies[0],p.A) ) );
        }
    else {
        res.first = 0.0;
        res.second = 0.0;
    }
    return res;
}

} // end of atima namespace
