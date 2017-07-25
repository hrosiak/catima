#include <math.h>
#include <iostream>
#include "storage.h"
namespace catima {
    Data _storage;
    
    
    bool operator==(const DataPoint &a, const DataPoint &b){
	if( (a.m == b.m) && (a.p == b.p) && (a.config == b.config)){
	    return true;
	}
	else{
	    return false;
        }
    }
    
    DataPoint::~DataPoint(){
	
    }
    
    Data::Data(){
        //storage.reserve(max_storage_data);
        storage.resize(max_storage_data);
        index = storage.begin();
    }
    
    Data::~Data(){
    }
    

void Data::Add(const Projectile &p, const Material &t, Config c){
	DataPoint dp(p,t,c);
	for(auto &e:storage){
	    if(e==dp)return; 
	}
    if(index==storage.end())index=storage.begin();
    #if(0)
    *index = dp;
    index->range = calculate_range(p,t,c);
	index->range_straggling = calculate_range_straggling(p,t,c);
	index->angular_variance = calculate_angular_variance(p,t,c);
    #else
    *index = calculate_DataPoint(p,t,c);
    #endif

    index++;
    }
    
DataPoint& Data::Get(const Projectile &p, const Material &t, Config c){
	for(auto &e:storage){
	    if( (e.p==p) && (e.m==t) && (e.config==c)){
		return e;
	    }
	}
    Add(p,t,c);
	//return storage.back();
    return *std::prev(index);
    }

//////////// Interpolator ////////////////////////////////
Interpolator::Interpolator(const double *x, const double *y, int num,interpolation_t type){
    acc = gsl_interp_accel_alloc ();

    if(type==cspline)
	spline = gsl_spline_alloc (gsl_interp_cspline, num);
    else
	spline = gsl_spline_alloc (gsl_interp_linear, num);

    gsl_spline_init (spline, x, y, num);
    min= x[0];
    max= x[num-1];

}
Interpolator::Interpolator(const std::vector<double>& x, const std::vector<double>& y,interpolation_t type){
    //Interpolator(x.data(),y.data(),x.size());
    acc = gsl_interp_accel_alloc ();
    if(type==cspline)
	spline = gsl_spline_alloc (gsl_interp_cspline, x.size());
    else
	spline = gsl_spline_alloc (gsl_interp_linear, x.size());

    gsl_spline_init (spline, x.data(), y.data(), x.size());
    min= x[0];
    max= x[x.size()-1];
}

Interpolator::~Interpolator(){
    gsl_interp_accel_free (acc);
    gsl_spline_free (spline);
}

double Interpolator::eval(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval(spline, x, acc);
}

double Interpolator::derivative(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval_deriv (spline, x, acc);
}

}
