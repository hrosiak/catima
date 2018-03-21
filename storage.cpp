/*
 *  Copyright(C) 2017
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <iostream>
#include "storage.h"
#include "catima/catima.h"
namespace catima {
    Data _storage;
    EnergyTable<max_datapoints> energy_table(logEmin,logEmax);
    
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
        //storage.reserve(max_storage_data); // disabled because of "circular" storage
        storage.resize(max_storage_data);
        index = storage.begin();
    }
    
    Data::~Data(){
    }
    

void Data::Add(const Projectile &p, const Material &t, const Config &c){
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
    
DataPoint& Data::Get(const Projectile &p, const Material &t, const Config &c){
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
InterpolatorGSL::InterpolatorGSL(const double *x, const double *y, int num,interpolation_t type){
    acc = gsl_interp_accel_alloc ();

    if(type==cspline)
	spline = gsl_spline_alloc (gsl_interp_cspline, num);
    else
	spline = gsl_spline_alloc (gsl_interp_linear, num);

    gsl_spline_init (spline, x, y, num);
    min= x[0];
    max= x[num-1];

}
InterpolatorGSL::InterpolatorGSL(const std::vector<double>& x, const std::vector<double>& y,interpolation_t type){
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

InterpolatorGSL::~InterpolatorGSL(){
    gsl_interp_accel_free (acc);
    gsl_spline_free (spline);
}

double InterpolatorGSL::eval(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval(spline, x, acc);
}

double InterpolatorGSL::derivative(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval_deriv (spline, x, acc);
}


//////////// Interpolator2 ////////////////////////////////
#ifdef BUILTIN_SPLINE
Interpolator2::Interpolator2(const double *x, const double *y, int num){
    ss.set_points(x,y,num);
    min= x[0];
    max= x[num-1];

}

double Interpolator2::eval(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return ss(x);
}

double Interpolator2::derivative(double x){
    if(x<min)x=min;
    if(x>max)x=max;
    return ss.deriv(1,x);
}
#endif

}
