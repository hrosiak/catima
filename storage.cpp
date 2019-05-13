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

#ifdef GSL_INTERPOLATION
//////////// Interpolator ////////////////////////////////
InterpolatorGSL::InterpolatorGSL(const EnergyTable<max_datapoints>& x, const std::vector<double>& y, interpolation_t type){
    acc = gsl_interp_accel_alloc ();
    const int num = y.size();
    if(type==cspline)
    spline = gsl_spline_alloc (gsl_interp_cspline, num);
    else
    spline = gsl_spline_alloc (gsl_interp_linear, num);

    gsl_spline_init (spline, x.values, y.data(), num);
    min= x[0];
    max= x[num-1];

}

InterpolatorGSL::~InterpolatorGSL(){
    gsl_interp_accel_free (acc);
    gsl_spline_free (spline);
}

double InterpolatorGSL::eval(double x) const{
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval(spline, x, acc);
}

double InterpolatorGSL::derivative(double x)const{
    if(x<min)x=min;
    if(x>max)x=max;
    return gsl_spline_eval_deriv (spline, x, acc);
}
#endif

#ifdef STORE_SPLINES
    const Interpolator& get_range_spline(const DataPoint &data){
        return data.range_spline;
    }

    const Interpolator& get_range_straggling_spline(const DataPoint &data){
        return data.range_straggling_spline;
    }

    const Interpolator& get_angular_variance_spline(const DataPoint &data){
        return data.angular_variance_spline;
    }
#else
    Interpolator get_range_spline(const DataPoint &data){
        //return Interpolator(energy_table.values,data.range);
        //return data.range_spline;
        return Interpolator(energy_table,data.range);
    }

    Interpolator get_range_straggling_spline(const DataPoint &data){
        //return Interpolator(energy_table.values,data.range_straggling);
        //return data.range_straggling_spline;
        return Interpolator(energy_table,data.range_straggling);
    }

    Interpolator get_angular_variance_spline(const DataPoint &data){
        //return Interpolator(energy_table.values,data.angular_variance);
        //return data.angular_variance_spline;
        return Interpolator(energy_table,data.angular_variance);
    }
#endif
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
    *index = calculate_DataPoint(p,t,c);
#ifdef STORE_SPLINES
    //index->range_spline = Interpolator(energy_table.values,index->range);
    //index->range_straggling_spline = Interpolator(energy_table.values,index->range_straggling);
    //index->angular_variance_spline = Interpolator(energy_table.values,index->angular_variance);
    index->range_spline = Interpolator(energy_table, index->range);
    index->range_straggling_spline = Interpolator(energy_table, index->range_straggling);
    index->angular_variance_spline = Interpolator(energy_table, index->angular_variance);
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

}
