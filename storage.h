/*
 *  Author: Andrej Prochazka
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

#ifndef STORAGE
#define STORAGE

#include <vector>
#include <iterator>
//#include <unordered_set>
#include <gsl/gsl_spline.h>
#include "catima/constants.h"
#include "catima/catima.h"


namespace catima{

	enum interpolation_t {cspline, linear};

    //inline double energy_function( int i ) { return exp(M_LN10*(logEmin + ((double)i)*(logEmax-logEmin)/(max_datapoints - 1.0))); }
  
    enum DataType{TYPE_RANGE,TYPE_LS};
    
    template<int N>
    struct EnergyTable{
	constexpr EnergyTable(double logmin, double logmax):values(),step(0.0),num(N){
		step = (logmax-logmin)/(N - 1.0);
	    for(auto i=0;i<N;i++){
		values[i]=exp(M_LN10*(logmin + ((double)i)*step));
		}
	    }
	double operator()(int i)const{return values[i];}
	double values[N];
	double step;
	std::size_t num;
    };

	template<int N>
	int EnergyTable_index(const EnergyTable<N> &table, double val){
		val = log(val)/M_LN10;
		if(val<table.values[0] || val>table.values[table.num-1])return -1;
		int i = (int)val/table.step;
		//double x = 1.0 - ((val - table.valuesp[i])/table.step);
		return i;
	}
	
	template<int N>
	double EnergyTable_interpolate(const EnergyTable<N> &table, double xval, double *y){
	    double r;
	    double lxval = log(xval)/M_LN10;
	    if(xval<table.values[0] || xval>table.values[table.num-1])return 0.0;
	    if(xval==table.values[table.num-1])return y[table.num-1];
	    int i = (int)(lxval/table.step);
		double linstep = table.values[i+1] - table.values[i];
	    double x = 1.0 - ((xval - table.values[i])/linstep);
	    r = (x*y[i]) + ((1-x)*y[i+1]);
	    return r;
	}
    /*
    template<int N>
    struct EnergyTableLinear{
	constexpr EnergyTableLinear():values(),num(N){
	    for(auto i=0;i<N;i++){
		values[i]=exp(M_LN10*(logEmin + ((double)i)*(logEmax-logEmin)/(N - 1.0)));
		}
	    }
	double operator()(int i)const{return values[i];}
	double values[N];
	std::size_t num;
    };
    */
    constexpr EnergyTable<max_datapoints> energy_table(logEmin,logEmax);
  
    class DataPoint{
	public:
	Projectile p;
	Material m;
	Config config;
	std::vector<double> range;
	std::vector<double> range_straggling;
	std::vector<double> angular_variance;

	DataPoint(){};
	DataPoint(const Projectile _p, const Material _m,Config _c=default_config):p(_p),m(_m),config(_c){};
	~DataPoint();
	friend bool operator==(const DataPoint &a, const DataPoint &b);
    };
    
    class Data{
	public:
	Data();
	~Data();
	void Add(const Projectile &p, const Material &t, Config c=default_config);
	int GetN() const {return storage.size();};
	void Reset(){storage.clear();storage.resize(max_storage_data);};
	DataPoint& Get(const Projectile &p, const Material &t, Config c=default_config);
	int get_index() {return std::distance(storage.begin(),index);}
	private:
	std::vector<DataPoint> storage;
	std::vector<DataPoint>::iterator index;
    };
    
/// Interpolation class, to store interpolated values
	class Interpolator{
        public:
        Interpolator(const double *x, const double *y, int num,interpolation_t type=cspline);
        Interpolator(const std::vector<double>& x, const std::vector<double>& y,interpolation_t type=cspline);
        ~Interpolator();
        double operator()(double x){return eval(x);};
        double eval(double x);
        double derivative(double x);
        double get_min(){return min;};
        double get_max(){return max;};

        private:
        double min=0;
        double max=0;
        gsl_interp_accel *acc;
        gsl_spline *spline;
    };
    
    extern Data _storage;
	
	inline DataPoint& get_data(const Projectile &p, const Material &t, Config c=default_config){
		return _storage.Get(p, t, c);
	}

    bool operator==(const DataPoint &a, const DataPoint &b);
}

#endif
