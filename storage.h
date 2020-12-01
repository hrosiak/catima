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
#include <array>
#include <iterator>
#include <cmath>
//#include <unordered_set>
#include "catima/build_config.h"
#include "catima/constants.h"
#include "catima/structures.h"
#include "catima/config.h"

#include "catima/spline.h"


namespace catima{

    /**
     * Class to store energy points, log spaced from logmin to logmax.
     */
	template<int N>
    struct EnergyTable{
		EnergyTable(double logmin, double logmax):values(),step(0.0),num(N){
		step = (logmax-logmin)/(N - 1.0);
	    for(auto i=0;i<N;i++){
			values[i]=exp(LN10*(logmin + ((double)i)*step));
		}
	    }
	double operator()(int i)const{return values[i];}
    double operator[](int i)const{return values[i];}
    static constexpr int size() {return N;};
	double values[N];
	double step;
	double* begin(){return values;}
    double* end(){return &values[num];}
    int index(double v)const noexcept{
        double lxval = (log(v/values[0])/LN10);
        if(v<values[0] || step==0.0)return -1;
        if(v>=values[N-1]-numeric_epsilon)return N-1;
        int i = static_cast<int> (std::floor(lxval/step));
        if(v >= values[i+1]-numeric_epsilon)i++; // this is correction for floating point precision
        return i;
    };
	std::size_t num;
    };

	extern EnergyTable<max_datapoints> energy_table;

	template<int N>
	int EnergyTable_index(const EnergyTable<N> &table, double val){
		if(val<table.values[0] || val>table.values[table.num-1])return -1;
		double lxval = (log(val/table.values[0])/LN10);
        int i = static_cast<int>( std::floor(lxval/table.step));
        if(val >= table.values[i+1]-numeric_epsilon)i++; // this is correction for floating point precision
        return i;
	}

	template<int N>
	double EnergyTable_interpolate(const EnergyTable<N> &table, double xval, double *y){
	    double r;
	    if(xval<table.values[0] || xval>table.values[table.num-1])return 0.0;
	    if(xval==table.values[table.num-1])return y[table.num-1];
        int i = EnergyTable_index(table, xval);
	    double linstep = table.values[i+1] - table.values[i];
        if(linstep == 0.0)return table.values[i];
	    double x = 1.0 - ((xval - table.values[i])/linstep);
	    r = (x*y[i]) + ((1-x)*y[i+1]);
	    return r;
	}

    //////////////////////////////////////////////////////////////////////////////////////
    #ifdef GSL_INTERPOLATION
    /// Interpolation class, to store interpolated values
        class InterpolatorGSL{
            public:
            InterpolatorGSL(){};
            InterpolatorGSL(const EnergyTable<max_datapoints>& x, const std::vector<double>& y, interpolation_t type=cspline);
            ~InterpolatorGSL();
            double operator()(double x)const{return eval(x);};
            double eval(double x) const;
            double derivative(double x) const;
            double get_min()const{return min;};
            double get_max()const{return max;};

            private:
            double min=0;
            double max=0;
            gsl_interp_accel *acc;
            gsl_spline *spline;
        };
    #endif

    class InterpolatorCSpline{
    public:
        using xtype = EnergyTable<max_datapoints>;
        InterpolatorCSpline()=default;
        InterpolatorCSpline(const xtype &table, const std::vector<double> &y):
            min(table.values[0]), max(table.values[max_datapoints-1]), ss(table,y){}
        double operator()(double x)const{return eval(x);}
        double eval(double x)const{return ss.evaluate(x);}
        double derivative(double x)const{return ss.deriv(x);}
        double get_min()const{return min;}
        double get_max()const{return max;}

    private:
        double min=0;
        double max=0;
        cspline_special<xtype> ss;
    };

#ifdef GSL_INTERPOLATION
using Interpolator = InterpolatorGSL;
#else
using Interpolator = InterpolatorCSpline;
#endif

#ifdef STORE_SPLINES
    using spline_type = const Interpolator&;
#else
    using spline_type = Interpolator;
#endif

// return vector with lineary spaced elements from a to b, num is number of elements

/**
 * @brief structure to store calculated data points and optionally also splines
 */
class DataPoint{
	public:
	Projectile p;
	Material m;
	Config config;

    std::vector<double> range;
    std::vector<double> range_straggling;
    std::vector<double> angular_variance;
#ifdef STORE_SPLINES
    Interpolator range_spline;
    Interpolator range_straggling_spline;
    Interpolator angular_variance_spline;
#endif
    DataPoint()=default;
    DataPoint(const Projectile _p, const Material _m,const Config &_c=default_config):p(_p),m(_m),config(_c){}
    DataPoint(const DataPoint&)=delete;
    DataPoint(DataPoint&&)=default;
    DataPoint& operator=(const DataPoint&)=default;
    DataPoint& operator=(DataPoint&&)=default;
	friend bool operator==(const DataPoint &a, const DataPoint &b);
    };

#ifdef STORE_SPLINES
    const Interpolator& get_range_spline(const DataPoint &data);
    const Interpolator& get_range_straggling_spline(const DataPoint &data);
    const Interpolator& get_angular_variance_spline(const DataPoint &data);
#else
    Interpolator get_range_spline(const DataPoint &data);
    Interpolator get_range_straggling_spline(const DataPoint &data);
    Interpolator get_angular_variance_spline(const DataPoint &data);
#endif

/**
 * @brief The Data class to store DataPoints
 */
    class Data{
    public:
        Data();
        ~Data();

        /**
         * @brief Add new DataPoint
         * @param p - Projectile
         * @param t - Material
         * @param c - Config
         */
        void Add(const Projectile &p, const Material &t, const Config &c=default_config);

        int GetN() const {return storage.size();};
        void Reset(){storage.clear();storage.resize(max_storage_data);index=storage.begin();};

        /**
         * @brief Get DataPoint reference for projectile-target-config combination
         * @param p - Projectile
         * @param t - Material
         * @param c - Config
         * @return reference to DataPoint
         */
        DataPoint& Get(const Projectile &p, const Material &t, const Config &c=default_config);
        DataPoint& Get(unsigned int i){return storage[i];};
        int get_index() {return std::distance(storage.begin(),index);}

    private:
        std::vector<DataPoint> storage;
        std::vector<DataPoint>::iterator index;
    };

    extern Data _storage;

    /**
     * @brief get_data - Get DataPoint from the global storage class
     * @param p - Projectile
     * @param t - Material
     * @param c - Config
     * @return const reference to DataPoint
     */
    inline const DataPoint& get_data(const Projectile &p, const Material &t, const Config &c=default_config){
		return _storage.Get(p, t, c);
	}

    bool operator==(const DataPoint &a, const DataPoint &b);

}

#endif
