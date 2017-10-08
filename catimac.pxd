"""
    catima cython
    ~~~~~~~~~~~~~~~~~
    :copyright: (c) 2017 by Andrej Prochazka
    :licence: GNU Affero General Public License, see LICENCE for more details
"""

from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "catima/structures.h" namespace "catima":
    cdef struct Target:
        double A
        int Z
    
    cdef struct Projectile:
        double A
        double Z
        double Q
        double T

    cdef struct Result:
        double Ein
        double Eout
        double Eloss
        double range
        double dEdxi
        double dEdxo
        double sigma_E
        double sigma_a
        double sigma_r
        double tof

    cdef cppclass MultiResult:
        vector[Result] results
        Result total_result

    cdef cppclass Material:
        Material() except +
        void add_element(double , int , double )
        pair[Target,double] get_element(int) 
        int ncomponents()
        double M()
        double density()
        void density(double val)
        double thickness()
        void thickness(double val)
    
    cdef cppclass Layers:
        Layers() except +
        const vector[Material]& get_materials() const
        void add(Material m)
        int num()const
        Material& operator[](int i)
        Layers& operator=(const Layers& other)

cdef extern from "catima/material_database.h" namespace "catima":
    cdef Material get_material(int)

cdef extern from "catima/config.h" namespace "catima":
    cdef struct Config:
        char z_effective;
        char skip;
        char dedx;

cdef extern from "catima/catima.h" namespace "catima":
    cdef double dedx(Projectile &p, double T, const Material &t,const Config &c)
    cdef double range(Projectile &p, double T, const Material &t, const Config &c);
    cdef double dedx_from_range(Projectile &p, double T, const Material &t, const Config &c);
    cdef double energy_out(Projectile &p, double T, const Material &t, const Config &c);

    cdef double domega2de(Projectile &p, double T, const Material &t, const Config &c);
    cdef double da2de(Projectile &p, double T, const Material &t, const Config &c);

    cdef double range_straggling(Projectile &p, double T, const Material &t, const Config &c);
    cdef double da2dx(Projectile &p, double T, const Material &t, const Config &c);
    cdef double dedx_rms(Projectile &p, Target &t, const Config &c);
    cdef double angular_variance(Projectile &p, double T, const Material &t, const  Config& c);

    cdef Result calculate(Projectile &p, const Material &t, const Config &c);
    cdef MultiResult calculate(Projectile &p, const Layers &layers, const Config &c);

cdef extern from "catima/calculations.h" namespace "catima":
    cdef double z_effective(const Projectile &p, const  Target &t, const  Config &c);
    cdef double z_eff_Pierce_Blann(double z, double beta);

cdef extern from "catima/constants.h" namespace "catima":        
    int max_datapoints "catima::max_datapoints"
    int logEmin "catima::logEmin"
    int logEmax "catima::logEmax"
    
cdef extern from "catima/storage.h" namespace "catima":
    cdef cppclass Interpolator:
        Interpolator(const double *x, const double *y, int num) except +
        double eval(double)
        double derivative(double)
        
    cdef cppclass DataPoint:
        vector[double] range
        vector[double] range_straggling
        vector[double] angular_variance
    
    cdef cppclass EnergyTableType "catima::EnergyTable[max_datapoints]":
        size_t num;
        double operator()(int i)
        
        
    cdef EnergyTableType energy_table;        
    cdef DataPoint& get_data(const Projectile &p, const Material &t, Config c);
