"""
    catima cython
    ~~~~~~~~~~~~~~~~~
    :copyright: (c) 2017 by Andrej Prochazka
    :licence: GNU Affero General Public License, see LICENCE for more details
"""

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp cimport bool

cdef extern from "catima/structures.h" namespace "catima":
    cdef struct Target:
        double A
        int Z
        double stn
    
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
        double sp

    cdef cppclass MultiResult:
        vector[Result] results
        Result total_result

    cdef cppclass Material:
        Material() except +
        void add_element(double , int , double )
        Target get_element(int) 
        int ncomponents()
        double M()
        double density()
        void density(double val)
        double thickness()
        void thickness(double val)
        void thickness_cm(double val)
        void calculate()
        double I()
        void I(double val)
    
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
        char dedx_straggling

cdef extern from "catima/catima.h" namespace "catima":
    cdef double dedx(Projectile &p, double T, const Material &t,const Config &c)
    cdef double domega2dx(Projectile &p, double T, const Material &mat, const Config &c)
    cdef double range(Projectile &p, double T, const Material &t, const Config &c);
    cdef double dedx_from_range(Projectile &p, double T, const Material &t, const Config &c);
    cdef vector[double] dedx_from_range(Projectile &p, vector[double] &T, const Material &t, const Config &c);
    cdef double energy_out(Projectile &p, double T, const Material &t, const Config &c);
    cdef vector[double] energy_out(Projectile &p, vector[double] &T, const Material &t, const Config &c);

    cdef double domega2de(Projectile &p, double T, const Material &t, const Config &c);
    cdef double da2de(Projectile &p, double T, const Material &t, const Config &c);

    cdef double range_straggling(Projectile &p, double T, const Material &t, const Config &c);
    cdef double da2dx(Projectile &p, double T, const Material &t, const Config &c);
    cdef double angular_variance(Projectile &p, double T, const Material &t, const  Config& c);

    cdef Result calculate(Projectile &p, const Material &t, const Config &c);
    cdef MultiResult calculate(Projectile &p, const Layers &layers, const Config &c);

    cdef pair[double,double] w_magnification(Projectile p, double E, const Material &t, const Config &c);

cdef extern from "catima/calculations.h" namespace "catima":
    cdef double bethek_lindhard(const Projectile &p);
    cdef double bethek_lindhard_X(const Projectile &p);
    cdef double bethek_dedx_e(Projectile &p,const Target &t, const Config &c, double I);
    cdef double sezi_dedx_e(const Projectile &p, const Target &t);
    cdef double z_effective(const Projectile &p, const  Target &t, const  Config &c);
    cdef double z_eff_Pierce_Blann(double z, double beta);
    cdef double z_eff_Anthony_Landford(double pz, double beta, double tz);
    cdef double z_eff_Hubert(double pz, double E, double tz);
    cdef double z_eff_Winger(double pz, double beta, double tz);
    cdef double z_eff_global(double pz, double E, double tz);
    cdef double z_eff_atima14(double pz, double E, double tz);
    cdef double z_eff_Schiwietz(double pz, double beta, double tz);
    cdef double gamma_from_T(double T);
    cdef double beta_from_T(double T);

cdef extern from "catima/constants.h" namespace "catima":        
    int max_datapoints "catima::max_datapoints"
    int max_storage_data "catima::max_storage_data"
    int logEmin "catima::logEmin"
    int logEmax "catima::logEmax"
    bool reactions "catima::reactions"

cdef extern from "catima/storage.h" namespace "catima":
    cdef cppclass Interpolator:
        Interpolator(const double *x, const double *y, int num) except +
        double eval(double)
        double derivative(double)
        
    cdef cppclass DataPoint:
        Projectile p
        Material m
        Config config
        vector[double] range
        vector[double] range_straggling
        vector[double] angular_variance
    
    cdef cppclass Data:
        Data() except +
        DataPoint& Get(unsigned int i)
        int GetN()
    
    
    cdef cppclass EnergyTableType "catima::EnergyTable[max_datapoints]":
        size_t num;
        double operator()(int i)
        
    cdef EnergyTableType energy_table;
    cdef Data _storage;
    cdef DataPoint& get_data(const Projectile &p, const Material &t, const Config c);
