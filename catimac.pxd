from libcpp.pair cimport pair
from libcpp.vector cimport vector

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

    cdef cppclass MultiResult

    cdef cppclass Material:
        Material() except +
        void add_element(double , int , double )
        pair[Target,double] getElement(int) 
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
    cdef int max_datapoints;
    
cdef extern from "catima/storage.h" namespace "catima":
    cdef cppclass Interpolator:
        double eval(double)
        double derivative(double)
