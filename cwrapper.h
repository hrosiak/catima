#ifndef CATIMA_CWRAPPER
#define CATIMA_CWRAPPER

#ifdef __cplusplus
extern "C" {
#endif

struct CatimaResult{
        double Ein;
        double Eout;
        double Eloss;
        double range;
        double dEdxi;
        double dEdxo;
        double sigma_E;
        double sigma_a;
        double sigma_r;
        double tof;   
};


    enum z_eff_type {
        none = 0,
        atima = 1
    };

    enum skip_calculation{
        skip_none = 0,
        skip_tof = 1,
        skip_sigma_a = 2,
        skip_sigma_r = 4
    };

struct CatimaConfig {
        char z_effective;
        char skip;
};

struct CatimaConfig catima_defaults = {none,skip_none};

typedef struct CatimaResult CatimaResult;

CatimaResult catima_calculate(double pa, int pz, double T, double ta, double tz, double thickness, double density);
double catima_angular_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz);
double catima_energy_straggling_from_E(double pa, int pz, double Tin, double Tout,double ta, double tz);


#ifdef __cplusplus
}
#endif

#endif