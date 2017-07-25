/// \file config.h
#ifndef CONFIG
#define CONFIG

#include <cstring>
namespace catima{
    
    /** 
      * \enum z_eff_type 
      * enum to select formulat to calculate effective charge of the Projectile
      */
    enum z_eff_type:char {
        none = 0,
        atima = 1,       // the same as Pierce Blann
        pierce_blann = 1,
        anthony_landorf = 2,
        hubert = 3
    };

    /**
      * enum to select which calculation to skip
      */
    enum skip_calculation:char{
        skip_none = 0,
        skip_tof = 1,
        skip_sigma_a = 2,
        skip_sigma_r = 4
    };

    /**
      * enum to select which dEdx correction to skip
      */
    enum corrections:char{
        no_barkas = 1,
        no_lindhard = 2,
        no_shell_correction = 4
    };

    /**
      * structure to store calculation configuration
      * each group of options are grouped and enum are suppose to use
      * see catima::z_eff_type, catima::skip_calculation, catima::corrections
      *
      * check catima::z_effective()
      * 
      */
    struct Config{
        char z_effective=z_eff_type::atima;
        char skip=skip_none;
        char dedx = 0;
    };
    
    extern Config default_config;
}

#endif
