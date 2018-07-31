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
        pierce_blann = 1,
        anthony_landorf = 2,
        hubert = 3,
        winger = 4,
        schiwietz = 5,
        global = 6,
        atima14 = 7
    };

    /**
      * enum to select which calculation to skip
      */
    enum skip_calculation:unsigned char{
        skip_none = 0,
        skip_tof = 1,
        skip_sigma_a = 2,
        skip_sigma_r = 4,
        skip_reactions = 128
    };

    /**
      * enum to select which dEdx correction to skip
      */
    enum corrections:unsigned char{
        no_barkas = 1,
        no_lindhard = 2,
        no_shell_correction = 4
    };

    /**
      * enum to select which dEdx straggling options
      */
    enum omega:unsigned char{
        atima = 0,
        bohr = 1,
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
        unsigned char z_effective=z_eff_type::pierce_blann;
        //char z_effective=z_eff_type::atima14;
        unsigned char dedx = 0;
        unsigned char dedx_straggling = omega::atima;
        #ifdef NUREX
        unsigned char skip=skip_none;
        #else
        unsigned char skip=skip_calculation::skip_reactions;
        #endif
    };
    
    extern Config default_config;
}

#endif
