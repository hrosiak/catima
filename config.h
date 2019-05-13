/// \file config.h
#ifndef CONFIG
#define CONFIG
#include <cstring>
namespace catima{
    
    /** 
      * \enum z_eff_type 
      * enum to select formulat to calculate effective charge of the Projectile
      */
    enum z_eff_type:unsigned char {
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
        no_shell_correction = 4,
        no_highenergy = 8
    };

    /**
      * enum to select which dEdx straggling options
      */
    enum omega_types:unsigned char{
        atima = 0,
        bohr = 1,
    };

    /**
      * enum to select which how low energy part is calculated
      */
    enum low_energy_types:unsigned char{
        srim_85 = 0,
        srim_95 = 1,
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
        #ifndef GLOBAL
        unsigned char z_effective=z_eff_type::pierce_blann;
        #else
        unsigned char z_effective=z_eff_type::atima14;
        #endif

        #ifdef REACTIONS
        unsigned char skip=skip_none;
        #else
        unsigned char skip=skip_calculation::skip_reactions;
        #endif

        unsigned char corrections = 0;
        unsigned char calculation = 1;
    };

    inline void set_config_lowenergy(Config c, low_energy_types lt){
      c.calculation = c.calculation & (lt<<2);
    }

    inline unsigned char config_lowenergy(const Config c){
      return (c.calculation>>2) & 0x7;
    }

    inline void set_config_omega(Config c, omega_types ot){
      c.calculation = c.calculation & ot;
    }

    inline unsigned char config_omega(const Config c){
      return c.calculation & 0x3;
    }

    
    extern Config default_config;
}

#endif
