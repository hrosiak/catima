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
      * enum to select which dEdx correction to skip
      */
    enum corrections:unsigned char{
        no_barkas = 1,
        no_lindhard = 2,
        no_shell_correction = 4,
        no_highenergy = 8
    };

    /**
      * enum to select dEdx straggling options
      */
    enum omega_types:unsigned char{
        atima = 0,
        bohr = 1,
    };

    /**
      * enum to select how low energy part is calculated
      */
    enum low_energy_types:unsigned char{
        srim_85 = 0,
        srim_95 = 1,
    };

    /**
      * enum to select angular scattering
      */
    enum scattering_types:unsigned char{
        fermi_rossi = 0,
        atima_scattering = 255,
    };

    /**
      * structure to store calculation configuration
      */
    struct Config{
        #ifndef GLOBAL
        unsigned char z_effective=z_eff_type::pierce_blann;
        #else
        unsigned char z_effective=z_eff_type::atima14;
        #endif

        unsigned char corrections = 0;
        unsigned char calculation = 1;
        unsigned char low_energy = low_energy_types::srim_85;
        unsigned char scattering = scattering_types::atima_scattering;
    };


    extern Config default_config;
}

#endif
