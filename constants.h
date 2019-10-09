#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <limits>
#include "catima/build_config.h"
namespace catima {

// config
constexpr double Ezero = 1E-3; // lowest E to calculate, below taken as 0
constexpr double logEmin = -3; // log of minimum energy
constexpr double logEmax = 7.0;  // log of max energy
constexpr int max_datapoints = 500; // how many datapoints between logEmin and logEmax
constexpr int max_storage_data = 60; // number of datapoints which can be stored in cache
constexpr double numeric_epsilon = 10*std::numeric_limits<double>::epsilon();

constexpr double thin_target_limit = 1 - 1e-3;

#ifdef REACTIONS
constexpr double emin_reaction = 30.0;
constexpr bool reactions = true;
#else
constexpr bool reactions = false;
#endif

// constants
constexpr double PI = 3.1415926535897932384626433832795;
constexpr double Avogadro = 6.022140857; // 10^23
constexpr double electron_mass = 0.510998928;   // MeV/c^2
constexpr double atomic_mass_unit = 931.4940954; // MeV/c^2
constexpr double classical_electron_radius = 2.8179403227; //fm
constexpr double fine_structure = 1/137.035999139;
constexpr double fine_structure_inverted = 1/fine_structure;
constexpr double c_light = 299.792458; //Mm/s
constexpr double bohr_velocity = 2.19 / c_light; // in c unit
constexpr double hbar = 6.582119514; // in eV*s * 10^-16

constexpr double dedx_constant = 0.3070749187;  //4*pi*Na*me*c^2*r_e^2  //MeV cm^2
constexpr double domega2dx_constant = dedx_constant*electron_mass;  //4*pi*Na*me*c^2*r_e^2  //MeV^2 cm^2


    // units //
    namespace units{
        constexpr double g = 1.0;
        constexpr double mg = 1000.0;
        constexpr double cm3 = 1.0;
        constexpr double cm = 1.0;
        constexpr double mm = 10.;
        constexpr double keV = 1000.0;
        constexpr double ns = 1.0;
    }

}


#endif
