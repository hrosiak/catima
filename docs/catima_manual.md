CATima library manual
=====================

Units
------
The following units are used for input and outputs:
  * projectile energy - MeV/u
  * density - g/cm^3
  * material thickness - g/cm^2
  * material length - cm
  * angle - rad
  * time of flight - ns


Projectile
----------
The __Projectile__ class is used to store projectile data.
Each projectile must provide A,Z combination, additionally charge state can be set as well.
The example of projectile definition:
```cpp
catima::Projectile p1(12,6); //12C projectile
catima::Projectile p2(12,6,5); //12C(5+) projectile
```

to set the energy of the projectile in MeV/u units:
```cpp
p1.T = 1000.0;
p2(1000.0);
```

Material Definition
-------------------
The material is defined and stored in __Material__ class.
The class stores the atom constituents, thickness and density of the material etc.

There are 2 ways to define materials: specifying all elements in constructor or using `add_element(double, int, double)` function.
The example of water definition:

```cpp

catima::Material carbon({0,6,1}); // carbon with elemental atomic weight

catima::Material water1({
	{1,1,2}, // {weight, Z, stn or weight fraction}
	{16,8,1}
	});
water1.density(1.0);
water1.thickness(2.0);
water1.I(78.); // set custom ionization potential in eV

catima::Material water2;
water2.add_element(1,1,2);
water2.add_element(16,8,1);
water2.density(1.0).thickness(2.0);
```
If mass number is equal to 0, the mass number of the element is taken as atomic weight of the element.
Compound elements can be defined either via stoichiometric number or via weight fraction. If the number is less than 1
it is assumed weight fraction is being used, otherwise stoichiometric number or molar fraction is being used.
```cpp
catima::Material air ({{0,7,0.755267},{0,8,0.231781},{0,18,0.012827},{0,6,0.000124}},0.001205); // weight fractions
catima::Material water ({{0,1,2},{0,8,1}},1); // mole fraction
```
### predefined materials ###
If the library is compiled with predefined materials database, the Material can be retrieved from the database as:
```cpp
using namespace catimal
Material water = get_material(material::WATER);
Material graphite = get_material(6);
```

The list of predefined material can be found at __material_database.h__ file


Calculation
-----------
To calculate all observable following function can be used:
```cpp
Result calculate(Projectile &p, const Material &t, double T, Config c)
Result calculate(Projectile &p, const Material &t, Config c)
```

Both function returns structure ___Result___ which contains all calculated variables.

```cpp
    struct Result{
        double Ein=0.0;
        double Eout=0.0;
        double Eloss = 0.0;
        double range=0.0;
        double dEdxi=0.0;
        double dEdxo=0.0;
        double sigma_E=0.0;
        double sigma_a=0.0;
        double sigma_r=0.0;
        double tof=0.0;
    };
```


If one is interested only in one of the variable, the following function can be used:

```cpp
double dedx(Projectile &p, double T, const Material &t, Config c=default_config);
double domega2dx(Projectile &p, double T, const Material &t, Config c=default_config);
double range(Projectile &p, double T, const Material &t, Config c=default_config);
double range_straggling(Projectile &p, double T, const Material &t, Config c=default_config);
double angular_straggling(Projectile &p, double T, const Material &t, Config c=default_config);
```

Example calculation:
```cpp
...
double T=1000;
auto result = catima::calculate(p1(1000),water1);
cout<<"T "<<T<<", dEdx = "<<result.dEdxi<<" MeV/g/cm2"<<", range = "<<result.range<<" g/cm2"<<endl;
```

Multilayer Material
-------------------
The layers of __Materials__ are stored in __catima::Layers__ class.

There are following ways to define Layers from catima::Material classes:


```cpp
catima::Material graphite({12,6,1});
catima::Material nitrogen({14,7,1});
...
catima::Layers matter1;
matter1.add(graphite);
matter1.add(nitrogen);
matter1.add(graphite);
cout<<"number of layers = "<<matter1.num()<<"\n"; // 3
```
Layers can be copied from existing Layes:
```cpp
catima::Layers matter2;
matter2 = matter1; //matter2 contain 3 layers
matter2.add(nitrogen);
matter2.add(graphite); //matter2 contains 5 layers now
```

Layers can be created as a combination of another Layers:
```cpp
catima::Layers matter3 = matter1 + matter2;
```

Config
------
The calculation configuration is set via catima::Config class. The dafault configuration is predefined in __catima::default_config__ variable.
This __default_config__ is supplied as default argument to functions like catima::calculate. If custom config is needed another configuration can be provided.

the structure Config is defined as:
```cpp
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
```

### effective charge calculation###
the following effective charge calculations are buit in:
```
    enum z_eff_type:char {
        none = 0,
        pierce_blann = 1,
        anthony_landorf = 2,
        hubert = 3,
        winger = 4,
        schiwietz = 5,
        global = 6,
        atima14 = 7
    }; hubert = 3
```
  * z_eff_type::none - the provided Projectile Q is used as a charge
  * z_eff_type::pierce_blann - Pierce Blann formula, using function: z_eff_Pierce_Blann()
  * z_eff_type::anthony_landorf - function: z_eff_Anthony_Landorf()
  * z_eff_type::hubert - function: z_eff_Hubert()
  * z_eff_type::winger - function: z_eff_Winger()
  * z_eff_type::global - function: z_eff_global()
  * z_eff_type::atima14 - function: z_eff_atima14()




All available switches are defined in __config.h__ file.



Using the library
=================
the include direcotry and LD_LIBRARY_PATH must be properly adjusted.
The app must be linked against catima library.
For example check examples directory and makefile inside to see how to link.

All functions and classes are inside __catima namespace__

Normally including main file is enough:
```cpp
#include "catima/catima.h"
```


Using with C
-------------
the C wrapper is provided in cwapper.h, this file can be included in C app. The C app must be then linked against catima library.
It provides only basic interface.



Compound Pre-defined Materials
==============================
The predefined compound material can be used using the __Z__ field as an ID of the compound.
The following are supported:
```
		 Plastics = 201,
		 Air = 202,
		 CH2 = 203,
		 LH2 = 204,
		 LD2 = 205,
		 Water = 206,
		 Diamond = 207,
		 Glass = 208,
		 ALMG3 = 209,
		 ArCO2_30 = 210,
		 CF4 = 211,
		 Isobutane = 212,
		 Kapton = 213,
		 Mylar = 214,
		 NaF = 215,
		 P10 = 216,
		 Polyolefin = 217,
		 CmO2 = 218,
		 Suprasil = 219,
		 HAVAR = 220,
		 Steel = 221,
		 CO2 = 222,
		 CH4 = 223,
		 Methanol = 224,
		 Acetone = 225,
		 Acetylene = 226,
		 Adenine = 227,
		 Adipose_Tissue = 228,
		 Alanine = 229,
		 Bakelite = 230,
		 AgBr = 231,
		 AgCl = 232,
		 AgI = 233,
		 Al2O3 = 234,
		 Amber = 235,
		 Ammonia = 236,
		 Aniline = 237,
		 Anthracene = 238,
		 A_150 = 239,
		 B_100 = 240,
		 BaF2 = 241,
		 BaSO4 = 242,
		 Benzene = 243,
		 BeO = 244,
		 BGO = 245,
		 Blood = 246,
		 Bone_Compact = 247,
		 Bone_Cortical = 248,
		 Brain_ICRP = 249,
		 B4C = 250,
		 BC400 = 251,
		 nButanol = 252,
		 C_552 = 253,
		 CdTe = 254,
		 CdWO4 = 255,
		 CaCO3 = 256,
		 CaF2 = 257,
		 CaO = 258,
		 CaWO4 = 259,
		 CsF = 260,
		 CsI = 261,
		 CCl4 = 262,
		 C2Cl4 = 263,
		 Cellophane = 264,
		 Chlorobenzene = 265,
		 Chloroform = 266,
		 Cyclohexane = 267,
		 Concrete = 268,
		 Diethyl_Ether = 269,
		 Ethane = 270,
		 Ethanol = 271,
		 Ethylene = 272,
		 Eye_lens = 273,
		 Fe2O3 = 274,
		 FeO = 275,
		 Freon12 = 276,
		 Freon12B2 = 277,
		 Freon13 = 278,
		 Freon13B1 = 279,
		 Freon13I1 = 280,
		 Gd2O2S = 281,
		 GaAs = 282,
		 Gel_Photo_Emulsion = 283,
		 Glass_Pyrex = 284,
		 Glass_Lead = 285,
		 Glucose = 286,
		 Glutamine = 287,
		 Glycerol = 288,
		 Guanine = 289,
		 Gypsum = 290,
		 nHeptane = 291,
		 nHexane = 292,
		 KI = 293,
		 K2O = 294,
		 LaBr3 = 295,
		 LaOBr = 296,
		 La2O2S = 297,
		 Lung = 298,
		 MgCO3 = 299,
		 MgF2 = 300,
		 MgO = 301,
		 MS20_Tissue = 302,
		 Muscle_skeletal = 303,
		 Muscle_strained = 304,
		 Muscle_sucrose = 305,
		 Muscle_no_sucrose = 306,
		 Na2CO3 = 307,
		 NaI = 308,
		 NaCl = 309,
		 Na2O = 310,
		 NaNO3 = 311,
		 Naphthalene = 312,
		 Nitrobenzene = 313,
		 N2O = 314,
		 Octane = 315,
		 Paraffin = 316,
		 nPentane = 317,
		 PhotoEmulsion = 318,
		 PuO2 = 319,
		 Polyacrylonitrile = 320,
		 Polycarbonate = 321,
		 PMMA = 322,
		 POM = 323,
		 Polypropylene = 324,
		 Polystyrene = 325,
		 Propane = 326,
		 nPropanol = 327,
		 PVC = 328,
		 Pyridine = 329,
		 SiO2 = 330,
		 Skin = 331,
		 Sucrose = 332,
		 Teflon = 333,
		 TlCl = 334,
		 Toluene = 335,
		 Trichloroethylene = 336,
		 WF6 = 337,
		 UC2 = 338,
		 UC = 339,
		 UO2 = 340,
		 Urea = 341,
		 Valine = 342
```
