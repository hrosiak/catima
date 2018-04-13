Catima Caluclator
================
Catima Caluclator is a command line application and interface to the catime library.

Usage
-----
The application is executed from the command line:

```
    catima_calculator config_file";
```

example
```
    catima_calculator c.json
```

Config File Format
------------------
The file must be a valid JSON formatted file.  

The json file should contain the following keys: "projectile", "material", "energy"

#### projectile
The __projectile__ keywords are:
  * array - 2 number array, 1st is mass number, 2nd is charge of the projectile

examples:
```
    "projectile":[11.997,6],
```

#### material
The __material__ keyword is array of object for multi layer material,
or single object defining the material.
The material object must contain __Z__ keyword defining proton number 
of the projectile or the compound material id. 
Optional material object keywords are:
  * __Z__ - proton number or compunds id, mandatory
  * __A__ - mass number of the material, if 0 or undefined elemental atomic weight is used
  * __density__ - density in g/cm3, if 0 or undefined the tabulated density will be used.
  * __thickness__ - material or layer thickness in g/cm2

#### energy
The __energy__ keyword can be 
1.a number specifying the kinetic energy:
```
    "energy":"500.0"
```

2. array of numbers for multiple energies:
```
    "energy":[100,200,500,1000]
```

3. Object specifying minimum energy, maximum energy and energy step, to calculate multiple energies:
```
    "energy":{
        "min": 100,
        "max": 1000,
        "step": 10
    }
```
instead of "step" key the "num" can be specified for integer number of steps between min and max energy.

#### config
The calculation configuration can be change using __config__ keyword. If 
not specified default will be used. The config keyword is expected to be one of the strings
  * "atimav1.3" - for Atima v1.3 setting
  * "atimav1.4" - for Atima v1.4 setting
```
"config":"atimav1.4"
```

Example Files
-------------------
```
{
"projectile":[11.99671, 6],
"energy": 1000,
"material":[{"A":12.0107,
            "Z":6,
            "thickness":1.0
            },
            {
            "A":55.845,
            "Z":26,
            "density":7.8,
            "thickness":0.05
            }
            ],
"config":"atimav1.4"
}
```

```
{
"projectile":[11.99671, 6],
"energy":{
    "min": 100,
    "max": 1000,
    "step": 100
    },
"material":[{"A":12.0107,
            "Z":6,
            "thickness":1.0
            },
            {
            "A":55.845,
            "Z":26,
            "density":7.8,
            "thickness":0.05
            }
            ]
}
```

```
{
"projectile":[11.99671, 6],
"energy":{
    "min": 100,
    "max": 1000,
    "step": 100
    },
"material":{"A":12,
            "Z":6,
            "density":2.0,
            "thickness":1.0
            }
}
```

#### Compound material
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
