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
The __energy__ keyword can be:

  * a number specifying the kinetic energy:
```
    "energy":"500.0"
```

  * array of numbers for multiple energies:
```
    "energy":[100,200,500,1000]
```
  * Object specifying minimum energy, maximum energy and energy step, to calculate multiple energies:
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
