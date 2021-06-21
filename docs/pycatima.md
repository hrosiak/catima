Usage
=====
Installation
------------
Easiest way is to install pycatima on Linux and Windows is using pip:
```
pip install pycatima
```
note: python 3.7-3.9 is required


Pojectile
---------

Projectile is defined by __pycatima.Projectile__ class. It is initialized using:

**Projectile(A, Z, Q=Z, T=0)**

where __A__ is mass in u units, __Z__ is proton number, __Q__ is charge state 
and __T__ is energy in Mev/u units.

```python
import pycatima
p = pycatima.Projectile(238.00032 ,92)  # 238U92+ at 0 Mev/u 
p = pycatima.Projectile(238.00702 ,92,90,1000)  # 238U90+ at 1000 Mev/u 

# following methods are defined:
mass = p.A()  # get mass of the nucleus us u units
z = p.Z()  # get proton number
q = p.Q()  # get charge state
energy = p.T()  # get energy
p.T(1000)  # set energy
```

Material
--------
Material is defined by __pycatima.Material__ class. The recommended way of 
initialization is usign the following init signature:

**Material(elements, density, thickness, i_potential, mass)**
  * elements - list of elements, where element is defined as list of [A, Z, STN]
  A is atomic mass of the element, if 0 natural abundance atomic mass is taken, 
  Z is the proton number of the element, STN is the stoichiometric if >=1.0 or
  weight fraction if < 1.0
  * density - optional, defaults to 0
  * thickness - optional, if not defined 0
  * i_potential - optional, if <=0 it will be calulated using Bragg rule from elemental ionization potentials
  * mass - optional, if <=0 it will be calculated from elements masses and STN number.

__Material__ class has following methods:

  * **add_element(a, z, stn)** - adds another element to the material, see initialization comments for details about, a,z, stn 
  * **ncomponents()** - returns number of elements in the material
  * **density()** - returns density 
  * **density(value)** - set density
  * **thickness()** - get thickness in g/cm^2
  * **thickness(value)** - set thickness in g/cm^2
  * **thickness_cm(value)** - set thickness in cm
  * **I()** - get mean ionization potential
  * **I(value)** - set mean ionization potential

### Default Materials
The material with predefined density and atomic weights can be obtained for elemental targets and some compounds using:

**pycatima.get_material(id)** 

id is integer which identifies material, For elements use 1-99 for elemental targets and >200 for compounds. See available pre-defined compounds in __predefined material__ section in manual./


### Example 

```python
import pycatima

#H2O with natural atomic masses and Ipot=78eV
h2o = pycatima.Material([[0,1,2],[0,8,1]],density=1.0, i_potential=78)
h2o.thickness_cm(1.0) # 1cm of water
# C target
c_mat = pycatima.get_material(6)
c_mat.thickness(0.1) # set to 0.1g/cm2 
```  

Layers
------
Layers are sequential layers of __Material__ which __Projectile__ passes.
The __pycatima.Layers__ class is defined using following signature:
**Layers()**
This creates empty layers which needs to be filled by __Material__ class.

__Layers__ clas has following methods:
  * **add(material)** - add material to the layers, material must be __Material__ class
  * **add_layers(other)** - add all material from __Layers__ class __other__
  * **num()** - returns number of layers
  * **__get_item__(i)** returns i-th Material class
  * **__get(i)** returns i-th Material class
  
### Example

```python
import pycatima
layers = pycatima.Layers()

# define some materials
graphite = pycatima.get_material(6)
graphite.thickness(0.2)
p10 = pycatima.get_material(pycatima.material.P10)
p10.thickness_cm(2.0)
air = pycatima.get_material(pycatima.material.Air)
air.thickness_cm(2.0)

# now add materials to layers
layers.add(graphite)
layers.add(air)
graphite.thickness(0.1) # change thickness for next layer
layers.add(graphite)
layers.add(p10)
```



Calculation
-----------
Calculations are done using mainly via functions:

* **calculate(Projectile, Material)**
* **calculate_layers(Projectile, Layers)**

The __calculate__ function returns __Result__ class and __calculate_layers__ returns __MultiResult_class__

### Example
```python
import pycatima
water = catima.get_material(catima.material.Water)
water.thickness(1.0)
p = catima.Projectile(1,1)
p.T(1000)  # set projectile energy to 1000MeV/u

res = catima.calculate(p,water) # now res contains results
d = res.get_dict() # get results as dictionary
```

```python
p = catima.Projectile(12,6)
water = catima.get_material(catima.material.Water)
water.thickness(10.0)
graphite = catima.get_material(6)
graphite.thickness(1.0)
graphite.density(2.0)

mat = catima.Layers()
mat.add(water)
mat.add(graphite)

# now calculate results for projectile at 1000MeV/u
res = catima.calculate_layers(p(1000),mat)
```

Results
-------
__Results__ class stores results for 1 layer of __Material__. It has following variables:

  * Ein - Energy at entrance of the material in MeV/u
  * Eloss - Energy loss in material in MeV
  * Eout - Energy at the end of material in MeV/i
  * dEdxi - Stopping power at entrance
  * dEdxo - Stopping power at exit
  * range - range in the material in g/cm^2
  * sigma_E - Energy straggling in MeV/u
  * sigma_a - Angular straggling in rad
  * sigma_r - range straggling in g/cm^2
  * sigma_x - position straggling in cm
  * sp - non-reaction probability
  * tof - time of flight through the material

__MultiResults__ class stores the results for multiple layers. It consists of following variables:
 
  * total_result - __Result__ class with total results for projectile passing all layers.
  * results - list of __Result__ classes, one for each layer in __Layers__



Config
------
TODO
