"""
    catima python module
    ~~~~~~~~~~~
    This module provides interface to the catima c++ library
    :copyright: (c) 2017 by Andrej Prochazka
    :licence: GNU Affero General Public License, see LICENCE for more details
"""

cimport catimac
from enum import IntEnum
import numpy

cdef class Material:
    cdef catimac.Material cbase
    def __cinit__(self, elements=None, thickness=None, density=None, i_potential=None):
        self.cbase = catimac.Material()
        if(elements and (isinstance(elements[0],float) or isinstance(elements[0],int))):
            self.cbase.add_element(elements[0],elements[1],elements[2])
        if(elements and isinstance(elements[0],list)):
            for e in elements:
                self.cbase.add_element(e[0],e[1],e[2])
        self.cbase.calculate()
        if(not thickness is None):
            self.thickness(thickness)
        if(not density is None):
            self.density(density)
        if(not i_potential is None):
            self.I(i_potential)

    cdef from_c(self, catimac.Material &other):
        self.cbase = other
    
    cdef catimac.Material getc(self):
        cdef catimac.Material res
        res = self.cbase 
        return res
    
    def copy(self):
        res = Material()
        res.cbase = self.cbase
        return res

    def add_element(self, a, z , s):
        self.cbase.add_element(a, z, s)

    def ncomponents(self):
        return self.cbase.ncomponents()

    def molar_mass(self):
        return self.cbase.M()

    def M(self):
        return self.cbase.M()

    def density(self, val=None):
        if(val is None):
            return self.cbase.density()
        else:
            return self.cbase.density(val)

    def thickness(self, val=None):
        if(val is None):
            return self.cbase.thickness()
        else:
            return self.cbase.thickness(val)
    def I(self, val=None):
        if(val is None):
            return self.cbase.I()
        else:
            return self.cbase.I(val)

class material(IntEnum):
            PLASTIC = 201
            AIR = 202
            CH2 = 203
            LH2 = 204
            LD2 = 205
            WATER = 206
            DIAMOND = 207
            GLASS = 208
            ALMG3 = 209
            ARCO2_30 = 210
            CF4 = 211
            ISOBUTANE = 212
            KAPTON = 213
            MYLAR = 214
            NAF = 215
            P10 = 216
            POLYOLEFIN = 217
            CMO2 = 218
            SUPRASIL = 219
            HAVAR = 220
            STEEL = 221
            METHANE = 222

def get_material(int matid):
    res = Material()
    cdef catimac.Material cres = catimac.get_material(matid);
    res.from_c(cres)
    return res


cdef class Target:
    cdef catimac.Target cbase

    def __cinit__(self,a,z,stn):
        self.cbase.A = a
        self.cbase.Z = z
        self.cbase.stn = stn
    
    def A(self):
        return self.cbase.A
    def Z(self):
        return self.cbase.Z
    def stn(self):
        return self.cbase.stn

cdef class Layers:
    cdef public:
        materials
    def __init__(self):
        self.materials=[]
        
    def add(self,Material m):
        self.materials.append(m.copy())
        
    def num(self):
        return len(self.materials)

    def get(self, key):
        return self.materials[key]
        
    def __getitem__(self, key):
        if(isinstance(key,int) and key<self.num()):
            return self.get(key)
        return None
    
    def __add__(self, other):
        res = Layers()
        for e in self.materials:
            res.add(e)
            
        if(isinstance(other,Layers)):
            for e in other.materials:
                res.add(e)
        if(isinstance(other,Material)):
            res.add(other.copy())
        return res
    
    cdef catimac.Layers getc(self):
        cdef catimac.Layers res
        #for l in self.materials:
        #    res.add(l.getc())
        return res

cdef class Projectile:
    cdef catimac.Projectile cbase
    def __cinit__(self, A, Z, Q=None,T=None):
        self.cbase.A = A
        self.cbase.Z = Z
        self.cbase.Q = Z
        if(Q):
            self.cbase.Q = Q
        if(T):
            self.cbase.T = T
    def T(self,val=None):
        if(val is None):
            return self.cbase.T
        self.cbase.T = val;
    def __call__(self,val=None):
        if(val is None):
            return self.cbase.T
        else:
            self.T(val)
            return self
    def A(self):
        return self.cbase.A
    def Z(self):
        return self.cbase.Z
    def Q(self):
        return self.cbase.Q

cdef class Result:
    cdef public double Ein
    cdef public double Eout
    cdef public double Eloss
    cdef public double range
    cdef public double dEdxi
    cdef public double dEdxo
    cdef public double sigma_E
    cdef public double sigma_a
    cdef public double sigma_r
    cdef public double tof
    def __init__(self):
        self.Ein=0.0
        self.Eout=0.0
        self.Eloss=0.0
        self.range=0.0
        self.dEdxi=0.0
        self.dEdxo=0.0
        self.sigma_E=0.0
        self.sigma_a=0.0
        self.sigma_r=0.0
        self.tof=0.0

    def get_dict(self):
        return {"Ein":self.Ein,
                "Eout":self.Eout,
                "Eloss":self.Eloss,
                "range":self.range,
                "dEdxi":self.dEdxi,
                "dEdxo":self.dEdxo,
                "sigma_E":self.sigma_E,
                "sigma_a":self.sigma_a,
                "sigma_r":self.sigma_r,
                "tof":self.tof,
                }

    def __getitem__(self,key):
        d = self.get_dict()
        if(key in d):
            return d[key]

    cdef setc(self,catimac.Result &val):
        self.Ein=val.Ein
        self.Eout=val.Eout
        self.Eloss=val.Eloss
        self.range=val.range
        self.dEdxi=val.dEdxi
        self.dEdxo=val.dEdxo
        self.sigma_E=val.sigma_E
        self.sigma_a=val.sigma_a
        self.sigma_r=val.sigma_r
        self.tof=val.tof

cdef class MultiResult:
    cdef public Result total_result
    cdef public results
    cdef public total

    def __init__(self):
        self.total_result = Result()
        self.results = []
        self.total = {}

    cdef setc(self, catimac.MultiResult &val):
        self.total_result.setc(val.total_result)
        for e in val.results:
            self.results.append(e)
        self.total = self.total_result.get_dict()

    def __getitem__(self,key):
        if(isinstance(key,int) and key<len(self.results)):
            return self.results[key]
        if(isinstance(key,str) and key in self.total):
            return self.total[key]
        return None
    
    def getJSON(self):
        res = {}
        res["result"] = self.total
        res["partial"] = []
        for r in self.results:
            res["partial"].append(r)
        return res

class z_eff_type(IntEnum):
    none = 0,
    pierce_blann = 1
    anthony_landorf = 2
    hubert = 3
    winger = 4
    schiwietz = 5
    global_code = 6
    atima14 = 7

class omega_type(IntEnum):
    atima = 0,
    bohr = 1

class skip_calculation(IntEnum):
    skip_none = 0
    skip_tof = 1
    skip_sigma_a = 2
    skip_sigma_r = 4

class corrections(IntEnum):
        no_barkas = 1
        no_lindhard = 2
        no_shell_correction = 4

cdef class Config:
    cdef catimac.Config cbase
    def __cinit__(self):
        #self.cbase = catimac.Config()
        self.cbase.z_effective = z_eff_type.pierce_blann
        self.cbase.skip = 0
        self.cbase.dedx = 0
        self.cbase.dedx_straggling = omega_type.atima
    def z_effective(self, val=None):
        if(val is None):
            return self.cbase.z_effective
        else:
            self.cbase.z_effective = val
    def skip_calculation(self, val=None):
        if(val is None):
            return self.cbase.skip
        else:
            self.cbase.skip = val
    def dedx(self, val=None):
        if(val is None):
            return self.cbase.dedx
        else:
            self.cbase.dedx = val
    def dedx_straggling(self, val=None):
        if(val is None):
            return self.cbase.dedx_straggling
        else:
            self.cbase.dedx_straggling = val
    def set(self,other):
        if("z_effective" in other):
            self.cbase.z_effective = other["z_effective"]
        if("dedx" in other):
            self.cbase.dedx = other["dedx"]
        if("dedx_straggling" in other):
            self.cbase.dedx_straggling = other["dedx_straggling"]
	
    def get(self):
        res = {}
        res["z_effective"] = self.cbase.z_effective
        res["dedx"] = self.cbase.dedx
        res["dedx_straggling"] = self.cbase.dedx_straggling
        res["skip"] = self.cbase.skip
        return res

    def print_info(self):
        print("z_effective = %s"%z_eff_type(self.cbase.z_effective))
        print("dedx_straggling = %s"%omega_type(self.cbase.dedx_straggling))

default_config = Config()

def calculate(Projectile projectile, material, energy = None, config=default_config):
    if(not energy is None):
        projectile.T(energy)
    if(isinstance(material,Material)):
        return calculate_material(projectile, material, config = config)
    if(isinstance(material,Layers)):
        return calculate_layers(projectile, material, config = config) 

def calculate_material(Projectile projectile, Material material, energy = None, Config config = default_config):
    if(not energy is None):
        projectile.T(energy)
    cdef catimac.Result cres = catimac.calculate(projectile.cbase,material.cbase,config.cbase)
    res = Result()
    res.setc(cres)
    return res

def calculate_layers(Projectile projectile, Layers layers, energy = None,  Config config = default_config):
    cdef catimac.Layers clayers
    clayers = catimac.Layers()
    clayers = get_clayers(layers)
    if(not energy is None):
        projectile.T(energy)
    cdef catimac.MultiResult cres = catimac.calculate(projectile.cbase, clayers, config.cbase)
    res = MultiResult()
    res.setc(cres)
    return res

cdef catimac.Layers get_clayers(Layers layers):
    cdef catimac.Layers res
    cdef catimac.Material m
    for l in layers.materials:
        m = get_cmaterial(l)
        res.add(m)
    return res

cdef catimac.Material get_cmaterial(Material material):
    cdef catimac.Material res
    res = material.cbase
    return res

def projectile_range(Projectile projectile, Material material, energy = None, Config config = default_config):  
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.range(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()    
    return catimac.range(projectile.cbase, energy, material.cbase, config.cbase);

def dedx_from_range(Projectile projectile, Material material, energy = None, Config config = default_config):  
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.dedx_from_range(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()    
    return catimac.dedx_from_range(projectile.cbase, energy, material.cbase, config.cbase);

def domega2de(Projectile projectile, Material material, energy = None, Config config = default_config):  
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.domega2de(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()    
    return catimac.domega2de(projectile.cbase, energy, material.cbase, config.cbase);

def da2de(Projectile projectile, Material material, energy = None, Config config = default_config):  
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.da2de(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()    
    return catimac.da2de(projectile.cbase, energy, material.cbase, config.cbase);

def dedx(Projectile projectile, Material material, energy = None, Config config = default_config):
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.dedx(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()
    return catimac.dedx(projectile.cbase, energy, material.cbase, config.cbase)
    
def domega2dx(Projectile projectile, Material material, energy = None, Config config = default_config):
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.domega2dx(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()
    return catimac.domega2dx(projectile.cbase, energy, material.cbase, config.cbase)

def energy_out(Projectile projectile, Material material, energy = None, Config config = default_config):
    if(isinstance(energy,numpy.ndarray)):
        res = numpy.empty(energy.size)
        for i,e in enumerate(energy):
            res[i] = catimac.energy_out(projectile.cbase, e, material.cbase, config.cbase)
        return res    
    if(energy is None):
        energy = projectile.T()
    return catimac.energy_out(projectile.cbase, energy, material.cbase, config.cbase)

def lindhard(Projectile projectile):
    return catimac.bethek_lindhard(projectile.cbase);

def lindhard_X(Projectile projectile):
    return catimac.bethek_lindhard_X(projectile.cbase);

def z_effective(Projectile p, Target t, Config c = default_config):
    return catimac.z_effective(p.cbase, t.cbase, c.cbase)

def z_eff_Pierce_Blann(double z, double beta):
    return catimac.z_eff_Pierce_Blann(z,beta)
    
def z_eff_Anthony_Landford(double pz, double beta, double tz):
    return catimac.z_eff_Anthony_Landford(pz, beta, tz);

def z_eff_Hubert(double pz, double E, double tz):
    return catimac.z_eff_Hubert(pz, E, tz);

def z_eff_Winger(double pz, double beta, double tz):
    return catimac.z_eff_Winger(pz, beta, tz);

def z_eff_global(double pz, double E, double tz):
    return catimac.z_eff_global(pz, E, tz);

def z_eff_atima14(double pz, double E, double tz):
    return catimac.z_eff_atima14(pz, E, tz);

def z_eff_Schiwietz(double pz, double beta, double tz):
    return catimac.z_eff_Schiwietz(pz, beta, tz);

def gamma_from_T(double T):
    return catimac.gamma_from_T(T);

def beta_from_T(double T):
    return catimac.beta_from_T(T);

def get_data(Projectile projectile, Material material, Config config = default_config):
    data = catimac.get_data(projectile.cbase, material.cbase, config.cbase)
    return [data.range,data.range_straggling,data.angular_variance]

# constants
max_datapoints = catimac.max_datapoints
max_storage_data = catimac.max_storage_data
logEmin = catimac.logEmin
logEmax = catimac.logEmax

def energy_table(unsigned int i):
    if(i<catimac.energy_table.num):
        return catimac.energy_table(i)
    else:
        return -1.0

def get_energy_table():
    r = [catimac.energy_table(x) for x in range(catimac.energy_table.num)]
    return r

def storage_info():
    res = []
    for i in range(catimac.max_storage_data):
        data = catimac._storage.Get(i)
        if(data.p.A>0 and data.p.Z>0 and data.m.ncomponents()>0):
            matter = []
            for j in range(data.m.ncomponents()):
                e = data.m.get_element(j)
                matter.append([e.A,e.Z,e.stn])
            res.append({"projectile":[data.p.A,data.p.Z],"matter":matter, "config":data.config})
    return res

def catima_info():
    print("CATIMA version = 1.1")
    print("number of energy points = %g"%max_datapoints)
    print("min energy point = 10^%g MeV/u"%logEmin)
    print("max energy point = 10^%g MeV/u"%logEmax)
