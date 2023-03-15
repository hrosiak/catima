import sys
sys.path.insert(0,"../build")
import unittest
import pycatima as catima
import math

class TestStructures(unittest.TestCase):

    def test_Projectile(self):
        p = catima.Projectile(238,92)
        self.assertEqual(p.A(),238)
        self.assertEqual(p.Z(),92)
        self.assertEqual(p.Q(),92)

        p = catima.Projectile(238,92,90)
        self.assertEqual(p.A(),238)
        self.assertEqual(p.Z(),92)
        self.assertEqual(p.Q(),90)
        p.T(1000)
        self.assertEqual(p.T(),1000)
        p(500)
        self.assertEqual(p.T(),500)

        p = catima.Projectile(238,92,90, T=100)
        self.assertEqual(p.T(),100)
    
    def test_Material(self):
        mat = catima.Material()
        mat.add_element(12,6,1)
        self.assertEqual(mat.ncomponents(),1)
        mat.add_element(1,1,2)
        self.assertEqual(mat.ncomponents(),2)

        mat2 = catima.Material(12.01,6)
        self.assertEqual(mat2.ncomponents(),1)
        self.assertEqual(mat2.molar_mass(),12.01)

        mat3 = catima.Material([[12,6,1]])
        self.assertEqual(mat3.ncomponents(),1)
        self.assertEqual(mat3.molar_mass(),12)

        water = catima.Material([[1,1,2],[16,8,1]])
        self.assertEqual(water.molar_mass(),18)

        mat2 = catima.Material(0,6)
        self.assertEqual(mat2.ncomponents(),1)
        self.assertAlmostEqual(mat2.molar_mass(),12,1)

        mat5 = catima.Material(0,6,density=1.9, thickness=0.5)
        self.assertEqual(mat5.ncomponents(),1)
        self.assertEqual(mat5.thickness(),0.5)
        self.assertEqual(mat5.density(),1.9)

        mat6 = catima.Material(0,6,density=1.9, thickness=0.5,i_potential=80.0)
        self.assertEqual(mat6.ncomponents(),1)
        self.assertEqual(mat6.thickness(),0.5)
        self.assertEqual(mat6.density(),1.9)
        self.assertEqual(mat6.I(),80.0)

        # copy
        mat3.density(1.8)
        matc = catima.Material(mat3)
        self.assertEqual(matc.ncomponents(),1)
        self.assertEqual(matc.molar_mass(),12)
        self.assertEqual(matc.density(),1.8)
        mat3.density(2.0)
        self.assertEqual(matc.density(),1.8)
        self.assertEqual(mat3.density(),2.0)

    
    def test_default_material(self):
        m1 = catima.get_material(6);
        self.assertAlmostEqual(m1.molar_mass(),12,1)
        self.assertEqual(m1.ncomponents(),1)
        self.assertAlmostEqual(m1.density(),2.0,1)

        m2 = catima.get_material(catima.material.Water)
        self.assertEqual(m2.ncomponents(),2)
        self.assertAlmostEqual(m2.molar_mass(),18,1)
        self.assertAlmostEqual(m2.density(),1.0,1)

        m3 = catima.get_material(3001)
        self.assertEqual(m3.ncomponents(),0)
        self.assertAlmostEqual(m3.molar_mass(),0,1)
        self.assertAlmostEqual(m3.density(),0.0,1)

    def test_layers(self):
        graphite = catima.get_material(6)
        graphite.thickness(0.5)
        p10 = catima.get_material(catima.material.P10)
        p10.thickness(0.01)
        n2 = catima.get_material(7)
        n2.thickness(0.02)

        mat= catima.Layers()
        self.assertEqual(mat.num(),0)
        mat.add(graphite)
        self.assertEqual(mat.num(),1)
        self.assertAlmostEqual(mat[0].molar_mass(),12,1)
        self.assertAlmostEqual(mat[0].thickness(),0.5,1)
        self.assertAlmostEqual(mat[0].density(),2.0,1)
        
        mat.add(p10)
        self.assertEqual(mat.num(),2)

        graphite.thickness(1.0)
        graphite.density(1.8)
        mat.add(graphite)
        self.assertEqual(mat.num(),3)
        self.assertAlmostEqual(mat[2].molar_mass(),12,1)
        self.assertAlmostEqual(mat[0].thickness(),0.5,1)
        self.assertAlmostEqual(mat[0].density(),2.0,1)
        self.assertAlmostEqual(mat[2].thickness(),1.0,1)
        self.assertAlmostEqual(mat[2].density(),1.8,1)
        a = mat[2]
        a.thickness(1.2)
        a.density(1.9)
        self.assertAlmostEqual(mat[2].thickness(),1.2,1)
        self.assertAlmostEqual(mat[2].density(),1.9,1)
        
        self.assertEqual(mat.num(),3)
        with self.assertRaises(ValueError):
            mat[3]
    
        mat2 = catima.Layers()
        mat2.add(n2)
        self.assertEqual(mat2.num(),1)

        mats = mat2 + mat
        self.assertEqual(mats.num(),4)
        self.assertAlmostEqual(mats[0].molar_mass(),14,1)
        self.assertEqual(mats[0].thickness(),0.02)
        self.assertAlmostEqual(mats[1].molar_mass(),12,1)
        self.assertAlmostEqual(mats[3].molar_mass(),12,1)

        n2.thickness(0.5)
        mats = mats + n2
        self.assertEqual(mats.num(),5)
        self.assertAlmostEqual(mats[0].molar_mass(),14,1)
        self.assertEqual(mats[0].thickness(),0.02)
        self.assertAlmostEqual(mats[4].molar_mass(),14,1)
        self.assertEqual(mats[4].thickness(),0.5)

    def test_material_calculation(self):
        water = catima.get_material(catima.material.Water)
        p = catima.Projectile(1,1)
        
        p(1000)
        res = catima.calculate(p,water)
        res2 = catima.dedx(p,water)
        self.assertAlmostEqual(res.dEdxi,2.23,1)
        self.assertAlmostEqual(res.dEdxi,res2,3)
        res = catima.calculate(p(500),water)
        res2 = catima.dedx(p,water)
        self.assertAlmostEqual(res.dEdxi,2.76,1)
        self.assertAlmostEqual(res.dEdxi,res2,3)

        res = catima.calculate(p(9),water)
        res2 = catima.dedx(p,water)
        self.assertAlmostEqual(res.dEdxi,51.17,1)
        self.assertAlmostEqual(res.dEdxi,res2,3)
        res = catima.calculate(p(9),water)
        res = catima.calculate(p(9),water)
        self.assertAlmostEqual(res.dEdxi,51.17,1)
        
        p(900000)
        res = catima.calculate(p,water)
        res2 = catima.dedx_from_range(p,water)
        self.assertAlmostEqual(res.dEdxi,res2,3)
        
    def test_config(self):
        water = catima.get_material(catima.material.Water)
        water.density(1.0)
        water.thickness(1.0)
        p = catima.Projectile(1,1)
        conf = catima.Config()
        conf.calculation = catima.omega_types.bohr
        conf2 = catima.Config()
        conf2.calculation = catima.omega_types.atima
        self.assertEqual(conf.calculation, 1)
        self.assertEqual(conf2.calculation, 0)
        p(1000)
        res = catima.calculate(p,water,config=conf)
        res2 = catima.calculate(p,water,config=conf2)
        self.assertAlmostEqual(res.dEdxi,res2.dEdxi,delta=1e-6)
        self.assertNotAlmostEqual(res.sigma_E,res2.sigma_E,delta=1e-4)
        self.assertNotAlmostEqual(res.sigma_r,res2.sigma_r,delta=1e-4)
        
    
    def test_eout(self):
        graphite = catima.get_material(6)
        graphite.thickness(0.5)
        p = catima.Projectile(12,6)
        res = catima.calculate(p(1000),graphite)
        res2 = catima.energy_out(p(1000),graphite)
        self.assertAlmostEqual(res.Eout,997.077,1)
        self.assertAlmostEqual(res.Eout,res2,3)
    
    def test_eout_list(self):
        graphite = catima.get_material(6)
        graphite.thickness(0.5)
        p = catima.Projectile(12,6)
        energies = [100,500,1000]
        res = catima.calculate(p(1000),graphite)
        self.assertAlmostEqual(res.Eout,997.077,1)
        res2 = catima.energy_out(p,energies, graphite)
        self.assertEqual(len(res2),len(energies))
        self.assertAlmostEqual(res2[2], 997.077,1)
        self.assertAlmostEqual(res2[0], catima.calculate(p(energies[0]),graphite).Eout ,1)
        self.assertAlmostEqual(res2[1], catima.calculate(p(energies[1]),graphite).Eout ,1)
    
    def test_dedx_from_range_list(self):
        graphite = catima.get_material(6)
        graphite.thickness(0.5)
        p = catima.Projectile(12,6)
        energies = [100,500,1000]
        res2 = catima.dedx_from_range(p,energies, graphite)
        self.assertEqual(len(res2),len(energies))
        self.assertEqual(len(res2),3)
        for i,e in enumerate(energies):
            r = catima.dedx_from_range(p(e), graphite)
            self.assertAlmostEqual(res2[i], r, 0.1)

    def test_layer_calculation(self):
        p = catima.Projectile(12,6)
        water = catima.get_material(catima.material.Water)
        water.thickness(10.0)
        graphite = catima.get_material(6)
        graphite.thickness(1.0)
        graphite.density(2.0)

        mat = catima.Layers()
        mat.add(water)
        mat.add(graphite)
        res = catima.calculate_layers(p(1000),mat)
        self.assertEqual(len(res.results),2)
        self.assertAlmostEqual(res.total_result.Eout,926.3,1)
        self.assertAlmostEqual(res.total_result.sigma_a,0.00269,1)
        self.assertAlmostEqual(res.Eout,926.3,1)
        self.assertAlmostEqual(res.sigma_a,0.00269,4)
        self.assertAlmostEqual(res.tof,0.402,2)
        self.assertAlmostEqual(res.Eloss,884,0)

        self.assertAlmostEqual(res[0].Eout,932.24,0)
        self.assertAlmostEqual(res[1].Eout,926.3,0)
        self.assertAlmostEqual(res[0].sigma_a,0.0025,3)
        self.assertAlmostEqual(res[1].sigma_a,0.000774,4)
        self.assertAlmostEqual(res[0].range,107.1,0)
        self.assertAlmostEqual(res[1].range,111.3,0)
    
    def test_energy_table(self):
        table = catima.get_energy_table()
        self.assertEqual(table[0],catima.energy_table(0))
        self.assertEqual(table[10],catima.energy_table(10))
        self.assertEqual(len(table),catima.max_datapoints)
        
    def test_storage(self):
        p = catima.Projectile(12,6)
        water = catima.get_material(catima.material.Water)
        water.thickness(10.0)
        graphite = catima.get_material(6)
        graphite.thickness(1.0)
        
        data = catima.get_data(p, water)
        et = catima.get_energy_table()
        
        self.assertEqual(len(data),3)
        self.assertEqual(len(data[0]),len(et))
        
        res = catima.calculate(p(et[10]),water)
        self.assertAlmostEqual(res.range,data[0][10],6)
        self.assertAlmostEqual(catima.range(p,water),data[0][10],6)
        #self.assertAlmostEqual(catima.domega2de(p,water),data[1][10],6)
        
        res = catima.calculate(p(et[100]),water)
        self.assertAlmostEqual(res.range,data[0][100],6)
        self.assertAlmostEqual(catima.range(p,water),data[0][100],6)
        #self.assertAlmostEqual(catima.domega2de(p,water),data[1][100],6)
        
        res = catima.calculate(p(et[200]),water)
        self.assertAlmostEqual(res.range,data[0][200],6)
        self.assertAlmostEqual(catima.range(p,water),data[0][200],6)
        #self.assertAlmostEqual(catima.domega2de(p,water),data[1][200],6)
        
        res = catima.calculate(p(et[401]),water)
        self.assertAlmostEqual(res.range,data[0][401],6)
        self.assertAlmostEqual(catima.range(p,water),data[0][401],6)
        #self.assertAlmostEqual(catima.domega2de(p,water),data[1][401],6)
        
    def test_python_storage_access(self):
        
        p = catima.Projectile(12,6)
        water = catima.get_material(catima.material.Water)
        water.thickness(10.0)
        graphite = catima.get_material(6)
        graphite.thickness(1.0)
        data = catima.get_data(p, water)
        self.assertEqual(catima.max_storage_data,60) # assuming 60, this has to be changed manually
        r = catima.storage_info()
        
        #self.assertAlmostEqual(catima.da2de(p,water,et[100]),data[2][100],6)
        #self.assertAlmostEqual(catima.da2de(p,water,et[400]),data[2][400],6)

if __name__ == "__main__":
    unittest.main()
