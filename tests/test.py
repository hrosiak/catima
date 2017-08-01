import sys
sys.path.insert(0,"../build")
import unittest
import catima

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
        self.assertEqual(p(),1000)
        p(500)
        self.assertEqual(p.T(),500)
        self.assertEqual(p(),500)

        p = catima.Projectile(238,92,90, T=100)
        self.assertEqual(p.T(),100)
    
    def test_Material(self):
        mat = catima.Material()
        mat.add_element(12,6,1)
        self.assertEqual(mat.ncomponents(),1)
        mat.add_element(1,1,2)
        self.assertEqual(mat.ncomponents(),2)

        mat2 = catima.Material([12.01,6,1])
        self.assertEqual(mat2.ncomponents(),1)
        self.assertEqual(mat2.molar_mass(),12.01)

        mat3 = catima.Material([12,6,1])
        self.assertEqual(mat3.ncomponents(),1)
        self.assertEqual(mat3.molar_mass(),12)

        water = catima.Material([[1,1,2],[16,8,1]])
        self.assertEqual(water.molar_mass(),18)

        mat2 = catima.Material([0,6,1])
        self.assertEqual(mat2.ncomponents(),1)
        self.assertAlmostEqual(mat2.molar_mass(),12,1)

        mat5 = catima.Material([0,6,1],density=1.9, thickness=0.5)
        self.assertEqual(mat5.ncomponents(),1)
        self.assertEqual(mat5.thickness(),0.5)
        self.assertEqual(mat5.density(),1.9)

        # copy
        mat3.density(1.8)
        matc = mat3.copy()
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

        m2 = catima.get_material(catima.material.WATER)
        self.assertEqual(m2.ncomponents(),2)
        self.assertAlmostEqual(m2.molar_mass(),18,1)
        self.assertAlmostEqual(m2.density(),1.0,1)

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
        mat[2].thickness(1.2)
        mat[2].density(1.9)
        self.assertAlmostEqual(mat.materials[2].thickness(),1.2,1)
        self.assertAlmostEqual(mat.materials[2].density(),1.9,1)
        #self.assertAlmostEqual(mat.materials[0].thickness(),0.5,1)
        #self.assertAlmostEqual(mat.materials[0].density(),2.0,1)
        self.assertEqual(mat[3],None)
        self.assertEqual(mat["a"],None)

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
        water = catima.get_material(catima.material.WATER)
        p = catima.Projectile(1,1)
        
        p(1000)
        res = catima.calculate(p,water)
        res2 = catima.dedx(p,water)
        self.assertAlmostEqual(res.dEdxi,2.23,1)
        self.assertAlmostEqual(res["dEdxi"],2.23,1)
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
    
    def test_eout(self):
        graphite = catima.get_material(6)
        graphite.thickness(0.5)
        p = catima.Projectile(12,6)
        res = catima.calculate(p(1000),graphite)
        res2 = catima.energy_out(p(1000),graphite)
        self.assertAlmostEqual(res.Eout,997.077,1)
        self.assertAlmostEqual(res["Eout"],997.077,1)
        self.assertAlmostEqual(res.Eout,res2,3)

    def test_layer_calculateion(self):
        p = catima.Projectile(12,6)
        water = catima.get_material(catima.material.WATER)
        water.thickness(10.0)
        graphite = catima.get_material(6)
        graphite.thickness(1.0)

        mat = catima.Layers()
        mat.add(water)
        mat.add(graphite)
        res = catima.calculate_layers(p(1000),mat)
        self.assertEqual(len(res.results),2)
        self.assertAlmostEqual(res.total_result.Eout,926.3,1)
        self.assertAlmostEqual(res.total_result.sigma_a,0.00269,1)
        self.assertAlmostEqual(res["Eout"],926.3,1)
        self.assertAlmostEqual(res["sigma_a"],0.00269,4)
        self.assertAlmostEqual(res["tof"],0.402,2)
        self.assertAlmostEqual(res["Eloss"],884,0)

        self.assertAlmostEqual(res[0]["Eout"],932,0)
        self.assertAlmostEqual(res[1]["Eout"],926,0)
        self.assertAlmostEqual(res[0]["sigma_a"],0.00258,4)
        self.assertAlmostEqual(res[1]["sigma_a"],0.000774,4)
        self.assertAlmostEqual(res[0]["range"],107.1,0)
        self.assertAlmostEqual(res[1]["range"],110.7,0)

if __name__ == "__main__":
    unittest.main()