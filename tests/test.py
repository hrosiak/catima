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
        graphite = catima.get_material(6);
        graphite.thickness(0.5)
        p10 = catima.get_material(catima.material.P10);
        p10.thickness(0.01)

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
        self.assertEqual(mat[3],None)
        self.assertEqual(mat["a"],None)

if __name__ == "__main__":
    unittest.main()