#include "lest.hpp"
#include <math.h>
using namespace std;

#include "catima/catima.h"
#include "catima/material_database.h"

bool rcompare(double a, double b,double eps){
    if(fabs((a-b)/fabs(b))<eps){
      return true;
    }
    else{
      std::cout<<"\033[1;31m"<<a<<" == "<<b<<"\033[0m"<<std::endl;
      return false;
    }
      
}

const lest::test specification[] =
{
    CASE("atima material basic tests"){
        SETUP(""){
            catima::Material water({
                {1,1,2},
                {16,8,1}
                });
                
            catima::Material water2;
            water2.add_element(1,1,2);
            water2.add_element(16,8,1);
                
            catima::Material graphite;
            graphite.add_element(12,6,1);
            graphite.density(1.8);
            
            SECTION("component number check"){
                EXPECT(graphite.ncomponents()==1);
                EXPECT(water.ncomponents()==2);
                
                graphite.add_element(18,40,1);
                EXPECT(graphite.ncomponents()==2);      
                EXPECT(water.M()==18);
            }
            SECTION("Molar mass check"){
                EXPECT(graphite.M()==12);
                EXPECT(water.ncomponents()==2);
                EXPECT(water.M()==18);
                EXPECT(water2.ncomponents()==2);
                EXPECT(water2.M()==18);
            }
            SECTION("equal operator check"){
              EXPECT(water==water2);
              EXPECT(!(water==graphite));
            }
        }
    },
    CASE("Material automatic atomic weight"){
        catima::Material water({{0,1,2},{0,8,1}});
        catima::Material graphite(0,6);
        EXPECT(water.get_element(0).first.A == 1.00794);
        EXPECT(water.get_element(1).first.A == 15.9994);
        EXPECT(graphite.get_element(0).first.A == 12.0107);
        EXPECT(water.M()>16);
        EXPECT(graphite.M()>12);
    },
    CASE("default materials"){
        catima::Material m = catima::get_material(6);
        EXPECT(m.get_element(0).first.A == 12.0107);
        EXPECT(m.get_element(0).first.Z == 6);
        EXPECT(m.density() == 2.0);
        EXPECT(m.M() == 12.0107);

        m = catima::get_material(catima::material::WATER);
        EXPECT(m.get_element(0).first.A == 1.00794);
        EXPECT(m.get_element(0).first.Z == 1);
        EXPECT(m.get_element(1).first.A == 15.9994);
        EXPECT(m.get_element(1).first.Z == 8);
        EXPECT(m.density() == 1.0);
        EXPECT(m.M() > 16.0);


    },
    CASE("Layers"){
        catima::Material water2;
        water2.add_element(1,1,2);
        water2.add_element(16,8,1);
        water2.density(1.0).thickness(2.0);
                
        catima::Material graphite;
        graphite.add_element(12,6,1);
        graphite.density(1.8).thickness(1.0);

        catima::Layers detector1;
        EXPECT(detector1.num() == 0);
        detector1.add(graphite);
        EXPECT(detector1.num() == 1);
        detector1.add(water2);
        detector1.add(graphite);
        EXPECT(detector1.num() == 3);
        // check correct density and thickness
        EXPECT(detector1[0].density()==1.8);
        EXPECT(detector1[0].thickness()==1.0);
        EXPECT(detector1[1].density()==1.0);
        EXPECT(detector1[1].thickness()==2.0);
        EXPECT(detector1[0].get_element(0).first.Z == 6);
        EXPECT(detector1[0].get_element(0).first.A == 12);
        EXPECT(detector1[1].get_element(0).first.A == 1);
        EXPECT(detector1[1].get_element(0).first.Z == 1);
        EXPECT(detector1[2].density() == 1.8);
        detector1[1].density(1.2);
        EXPECT(detector1[1].density()==1.2);

        catima::Layers detector2;
        detector2 = detector1;
        EXPECT(detector2.num() == 3);

        detector2.add(water2);
        detector2.add(graphite);
        EXPECT(detector2.num() == 5);

        catima::Layers focal_material = detector1 + detector2;
        EXPECT(focal_material.num() == 8);
    
    },
    CASE("basic projectile tests"){
      catima::Projectile p{12,6,6,1000};
      EXPECT(p.A==12);
      EXPECT(p.Z==6);
      EXPECT(p.Q==6);
      EXPECT(p.T==1000);
      
      catima::Projectile p2(12,6);
      EXPECT(p.A==12);
      EXPECT(p2.Z==6);
      EXPECT(p2.Q==6);
      EXPECT(p2.T==0);
      p2(1000);
      EXPECT(p2.T==1000);
      p2(1000)(500);
      EXPECT(p2.T==500);

      catima::Projectile p3(12,6,5);
      
      EXPECT(p==p2);
      EXPECT( !(p==p3));
    },
    CASE("basic config test"){
      catima::Config c1;
      catima::Config c2;
      catima::Config c3{catima::z_eff_type::none};
      
      
      EXPECT(c1==c2);
      EXPECT( !(c1==c3));
    },
    CASE("constructors test"){
        catima::Material mat2(12,6,2.5,0.1);
        catima::Material mat3(12.01,6);
        catima::Material mat4({{12.01,6,1.0}});
        catima::Material mat5({
                        {12.01, 6, 1},
                        {16.00, 8, 2}
                        });
        EXPECT(mat2.ncomponents()==1);
        EXPECT(mat3.ncomponents()==1);
        EXPECT(mat3.get_element(0).first.A==12.01);
        EXPECT(mat4.ncomponents()==1);
        EXPECT(mat4.get_element(0).first.A==12.01);
        EXPECT(mat5.ncomponents()==2);
        EXPECT(mat5.get_element(0).first.A==12.01);
        EXPECT(mat5.get_element(0).first.Z==6);
        EXPECT(mat5.get_element(1).first.A==16.0);
        EXPECT(mat5.get_element(1).second==2);
        
        catima::Material mat6;
        mat6 = mat5;
        EXPECT(mat5==mat6);
        EXPECT(mat5.ncomponents()==mat6.ncomponents());
        EXPECT(mat5.get_element(0).first.A==mat6.get_element(0).first.A);
        EXPECT(mat5.get_element(1).first.A==mat6.get_element(1).first.A);
        mat5.add_element(12,6,1);
        EXPECT(mat5.ncomponents()==mat6.ncomponents()+1);

    },
    CASE("int and double stn check"){
      catima::Projectile p{1,1,1,1000};
      catima::Material mat1({
                        {12.01, 6, 1},
                        {16.00, 8, 1}
                        });
      catima::Material mat2({
                        {12.01, 6, 0.5},
                        {16.00, 8, 0.5}
                        });
      mat1.thickness(1.0);
      mat2.thickness(1.0);
      auto res1 = catima::calculate(p(1000),mat1);
      auto res2 = catima::calculate(p(1000),mat2);
      EXPECT(res1.dEdxi == res2.dEdxi);
      EXPECT(res1.range == res2.range);
      EXPECT(res1.sigma_a == res2.sigma_a);
      EXPECT(res1.sigma_r == res2.sigma_r);
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


