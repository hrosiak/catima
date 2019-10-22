#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
using namespace std;
#include "catima/catima.h"
#include "catima/material_database.h"

    TEST_CASE("atima material basic tests"){
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

            CHECK(graphite.ncomponents()==1);
            CHECK(water.ncomponents()==2);

            CHECK(graphite.M()==12);
            CHECK(water.ncomponents()==2);
            CHECK(water.M()==18);
            CHECK(water2.ncomponents()==2);
            CHECK(water2.M()==18);
            CHECK(graphite.I()==0.0);

            graphite.add_element(18,40,1);
            CHECK(graphite.ncomponents()==2);
            CHECK_FALSE(graphite.M()==approx(12,0.1));

            CHECK(water==water2);
            CHECK(!(water==graphite));
            water.density(1.0);
            water.thickness(1.0);
            CHECK(water.thickness()==approx(1.0,0.0001));
            water.thickness_cm(1.0);
            CHECK(water.thickness()==approx(1.0,0.0001));
            water.thickness_cm(2.0);
            CHECK(water.thickness()==approx(2.0,0.0001));
    }
    TEST_CASE("Material automatic atomic weight"){
        catima::Material water({{0,1,2},{0,8,1}});
        catima::Material graphite(0,6);
        CHECK(water.get_element(0).A == 1.00794);
        CHECK(water.get_element(1).A == 15.9994);
        CHECK(graphite.get_element(0).A == 12.0107);
        CHECK(water.M()>16);
        CHECK(graphite.M()>12);
    }
    TEST_CASE("default materials"){
        catima::Material m = catima::get_material(6);
        CHECK(m.get_element(0).A == 12.0107);
        CHECK(m.get_element(0).Z == 6);
        CHECK(m.density() == 2.0);
        CHECK(m.M() == 12.0107);

        m = catima::get_material(catima::material::Water);
        CHECK(m.get_element(0).A == 1.00794);
        CHECK(m.get_element(0).Z == 1);
        CHECK(m.get_element(1).A == 15.9994);
        CHECK(m.get_element(1).Z == 8);
        CHECK(m.density() == 1.0);
    }
    TEST_CASE("Layers"){
        catima::Material water2;
        water2.add_element(1,1,2);
        water2.add_element(16,8,1);
        water2.density(1.0).thickness(2.0);

        catima::Material graphite;
        graphite.add_element(12,6,1);
        graphite.density(1.8).thickness(1.0);

        catima::Layers detector1;
        CHECK(detector1.num() == 0);
        detector1.add(graphite);
        CHECK(detector1.num() == 1);
        detector1.add(water2);
        detector1.add(graphite);
        CHECK(detector1.num() == 3);
        // check correct density and thickness
        CHECK(detector1[0].density()==1.8);
        CHECK(detector1[0].thickness()==1.0);
        CHECK(detector1[1].density()==1.0);
        CHECK(detector1[1].thickness()==2.0);
        CHECK(detector1[0].get_element(0).Z == 6);
        CHECK(detector1[0].get_element(0).A == 12);
        CHECK(detector1[1].get_element(0).A == 1);
        CHECK(detector1[1].get_element(0).Z == 1);
        CHECK(detector1[2].density() == 1.8);
        detector1[1].density(1.2);
        CHECK(detector1[1].density()==1.2);

        catima::Layers detector2;
        detector2 = detector1;
        CHECK(detector2.num() == 3);

        detector2.add(water2);
        detector2.add(graphite);
        CHECK(detector2.num() == 5);

        catima::Layers focal_material = detector1 + detector2;
        CHECK(focal_material.num() == 8);

    }
    TEST_CASE("basic projectile tests"){
      catima::Projectile p{12,6,6,1000};
      CHECK(p.A==12);
      CHECK(p.Z==6);
      CHECK(p.Q==6);
      CHECK(p.T==1000);

      catima::Projectile p2(12,6);
      CHECK(p.A==12);
      CHECK(p2.Z==6);
      CHECK(p2.Q==6);
      CHECK(p2.T==0);
      p2(1000);
      CHECK(p2.T==1000);
      p2(1000)(500);
      CHECK(p2.T==500);

      catima::Projectile p3(12,6,5);

      CHECK(p==p2);
      CHECK( !(p==p3));
    }
    TEST_CASE("basic config test"){
      catima::Config c1;
      catima::Config c2;
      catima::Config c3{catima::z_eff_type::none};
      catima::Config c4;
      CHECK(c1.z_effective == catima::default_config.z_effective);

      CHECK(c1==c2);
      CHECK( !(c1==c3));
      CHECK(c1==c4);

      c4.z_effective = catima::z_eff_type::global;
      CHECK(!(c1==c4));
      auto c5 = c4;
      CHECK(c4==c5);

      c4.z_effective = catima::z_eff_type::hubert;
      CHECK(!(c4==c5) );
      CHECK(!(c4==c1));
      c4.z_effective = catima::default_config.z_effective;
      CHECK(!(c5==c4));
      CHECK((c1==c4));
    }
    TEST_CASE("constructors test"){
        catima::Material mat2(12,6,2.5,0.1);
        catima::Material mat3(12.01,6);
        catima::Material mat4({{12.01,6,1.0}});
        catima::Material mat5({
                        {12.01, 6, 1},
                        {16.00, 8, 2}
                        });
        CHECK(mat2.ncomponents()==1);
        CHECK(mat3.ncomponents()==1);
        CHECK(mat3.get_element(0).A==12.01);
        CHECK(mat4.ncomponents()==1);
        CHECK(mat4.get_element(0).A==12.01);
        CHECK(mat5.ncomponents()==2);
        CHECK(mat5.get_element(0).A==12.01);
        CHECK(mat5.get_element(0).Z==6);
        CHECK(mat5.get_element(1).A==16.0);
        CHECK(mat5.get_element(1).stn==2);

        catima::Material mat6;
        mat6 = mat5;
        CHECK(mat5==mat6);
        CHECK(mat5.density() == mat6.density());
        CHECK(mat5.ncomponents()==mat6.ncomponents());
        CHECK(mat5.get_element(0).A==mat6.get_element(0).A);
        CHECK(mat5.get_element(1).A==mat6.get_element(1).A);
        mat5.add_element(12,6,1);
        CHECK(mat5.ncomponents()==mat6.ncomponents()+1);

        // constructor with custom Ipot
        catima::Material water1({
                {1,1,2},
                {16,8,1}
                },1.0);
        catima::Material water2({
                {1,1,2},
                {16,8,1}
                },1.0, 78.0);
        CHECK(water1.ncomponents()==2);
        CHECK(water2.ncomponents()==2);
        CHECK(water1.density()==1.0);
        CHECK(water2.density()==1.0);
        CHECK(water1.I()==0.0);
        CHECK(water2.I()==78.0);
        CHECK_FALSE(water1==water2);

    }
    TEST_CASE("fraction vs stn init"){
      catima::Projectile p{12,6};
      catima::Material water1({
                        {0, 1, 0.111894},
                        {0, 8, 0.888106}
                        });
      catima::Material water2({
                        {0, 1, 2},
                        {0, 8, 1}
                        });
      water1.thickness(1.0);
      water2.thickness(1.0);
      auto res1 = catima::calculate(p(600),water1);
      auto res2 = catima::calculate(p(600),water2);
      CHECK(res1.dEdxi == approx(res2.dEdxi,0.001));
      CHECK(res1.range == approx(res2.range).R(1e-2));
      CHECK(res1.sigma_a == approx(res2.sigma_a).R(1e-2));
      CHECK(res1.sigma_r == approx(res2.sigma_r).R(1e-2));
    }
    TEST_CASE("fraction calculation"){
        catima::Material water1({
                        {0, 1, 0.111898},
                        {0, 8, 0.888102}
                        });
        catima::Material water2({
                        {0, 1, 2},
                        {0, 8, 1}
                        });

        auto air = catima::get_material(catima::material::Air);

        CHECK(water1.weight_fraction(0)==0.111898);
        CHECK(water2.weight_fraction(0)==approx(water1.weight_fraction(0)).R(1e-5));
        CHECK(water1.weight_fraction(1)==0.888102);
        CHECK(water2.weight_fraction(1)==approx(water1.weight_fraction(1)).R(1e-5));
        CHECK(water2.M()==approx(18).epsilon(0.1));

        CHECK(water1.M()==approx(6.0,0.1));
        CHECK(water2.M()==approx(18,0.1));

        CHECK(water1.molar_fraction(0)==approx(2.0/3.0).R(1e-5));
        CHECK(water2.molar_fraction(0)==approx(2.0).R(1e-5));
        CHECK(water1.molar_fraction(1)==approx(1.0/3.0).R(1e-5));
        CHECK(water2.molar_fraction(1)==approx(1.0).R(1e-5));
        CHECK(water1.molar_fraction(0)/water1.molar_fraction(1)==approx(2.0).R(1e-5));
        CHECK(water2.molar_fraction(0)/water2.molar_fraction(1)==approx(2.0).R(1e-5));


        catima::Material mat({12.0,6,1});
        CHECK(mat.M()==approx(12.0,0.001));
        CHECK(mat.weight_fraction(0)==approx(1.0).R(1e-6));

        //CHECK(air.M() == approx(28.97,0.1));

    }
    TEST_CASE("number density"){
        catima::Material c({12.0,6,1});
        auto water = catima::get_material(catima::material::Water);
        auto air = catima::get_material(catima::material::Air);
        water.density(0.9982);
        c.density(3.513);
        air.density(1.2041e-3);

        c.thickness_cm(1.0);
        CHECK(c.number_density()==approx(1.7662,0.01));
        CHECK(c.number_density_cm2()==approx(1.7662,0.01));
        CHECK(c.number_density(0)==approx(1.7662,0.01));
        CHECK(c.number_density_cm2(0)==approx(1.7662,0.01));
        CHECK(c.number_density(1)==0.0);
        CHECK(c.number_density_cm2(1)==0.0);
        c.thickness_cm(2.0);
        CHECK(c.number_density()==approx(1.7662,0.01));
        CHECK(c.number_density_cm2()==approx(2.0*1.7662,0.01));

        water.thickness_cm(1.0);
        CHECK(water.number_density()==approx(0.3336,0.001));
        CHECK(water.number_density_cm2()==approx(0.3336,0.001));
        CHECK(water.number_density(0)==approx(2*0.3336,0.001));
        CHECK(water.number_density_cm2(0)==approx(2*0.3336,0.001));
        CHECK(water.number_density(1)==approx(0.3336,0.001));
        CHECK(water.number_density_cm2(1)==approx(0.3336,0.001));
        water.thickness_cm(3.0);
        CHECK(water.number_density_cm2()==approx(3.0*0.3336,0.001));

        air.thickness_cm(1.0);
        CHECK(air.number_density(0)==approx(air.molar_fraction(0)*2*0.0002504,0.00001));
        CHECK(air.number_density(1)==approx(air.molar_fraction(1)*2*0.0002504,0.00001));
        CHECK(air.number_density(2)==approx(air.molar_fraction(2)*1*0.0002504,0.00001));
    }
