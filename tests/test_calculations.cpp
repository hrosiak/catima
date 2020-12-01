#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS
#include "doctest.h"
#include "testutils.h"
#include <math.h>
#include "catima/catima.h"
#include "catima/calculations.h"

using namespace std;

    TEST_CASE("nuclear stopping power"){
      catima::Target carbon{12.0107,6};
      catima::Projectile p{4.00151,2,2,1};
      
      double dif;
      p.T = 0.1/p.A; //0.1MeV
      dif = catima::dedx_n(p,carbon) - 14.27;
      CHECK( fabs(dif)< 1);
      
      p.T = 1/p.A; //1MeV
      dif = catima::dedx_n(p,carbon) - 2.161;
      CHECK( fabs(dif)< 0.1);
      
      p.T = 10/p.A; //10MeV
      dif = catima::dedx_n(p,carbon) - 0.2874;
      CHECK( fabs(dif) < 0.01);
      
      p.T = 100/p.A; //100MeV
      dif = catima::dedx_n(p,carbon) - 0.03455;
      CHECK( fabs(dif) < 0.001);
    }
    
    TEST_CASE("proton stopping power from srim"){
        catima::Projectile p{1,1,1,1};
        auto he = catima::get_material(2);
        auto carbon = catima::get_material(6);

        p.T = 1;
        CHECK( catima::sezi_dedx_e(p,he) == approx(283,1));
        p.T = 10;
        CHECK( catima::sezi_dedx_e(p,he) == approx(45.6,1));
        
        p.T = 30;        
        CHECK( catima::sezi_dedx_e(p,he) == approx(18.38,1));

        p.T = 1;        
        CHECK( catima::sezi_dedx_e(p,carbon) == approx(229.5,1));
        p.T = 10;
        CHECK( catima::sezi_dedx_e(p,carbon) == approx(40.8,1));
        p.T = 30;
        CHECK( catima::sezi_dedx_e(p,carbon) == approx(16.8,1));        
    }
    TEST_CASE("dedx, low energy, from sezi"){
        catima::Projectile p{4,2,2,1};
        auto carbon = catima::get_material(6);
        
        // He projectile TEST_CASE
        p.T = 1;
        CHECK( catima::sezi_dedx_e(p,carbon)+catima::dedx_n(p,carbon) == approx(922.06).R(0.0001) );
        p.T = 3;
        CHECK( catima::sezi_dedx_e(p,carbon)+catima::dedx_n(p,carbon) == approx(433.09).R(0.0001) );
        
        // C projectile TEST_CASE
        p.A = 12;
        p.Z = 6;
        p.T = 1;
        CHECK( catima::sezi_dedx_e(p,carbon)+catima::dedx_n(p,carbon) == approx( 5792.52).R(0.0001) );
        p.T = 9.9;
        CHECK( catima::sezi_dedx_e(p,carbon)+catima::dedx_n(p,carbon) == approx(1485.36).R(0.0001) );
        
    }
    TEST_CASE("LS check: deltaL values"){
        catima::Projectile p{238,92,92,1};
        
        p.T = 93.1494;
        CHECK(catima::bethek_lindhard(p)== approx(-0.5688,0.0001));
        
        p.T = 380.9932;
        CHECK(catima::bethek_lindhard(p)== approx(0.549199,0.0001));
        
        p.T = 995.368987;
        CHECK(catima::bethek_lindhard(p)== approx(1.106649).R(0.001) );
        
        p.T = 2640.032566;
        CHECK(catima::bethek_lindhard(p)== approx(1.35314).R(0.001) );

        p.T = 6091.392448;
        CHECK(catima::bethek_lindhard(p)== approx(1.365643).R(0.001) );

        p.T = 37277.695445;
        CHECK(catima::bethek_lindhard(p)== approx(0.689662).R(0.001) );
    }
    
    TEST_CASE("LS check: straggling values"){
        catima::Projectile p{238,92,92,1};
        
        auto f = [&](){return catima::bethek_lindhard_X(p);};
        
        p.T = 93.1494;
        CHECK( f() == approx(1.56898).R(0.01) );
        
        p.T = 380.9932;
        CHECK( f() == approx(1.836008).R(0.01) );
        
        p.T = 996.9855;
        CHECK( f() == approx(1.836528).R(0.01) );
        
        p.T = 2794.4822;
        CHECK( f()== approx(1.768018).R(0.01) );

        for(double e:{2000,20000,200000, 9000000, 50000000})
            CHECK(catima::precalculated_lindhard_X(p(e)) >= 0.0);
    }
    TEST_CASE("ultrarelativistic corrections"){
        catima::Projectile p{238,92};
        catima::Target t{27,13};
        CHECK(catima::pair_production(p(1e3),t) == approx(0.0,1e-3));
        CHECK(catima::bremsstrahlung(p(1e3),t) == approx(0.0,1e-3));

        CHECK(catima::pair_production(p(1e6),t) == approx(1900,300));
        CHECK(catima::bremsstrahlung(p(1e6),t) == approx(170,20));
        CHECK(catima::pair_production(p(7e6),t) == approx(19000,5000));
        CHECK(catima::bremsstrahlung(p(7e6),t) == approx(6000,500));
    }
    TEST_CASE("dEdx for compounds"){
        catima::Projectile p{1,1,1,1000};
        catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });

        CHECK( catima::dedx(p(1000), water) == approx(2.23).R(5e-3));
        CHECK( catima::dedx(p(500), water) == approx(2.76).R(5e-3));
        CHECK( catima::dedx(p(9), water) == approx(51.17).R(5e-3));
    }
    TEST_CASE("dEdx from spline vs dEdx"){
        catima::Projectile p{238,92,92,1000};
        catima::Material graphite({
                {12.011,6,1},
                });
          graphite.density(2);

        auto res = catima::calculate(p(1000),graphite);
        CHECK(catima::dedx(p(1000), graphite) == approx(res.dEdxi).R(0.001) );
        
        res = catima::calculate(p,graphite,500);
        CHECK(catima::dedx(p(500), graphite) == approx(res.dEdxi).R(0.001) );
        
        res = catima::calculate(p,graphite,9);
        CHECK(catima::dedx(p(9), graphite) == approx(res.dEdxi).R(0.001) );
    }
   

    TEST_CASE("Eout test"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      catima::Material graphite;
            graphite.add_element(12,6,1);
            graphite.density(2.0);
            graphite.thickness(0.5);
      
      auto res = catima::calculate(p,graphite);
      CHECK( res.Eout == approx(997.07,01));
    }
  TEST_CASE("TOF test"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(1.0);

      catima::Material graphite;
            graphite.add_element(12,6,1);
            graphite.density(2.0);
            graphite.thickness(0.5);
      double dif;
      
      auto res = catima::calculate(p,water);
      dif = res.tof - 0.038;
      CHECK( fabs(dif) < 0.01);
    }
    TEST_CASE("angular scattering"){
      catima::Projectile p{1,1,1,158.6};
      catima::Material cu = catima::get_material(29);
          
      
      cu.thickness_cm(0.02963);
      auto res = catima::calculate(p,cu);
      CHECK( 1000*res.sigma_a == approx(7.2,0.5));

      cu.thickness_cm(0.2963);
      res = catima::calculate(p,cu);
      CHECK( 1000*res.sigma_a == approx(23.5,3));
      
      catima::Layers ll;
      cu.thickness_cm(0.2963/2.);
      ll.add(cu);
      ll.add(cu);
      auto res2 = catima::calculate(p,ll);
      CHECK( 1000*res2.total_result.sigma_a == approx(23.5,3));
      CHECK( 1000*res2.total_result.sigma_a == approx(1000*res.sigma_a,0.05));
      }

    TEST_CASE("displacement test"){
      catima::Projectile p{1,1,1,215};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      catima::Material water2({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water2.thickness_cm(9.6);
      water2.density(2.0);

      water.thickness_cm(9.6);
      auto res = catima::calculate(p,water);
      CHECK( res.sigma_x == approx(0.1,0.03));
      auto resb = catima::calculate(p,water2);
      CHECK( res.sigma_x > resb.sigma_x);

      water.thickness_cm(17.0);
      auto res17 = catima::calculate(p,water);
      CHECK( res17.sigma_x == approx(0.25,0.05));

      water.thickness_cm(29.4);
      auto res29 = catima::calculate(p,water);
      CHECK( res29.sigma_x == approx(0.66,0.08));

      catima::Layers ll;
      water.thickness_cm(9.6);
      ll.add(water);
      auto res2 = catima::calculate(p,ll);
      CHECK(res2.total_result.sigma_x == approx(0.1,0.03));
      water.thickness_cm(7.4);      
      ll.add(water);
      res2 = catima::calculate(p,ll);
      CHECK(ll.thickness_cm() == approx(17,0.00001));
      CHECK(res2.total_result.sigma_x == approx(0.25,0.05));
      CHECK(res2.total_result.sigma_x == approx(res17.sigma_x).R(0.03));

      catima::Layers l;
      water.thickness_cm(9.6/3);
      l.add(water);
      l.add(water);
      l.add(water);
      res2 = catima::calculate(p(215),l);
      CHECK(res2.total_result.sigma_x == approx(res.sigma_x).R(0.03));

      for(int n=2;n<10;n++){
        catima::Layers lll;
        water.thickness_cm(29.4/n);
        for(int i=0;i<n;i++)lll.add(water);        
        res2 = catima::calculate(p(215),lll);
        CHECK(res2.total_result.sigma_x == approx(res29.sigma_x).R(0.01));
        }
      

      /// position sigma is taken from range thickness if particle stoppped inside
      water.thickness_cm(30);
      auto r1 = catima::calculate(p,water);
      water.thickness_cm(40);
      auto r2 = catima::calculate(p,water);
      CHECK(r1.sigma_x == approx(r2.sigma_x).R(0.01));

    }
    TEST_CASE("result from stopped material"){
        catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(1000.0);
      auto res = catima::calculate(p,water);
      CHECK(res.Eout == 0.0);
      CHECK(res.Eloss == 1000*12);
      CHECK(res.sigma_E == 0.0);
      CHECK(res.sigma_a == 0.0);
      CHECK(res.sigma_r > 0.0);
      CHECK(res.dEdxo == 0.0);
      CHECK(res.tof == 0.0);

      catima::Layers mat;
      mat.add(water);
      auto res2= catima::calculate(p,mat);
      CHECK(res2.results.size() == 1);
      CHECK(res2.total_result.Eout == res2.results[0].Eout);
      CHECK(res2.total_result.Eout == 0.0);
      CHECK(res2.total_result.Eloss == 1000*12);
      CHECK(res2.total_result.sigma_E == 0.0);
      CHECK(res2.total_result.sigma_a == 0.0);
      CHECK(res2.total_result.tof == 0.0);
    }
    TEST_CASE("constant results from material"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(10.0);
      auto res = catima::calculate(p,water);
      auto res2 = catima::calculate(p,water);
      CHECK(res.Eout == res2.Eout);
      CHECK(res.Eloss == res2.Eloss);
      CHECK(res.sigma_E == res2.sigma_E);
      CHECK(res.sigma_a == res2.sigma_a);
      CHECK(res.sigma_r == res2.sigma_r);
      CHECK(res.dEdxo == res2.dEdxo);
      CHECK(res.tof == res2.tof);
    }
    TEST_CASE("simplified calculation"){
        catima::Projectile p{12,6,6,1000};  
        catima::Material graphite({
                {12.011,6,1},
                });
        graphite.density(2.0).thickness(1.0);
        auto res1 = catima::calculate(p,graphite);
        auto res2 = catima::calculate(12,6,1000,12.011,6,1.0,2.0);
        CHECK(res1.Eout == res2.Eout);
        CHECK(res1.Eloss == res2.Eloss);
        CHECK(res1.sigma_E == res2.sigma_E);
        CHECK(res1.sigma_a == res2.sigma_a);
        CHECK(res1.sigma_r == res2.sigma_r);
        CHECK(res1.dEdxo == res2.dEdxo);
        CHECK(res1.tof == res2.tof);

        auto ra = catima::angular_straggling_from_E(p,res1.Ein,res1.Eout,graphite);
        CHECK(res1.sigma_a == ra);

        auto re = catima::energy_straggling_from_E(p,res1.Ein,res1.Eout,graphite);
        CHECK(res1.sigma_E == re);

        auto eo1 = catima::energy_out(p(1000),graphite);
        CHECK(res1.Eout == eo1);

        auto de1 = catima::dedx_from_range(p(1000),graphite);
        CHECK(res1.dEdxi == de1);

    }
    TEST_CASE("multilayer basic"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(10.0);

      catima::Material graphite({
                {12.011,6,1},
                });
      graphite.density(2.0).thickness(1.0);

        catima::Layers mat;
        mat.add(water);
        mat.add(graphite);
        
        auto res = catima::calculate(p(1000),mat);
        CHECK(res.total_result.Eout == approx(926.3,0.1));
        CHECK(res.total_result.sigma_a == approx(0.00269).R(0.05));
        CHECK(res.total_result.tof == approx(0.402).R(0.001));
        CHECK(res.total_result.Eloss == approx(884.2,1.0));
        //CHECK(rcompare(res.total_result.sigma_E,0.7067,0.2));
        CHECK(res.results[0].Eout == approx(932.24,0.1));
        CHECK(res.results[0].sigma_a == approx(0.00258).R(0.05));
        CHECK(res.results[0].range == approx(107.163,0.1));
        CHECK(res.results[1].Eout == approx(926.3,0.1));
        CHECK(res.results[1].sigma_a == approx(0.000774).R(0.05));
        CHECK(res.results[1].range == approx(111.3,0.1));

        auto res0 = catima::calculate(p(1000),water);
        CHECK(res0.Eout == res.results[0].Eout);
        CHECK(res0.sigma_a == res.results[0].sigma_a);
        CHECK(res0.sigma_E == res.results[0].sigma_E);
        CHECK(res0.sigma_r == res.results[0].sigma_r);
        CHECK(res0.tof == res.results[0].tof);

    }
    TEST_CASE("default material calculations"){
        catima::Projectile p{12,6,6,350};
        auto air = catima::get_material(catima::material::Air);
        air.thickness(0.500);
        auto res = catima::calculate(p(350),air);
        CHECK(res.Eout == approx(345.6).epsilon(1.0));
        CHECK(res.sigma_a == approx(0.0013).epsilon(1e-4));
        CHECK(res.sigma_E == approx(0.12).epsilon(1e-3));
        CHECK(res.dEdxi == approx(103.5).epsilon(1e-1));

        res = catima::calculate(p(150),air);
        CHECK(res.dEdxi == approx(173.6).epsilon(1e0));
        res = catima::calculate(p(1000),air);
        CHECK(res.dEdxi == approx(70.69).epsilon(1e-0));

        
        auto water = catima::get_material(catima::material::Water);
        auto res2 = catima::calculate(p(600),water,600);
        CHECK(res2.dEdxi == approx(92.5).epsilon(2));
    }
    TEST_CASE("z_eff"){
        using namespace catima;
        Projectile p_u(238,92);
        Target t;
        t.Z = 13;
        Config c;

        c.z_effective = z_eff_type::pierce_blann;
        CHECK(z_eff_Pierce_Blann(92,beta_from_T(5000.)) == approx(91.8).epsilon(0.2));
        CHECK(z_eff_Pierce_Blann(92,beta_from_T(5000.)) == z_effective(p_u(5000.),t,c));
        
        CHECK(z_eff_Winger(92,0.99,6) == approx(91.8).epsilon(0.5));
        CHECK(z_eff_Winger(92,beta_from_T(5000.),13) == approx(91.8).epsilon(0.2));
        c.z_effective = z_eff_type::winger;
        CHECK(z_eff_Winger(92,beta_from_T(5000.),13) == z_effective(p_u(5000.),t,c));
        
        CHECK(z_eff_Schiwietz(92,0.99,6) == approx(91.8).epsilon(0.5));
        c.z_effective = z_eff_type::schiwietz;
        CHECK(z_eff_Schiwietz(92,beta_from_T(5000.),13) == z_effective(p_u(5000.),t,c));

        CHECK(z_eff_Hubert(92,1900,13) == approx(91.88).epsilon(0.1));
        c.z_effective = z_eff_type::hubert;
        CHECK(z_eff_Hubert(92,1900,13) == z_effective(p_u(1900.),t,c));

        #ifdef GLOBAL
        CHECK(z_eff_global(92,1900,13) == approx(91.88).epsilon(0.05));
        c.z_effective = z_eff_type::global;
        CHECK(z_eff_global(92,1900,13) == z_effective(p_u(1900.),t,c));
        CHECK(z_eff_global(92,1000,13) == approx(91.71).epsilon(0.05));
        CHECK(z_eff_global(92,500,13) == approx(91.22).epsilon(0.1));
        CHECK(z_eff_global(92,100,6) == approx(89.61).epsilon(0.2));
        //CHECK(z_eff_global(92,100,13) == approx(89.42).epsilon(0.1));
        //CHECK(z_eff_global(92,100,29) == approx(88.37).epsilon(0.1));
        //CHECK(z_eff_global(92,50,13) == approx(85.94).epsilon(0.1));
        CHECK(z_eff_global(92,2001,13) == approx(92.0).epsilon(0.01));
        CHECK(z_eff_global(92,2000,13) == approx(92.0).epsilon(0.2));

        CHECK(z_eff_atima14(92,1900,13) == approx(91.88).epsilon(0.05));
        c.z_effective = z_eff_type::atima14;
        CHECK(z_eff_atima14(92,1900,13) == z_effective(p_u(1900.),t,c));
        #endif
    }
    TEST_CASE("vector_inputs"){
        catima::Projectile p{12,6,6,1000};
        catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
        catima::Material graphite;
        graphite.add_element(12,6,1);
        graphite.density(2.0);
        graphite.thickness(0.5);
      
      auto res = catima::calculate(p,graphite);
      CHECK( res.Eout == approx(997.07,01));

      std::vector<double> energies{100,500,1000};
      auto res2 = catima::energy_out(p,energies, graphite);
      CHECK(res2.size()==energies.size());
      CHECK(res2[2] == approx(997.07,01));
      CHECK(res2[0] == approx(catima::energy_out(p(energies[0]),graphite),0.1));
      CHECK(res2[1] == approx(catima::energy_out(p(energies[1]),graphite),0.1));
      CHECK(res2[2] == approx(catima::energy_out(p(energies[2]),graphite),0.1));

      auto res3 = catima::dedx_from_range(p,energies,graphite);
      CHECK(res3.size()==energies.size());
      CHECK(res3[0] == approx(catima::dedx_from_range(p(energies[0]),graphite),0.1));
      CHECK(res3[1] == approx(catima::dedx_from_range(p(energies[1]),graphite),0.1));
      CHECK(res3[2] == approx(catima::dedx_from_range(p(energies[2]),graphite),0.1));
    }
    TEST_CASE("constants"){
        using namespace catima;
        CHECK(0.1*hbar*c_light/atomic_mass_unit == approx(0.21183,0.0001));
        CHECK(16.0*dedx_constant*electron_mass*fine_structure/(atomic_mass_unit*3.0*4.0*PI) == approx(5.21721169334564e-7).R(1e-3));
        }
    
