#include "lest.hpp"
#include "testutils.h"
#include <math.h>
using namespace std;
using catima::approx;
#include "catima/catima.h"
#include "catima/calculations.h"

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
    
    CASE("nuclear stopping power"){
      catima::Target carbon{12.0107,6};
      catima::Projectile p{4.00151,2,2,1};
      
      double dif;
      p.T = 0.1/p.A; //0.1MeV
      dif = catima::dedx_n(p,carbon) - 14.27;
      EXPECT( fabs(dif)< 1);
      
      p.T = 1/p.A; //1MeV
      dif = catima::dedx_n(p,carbon) - 2.161;
      EXPECT( fabs(dif)< 0.1);
      
      p.T = 10/p.A; //10MeV
      dif = catima::dedx_n(p,carbon) - 0.2874;
      EXPECT( fabs(dif) < 0.01);
      
      p.T = 100/p.A; //100MeV
      dif = catima::dedx_n(p,carbon) - 0.03455;
      EXPECT( fabs(dif) < 0.001);
    },
    
    CASE("proton stopping power from srim"){
        catima::Projectile p{1,1,1,1};
        catima::Target he{4.002600,2};
        catima::Target carbon{12.0107,6};
        double dif,dif2;
        
        
        dif = catima::sezi_p_se(1,he) - 283;
        EXPECT( fabs(dif)< 1);
        p.T = 1;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
        
        dif = catima::sezi_p_se(10,he) - 45.6;
        EXPECT( fabs(dif)< 1);
        p.T = 10;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
        
        dif = catima::sezi_p_se(30,he) - 18.38;
        EXPECT( fabs(dif)< 1);
        p.T = 30;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
        
        dif = catima::sezi_p_se(1,carbon) - 229.5;
        EXPECT( fabs(dif)< 1);
        p.T = 1;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
        
        dif = catima::sezi_p_se(10,carbon) - 40.8;
        EXPECT( fabs(dif)< 1);
        p.T = 10;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
        
        dif = catima::sezi_p_se(30,carbon) - 16.8;
        EXPECT( fabs(dif)< 1);
        p.T = 30;
        dif2 = catima::sezi_p_se(p.T,he) - catima::sezi_dedx_e(p,he);
        EXPECT( fabs(dif2)< 0.000001);
    },
    CASE("dedx, low energy, from sezi"){
        catima::Projectile p{4,2,2,1};
        catima::Target carbon{12.0107,6};
        double dif;
        double exp;
        
        // He projectile case
        p.T = 1;
        EXPECT( rcompare( catima::dedx(p,carbon), 922.06, 0.0001) );
        p.T = 3;
        EXPECT( rcompare( catima::dedx(p,carbon), 433.09, 0.0001) );
        
        // C projectile case
        p.A = 12;
        p.Z = 6;
        p.T = 1;
        EXPECT( rcompare( catima::dedx(p,carbon), 5792.52, 0.0001) );
        p.T = 9.9;
        EXPECT( rcompare( catima::dedx(p,carbon), 1485.36, 0.0001) );
        
    },
    CASE("LS check: deltaL values"){
        catima::Projectile p{238,92,92,1};
        catima::Target carbon{12.0107,6}; 
        
        p.T = 93.1494;
        EXPECT( rcompare( catima::bethek_lindhard(p), -0.5688, 0.01) );
        
        p.T = 380.9932;
        EXPECT( rcompare( catima::bethek_lindhard(p), 0.549199, 0.02) );
        
        p.T = 996.9855;
        EXPECT( rcompare( catima::bethek_lindhard(p), 1.0732, 0.03) );
        
        p.T = 2794.4822;
        EXPECT( rcompare( catima::bethek_lindhard(p), 1.358964, 0.02) );
    },
    
    CASE("LS check: straggling values"){
        catima::Projectile p{238,92,92,1};
        catima::Target carbon{12.0107,6}; 
        
        auto f = [&](){return catima::bethek_lindhard_X(p);};
        
        p.T = 93.1494;
        EXPECT( rcompare( f(), 1.56898, 0.01) );
        
        p.T = 380.9932;
        EXPECT( rcompare( f(), 1.836008, 0.01) );
        
        p.T = 996.9855;
        EXPECT( rcompare( f(), 1.836528, 0.03) );
        
        p.T = 2794.4822;
        EXPECT( rcompare( f(), 1.768018, 0.02) );
    },
    CASE("dEdx for compounds"){
        catima::Projectile p{1,1,1,1000};
        catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
        double dif;
        dif = catima::dedx(p,1000, water) - 2.23;
        EXPECT( fabs(dif) < 0.002);
        
        dif = catima::dedx(p,500,water) - 2.76;
        EXPECT( fabs(dif) < 0.005);
        
        dif = catima::dedx(p,9,water) - 51.17;
        EXPECT( fabs(dif) < 0.005);
    },
    CASE("dEdx from spline vs dEdx"){
        catima::Projectile p{238,92,92,1000};
        catima::Material graphite({
                {12.011,6,1},
                });
        
        double dif;
        auto res = catima::calculate(p(1000),graphite);
        dif = (catima::dedx(p,1000, graphite) - res.dEdxi)/res.dEdxi;
        EXPECT(dif == approx(0).epsilon(0.001) );
        
        res = catima::calculate(p,graphite,500);
        dif = catima::dedx(p,500,graphite) - res.dEdxi;
        EXPECT(dif/res.dEdxi == approx(0).epsilon(0.001) );
        
        res = catima::calculate(p,graphite,9);
        dif = catima::dedx(p,9,graphite) - res.dEdxi;
        EXPECT(dif/res.dEdxi == approx(0).epsilon(0.001) );
    },
    
//    CASE("dOmega2dx Ferisov test"){
    
    //},
    
    CASE("Eout test"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      catima::Material graphite;
            graphite.add_element(12,6,1);
            graphite.density(2.0);
            graphite.thickness(0.5);
      double dif;
      
      auto res = catima::calculate(p,graphite);
      dif = res.Eout - 997.077;
      EXPECT( fabs(dif) < 0.01);
    },
  CASE("TOF test"){
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
      EXPECT( fabs(dif) < 0.01);
    },
    CASE("result from stopped material"){
        catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(1000.0);
      auto res = catima::calculate(p,water);
      EXPECT(res.Eout == 0.0);
      EXPECT(res.Eloss == 1000*12);
      EXPECT(res.sigma_E == 0.0);
      EXPECT(res.sigma_a == 0.0);
      EXPECT(res.sigma_r > 0.0);
      EXPECT(res.dEdxo == 0.0);
      EXPECT(res.tof == 0.0);

      catima::Layers mat;
      mat.add(water);
      auto res2= catima::calculate(p,mat);
      EXPECT(res2.results.size() == 1);
      EXPECT(res2.total_result.Eout == res2.results[0].Eout);
      EXPECT(res2.total_result.Eout == 0.0);
      EXPECT(res2.total_result.Eloss == 1000*12);
      EXPECT(res2.total_result.sigma_E == 0.0);
      EXPECT(res2.total_result.sigma_a == 0.0);
      EXPECT(res2.total_result.tof == 0.0);
    },
    CASE("constant results from material"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1.00794,1,2},
                {15.9994,8,1}
                });
      water.density(1.0);
      water.thickness(10.0);
      auto res = catima::calculate(p,water);
      auto res2 = catima::calculate(p,water);
      EXPECT(res.Eout == res2.Eout);
      EXPECT(res.Eloss == res2.Eloss);
      EXPECT(res.sigma_E == res2.sigma_E);
      EXPECT(res.sigma_a == res2.sigma_a);
      EXPECT(res.sigma_r == res2.sigma_r);
      EXPECT(res.dEdxo == res2.dEdxo);
      EXPECT(res.tof == res2.tof);
    },
    CASE("simplified calculation"){
        catima::Projectile p{12,6,6,1000};  
        catima::Material graphite({
                {12.011,6,1},
                });
        graphite.density(2.0).thickness(1.0);
        auto res1 = catima::calculate(p,graphite);
        auto res2 = catima::calculate(12,6,1000,12.011,6,1.0,2.0);
        EXPECT(res1.Eout == res2.Eout);
        EXPECT(res1.Eloss == res2.Eloss);
        EXPECT(res1.sigma_E == res2.sigma_E);
        EXPECT(res1.sigma_a == res2.sigma_a);
        EXPECT(res1.sigma_r == res2.sigma_r);
        EXPECT(res1.dEdxo == res2.dEdxo);
        EXPECT(res1.tof == res2.tof);

        auto ra = catima::angular_straggling_from_E(p,res1.Ein,res1.Eout,graphite);
        EXPECT(res1.sigma_a == ra);

        auto re = catima::energy_straggling_from_E(p,res1.Ein,res1.Eout,graphite);
        EXPECT(res1.sigma_E == re);

        auto eo1 = catima::energy_out(p,1000,graphite);
        EXPECT(res1.Eout == eo1);

        auto de1 = catima::dedx_from_range(p,1000,graphite);
        EXPECT(res1.dEdxi == de1);

    },
    CASE("multilayer basic"){
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
        EXPECT(rcompare(res.total_result.Eout,926.3,0.01));
        EXPECT(rcompare(res.total_result.sigma_a,0.00269,0.1));
        EXPECT(rcompare(res.total_result.tof,0.402,0.01));
        EXPECT(rcompare(res.total_result.Eloss,884.218,0.01));
        //EXPECT(rcompare(res.total_result.sigma_E,0.7067,0.2));
        EXPECT(rcompare(res.results[0].Eout,932.24,0.01));
        EXPECT(rcompare(res.results[0].sigma_a,0.00258,0.1));
        EXPECT(rcompare(res.results[0].range,107.163,0.01));
        EXPECT(rcompare(res.results[1].Eout,926.3,0.01));
        EXPECT(rcompare(res.results[1].sigma_a,0.000774,0.1));
        EXPECT(rcompare(res.results[1].range,110.715,0.01));

        auto res0 = catima::calculate(p(1000),water);
        EXPECT(res0.Eout == res.results[0].Eout);
        EXPECT(res0.sigma_a == res.results[0].sigma_a);
        EXPECT(res0.sigma_E == res.results[0].sigma_E);
        EXPECT(res0.sigma_r == res.results[0].sigma_r);
        EXPECT(res0.tof == res.results[0].tof);

    },
    CASE("default material calculations"){
        catima::Projectile p{12,6,6,350};
        auto air = catima::get_material(catima::material::Air);
        air.thickness(0.500);
        auto res = catima::calculate(p(350),air);
        EXPECT(res.Eout == approx(345.6).epsilon(1.0));
        EXPECT(res.sigma_a == approx(0.0013).epsilon(1e-4));
        EXPECT(res.sigma_E == approx(0.12).epsilon(1e-3));
        
        auto water = catima::get_material(catima::material::Water);
        auto res2 = catima::calculate(p(600),water,600);
        EXPECT(res2.dEdxi == approx(92.5).epsilon(2));
    },
    CASE("z_eff"){
        using namespace catima;
        Projectile p_u(238,92);
        Target t;
        t.Z = 13;
        Config c;

        EXPECT(z_eff_Pierce_Blann(92,beta_from_T(5000.)) == approx(91.8).epsilon(0.2));
        EXPECT(z_eff_Pierce_Blann(92,beta_from_T(5000.)) == z_effective(p_u(5000.),t,c));
        
        EXPECT(z_eff_Winger(92,0.99,6) == approx(91.8).epsilon(0.5));
        EXPECT(z_eff_Winger(92,beta_from_T(5000.),13) == approx(91.8).epsilon(0.2));
        c.z_effective = z_eff_type::winger;
        EXPECT(z_eff_Winger(92,beta_from_T(5000.),13) == z_effective(p_u(5000.),t,c));
        
        EXPECT(z_eff_Schiwietz(92,0.99,6) == approx(91.8).epsilon(0.5));
        c.z_effective = z_eff_type::schiwietz;
        EXPECT(z_eff_Schiwietz(92,beta_from_T(5000.),13) == z_effective(p_u(5000.),t,c));

        EXPECT(z_eff_Hubert(92,1900,13) == approx(91.88).epsilon(0.1));
        c.z_effective = z_eff_type::hubert;
        EXPECT(z_eff_Hubert(92,1900,13) == z_effective(p_u(1900.),t,c));

        #ifdef GLOBAL
        EXPECT(z_eff_global(92,1900,13) == approx(91.88).epsilon(0.05));
        c.z_effective = z_eff_type::global;
        EXPECT(z_eff_global(92,1900,13) == z_effective(p_u(1900.),t,c));
        EXPECT(z_eff_global(92,1000,13) == approx(91.71).epsilon(0.05));
        EXPECT(z_eff_global(92,500,13) == approx(91.22).epsilon(0.1));
        EXPECT(z_eff_global(92,100,6) == approx(89.61).epsilon(0.2));
        //EXPECT(z_eff_global(92,100,13) == approx(89.42).epsilon(0.1));
        //EXPECT(z_eff_global(92,100,29) == approx(88.37).epsilon(0.1));
        //EXPECT(z_eff_global(92,50,13) == approx(85.94).epsilon(0.1));
        EXPECT(z_eff_global(92,2001,13) == approx(92.0).epsilon(0.01));
        EXPECT(z_eff_global(92,2000,13) == approx(92.0).epsilon(0.2));

        EXPECT(z_eff_atima14(92,1900,13) == approx(91.88).epsilon(0.05));
        c.z_effective = z_eff_type::atima14;
        EXPECT(z_eff_atima14(92,1900,13) == z_effective(p_u(1900.),t,c));
        #endif
    }
    
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


