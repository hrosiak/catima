#include "lest.hpp"
#include <math.h>
#include "catima/catima.h"
//#include "nucdata.h"

using namespace std;
using lest::approx;
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
    
    CASE("LS generated is equal to calculated"){
      catima::Projectile p;
      double a,b;
      
      p.A = 238;
      p.Z = 92;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard(p);
          b = catima::precalculated_lindhard(p);
          EXPECT(a==approx(b).epsilon(0.001));
      }
      
      p.A = 220;
      p.Z = 92;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard(p);
          b = catima::precalculated_lindhard(p);
          EXPECT(a==approx(b).epsilon(0.01));
      }
      
      p.A = 250;
      p.Z = 92;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard(p);
          b = catima::precalculated_lindhard(p);
          EXPECT(a==approx(b).epsilon(0.01));
      }
      
      
      p.A = 200;
      p.Z = 76;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard(p);
          b = catima::precalculated_lindhard(p);
          EXPECT(a==approx(b).epsilon(0.01));
      }
      
      p.A = 100;
      p.Z = 50;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard(p);
          b = catima::precalculated_lindhard(p);
          EXPECT(a==approx(b).epsilon(0.01));
      }
      
    },
    
    CASE("LS X generated is equal to calculated"){
      catima::Projectile p;
      double a,b;
      
      p.A = 238;
      p.Z = 92;
      for(double e:{90.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard_X(p);
          b = catima::precalculated_lindhard_X(p);
          EXPECT(a==approx(b).epsilon(0.001));
      }
      
      p.A = 220;
      p.Z = 92;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard_X(p);
          b = catima::precalculated_lindhard_X(p);
          EXPECT(a==approx(b).epsilon(0.02));
      }
      
      p.A = 250;
      p.Z = 92;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard_X(p);
          b = catima::precalculated_lindhard_X(p);
          EXPECT(a==approx(b).epsilon(0.03));
      }
      
      
      p.A = 200;
      p.Z = 76;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard_X(p);
          b = catima::precalculated_lindhard_X(p);
          EXPECT(a==approx(b).epsilon(0.03));
      }
      
      p.A = 100;
      p.Z = 50;
      for(double e:{100.0,1000.0,5000.0,30000.0}){
          p.T = e;
          a = catima::bethek_lindhard_X(p);
          b = catima::precalculated_lindhard_X(p);
          EXPECT(a==approx(b).epsilon(0.03));
      }
      
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


