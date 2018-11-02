#include "lest.hpp"
#include <math.h>
#include "testutils.h"
using namespace std;

#include "catima/catima.h"   
#include "catima/storage.h"   

using catima::approx;
const lest::test specification[] =
{
    
    CASE("datapoints equal operator"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1,1,2},
                {16,8,1}
                });
                
      catima::Material graphite({
                {12,6,1}
                });
      catima::Config c2;
      c2.z_effective = catima::z_eff_type::winger;
      catima::DataPoint a(p,water);
      catima::DataPoint b(p,water);
      catima::DataPoint c(p,graphite);
      catima::DataPoint d;
      catima::DataPoint e(p,water,c2);
      d = c;
      EXPECT(a == b);
      EXPECT(!(a==c));
      EXPECT(d == c);
      EXPECT(!(d==b));
      EXPECT(!(e==a));
      d = a;
      EXPECT(!(d==c));
      EXPECT(d==b);
      
    },
    
    CASE("storage add"){
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1,1,2},
                {16,8,1}
                });
                
      catima::Material graphite({
                {12,6,1}
                });
      
      catima::Config c1;
      c1.z_effective = catima::z_eff_type::global;
      
      catima::_storage.Reset();
      EXPECT(catima::_storage.get_index()==0);
      
      catima::_storage.Add(p,water);
      auto& dp = catima::_storage.Get(0);
      EXPECT(catima::_storage.get_index()==1);
      EXPECT(dp.p.A==12);
      EXPECT(dp.m.ncomponents()==2);
      catima::_storage.Add(p,water);
      auto& dp2 = catima::_storage.Get(1);
      EXPECT(catima::_storage.get_index()==1);
      EXPECT(dp2.p.A==0);
      EXPECT(dp2.m.ncomponents()==0);
      
      catima::_storage.Add(p,graphite);
      auto&  dp3 = catima::_storage.Get(1);
      EXPECT(catima::_storage.get_index()==2);
      EXPECT(dp3.p.A==12);
      EXPECT(dp3.m.ncomponents()==1);
      
      catima::_storage.Add(p,graphite);
      EXPECT(catima::_storage.get_index()==2);
      
      catima::_storage.Add(p,graphite, c1);
      EXPECT(catima::_storage.get_index()==3);
      
      catima::_storage.Add(p,graphite);
      EXPECT(catima::_storage.get_index()==3);
      c1.z_effective = catima::z_eff_type::atima14;
      catima::_storage.Add(p,graphite ,c1);
      EXPECT(catima::_storage.get_index()==4);
      
    },
    CASE("test maximum storage"){
      auto maxdata = catima::max_storage_data;
      catima::Projectile p{12,6,6,1000};
      catima::Material water({
                {1,1,2},
                {16,8,1}
                });
                
      catima::Material graphite({
                {12,6,1}
                });
      catima::_storage.Reset();      
      EXPECT(catima::_storage.get_index()==0);
      for(int i=1;i<maxdata+1;i++){
          catima::Projectile p1{2.0*i,(double)i,(double)i,1000};
          catima::_storage.Add(p1,graphite);
          EXPECT(catima::_storage.get_index()==i);
          EXPECT(catima::_storage.GetN()==maxdata);
      }
      EXPECT(catima::_storage.get_index()==maxdata);
      for(int i=1;i<maxdata-1;i++){
          catima::Projectile p1{2.0*i,(double)i,(double)i,1000};
          catima::_storage.Add(p1,water);
          EXPECT(catima::_storage.get_index()==i);
          EXPECT(catima::_storage.GetN()==maxdata);
      }

  
    },
    CASE("energy table"){
      double step = (catima::logEmax - catima::logEmin)/(catima::max_datapoints-1);
      EXPECT(catima::energy_table.step==step);
      EXPECT(catima::energy_table.values[0]==exp(M_LN10*(catima::logEmin)));
      EXPECT(catima::energy_table.values[1]==exp(M_LN10*(catima::logEmin+step)));
      EXPECT(catima::energy_table.values[2]==exp(M_LN10*(catima::logEmin+2.0*step)));
      EXPECT(catima::energy_table.values[3]==exp(M_LN10*(catima::logEmin+3.0*step)));
      EXPECT(catima::energy_table.values[4]==exp(M_LN10*(catima::logEmin+4.0*step)));
      EXPECT(catima::energy_table.values[5]==exp(M_LN10*(catima::logEmin+5.0*step)));
      EXPECT(catima::energy_table.values[catima::max_datapoints-1]==approx(exp(M_LN10*(catima::logEmax))).epsilon(1e-6));
    },
    CASE("indexing"){
      double val, dif;
      
      EXPECT(EnergyTable_index(catima::energy_table, 0.0)==-1);

      for(int i=0;i<catima::max_datapoints-1;i++){
          val = catima::energy_table.values[i];
          dif = catima::energy_table.values[i+1] - val;
          EXPECT(EnergyTable_index(catima::energy_table, val)==i);
          EXPECT(EnergyTable_index(catima::energy_table, val+0.5*dif)==i);
          EXPECT(catima::energy_table.index(val)==i);
          EXPECT(catima::energy_table.index(val+0.5*dif)==i);
      }
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


