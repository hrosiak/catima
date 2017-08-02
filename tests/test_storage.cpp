#include "lest.hpp"
#include <math.h>
using namespace std;

#include "catima/catima.h"   
#include "catima/storage.h"   
#include "catima/material_database.h"

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
      
      catima::DataPoint a(p,water);
      catima::DataPoint b(p,water);
      catima::DataPoint c(p,graphite);
      catima::DataPoint d;
      d = c;
      EXPECT(a == b);
      EXPECT(!(a==c));
      EXPECT(d == c);
      EXPECT(!(d==b));
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
      
      catima::_storage.Reset();
      EXPECT(catima::_storage.get_index()==0);
      
      catima::_storage.Add(p,water);
      EXPECT(catima::_storage.get_index()==1);
      
      catima::_storage.Add(p,water);
      EXPECT(catima::_storage.get_index()==1);
      
      catima::_storage.Add(p,graphite);
      EXPECT(catima::_storage.get_index()==2);
      
      catima::_storage.Add(p,graphite);
      EXPECT(catima::_storage.get_index()==2);
    },
    CASE("storage limit"){
      for(int i=0;i<100;i++){
        auto m = catima::get_material(i);
      }
    },
    CASE("energy table"){
      double step = (catima::logEmax - catima::logEmin)/(catima::max_datapoints-1);
      EXPECT(catima::energy_table.step==step);
      EXPECT(catima::energy_table.values[0]==exp(M_LN10*(catima::logEmin)));
      EXPECT(catima::energy_table.values[1]==exp(M_LN10*(catima::logEmin+step)));
      EXPECT(catima::energy_table.values[catima::max_datapoints-1]==exp(M_LN10*(catima::logEmax)));
    }
};

int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}


