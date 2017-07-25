#include "catima/catima.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
    catima::Material graphite;
    graphite.add_element(12,6,1); // arguments are A,Z, stoichiometric number
    graphite.density(1.8);        // density in g/cm3
    graphite.thickness(2.0);      // thickness in g/cm2
       
    catima::Material water({ // material with 2 atoms
        {1,1,2}, // 1H - two atoms
        {16,8,1} // 16O - 1 atom
    });
    water.density(1.0).thickness(2.0);
    
    catima::Projectile p(12,6); // define projectile, ie 12C

    catima::Layers layer1; // create layers from materials defined above
    layer1.add(graphite);
    layer1.add(water);
    layer1.add(graphite);
    
    auto results = catima::calculate(p(1000),layer1); //calculate rtansport through layers with initial energy 1000 MeV/u


    return 0;
}
