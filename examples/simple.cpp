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
    
    cout<<"C->C\n";
    for(double T:{20,50,100,600,1000}){
        auto result = catima::calculate(p(T),graphite);
	    cout<<"T "<<T<<", dEdx = "<<result.dEdxi<<" MeV/g/cm2"<<", range = "<<result.range<<" g/cm2"<<endl;
	}
    
    cout<<"C->H2O\n";
    for(double T:{20,50,100,600,1000}){
        auto result = catima::calculate(p,water,T);
	    cout<<"T "<<T<<", dEdx = "<<result.dEdxi<<" MeV/g/cm2"<<", range = "<<result.range<<" g/cm2"<<endl;
	}

    return 0;
}
