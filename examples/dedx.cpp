#include "catima/catima.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){       
    catima::Material water({ // material with 2 atoms
        {1,1,2}, // 1H - two atoms
        {16,8,1} // 16O - 1 atom
    });
    water.density(1.0).thickness(2.0);
    
    catima::Material water2({ // material with 2 atoms
        {1,1,2}, // 1H - two atoms
        {16,8,1} // 16O - 1 atom
    },1.0,78);
    
    water2.thickness(2.0);
    
    catima::Projectile p(12,6); // define projectile, ie 12C
        
    cout<<"C->H2O\n";
    for(double T=0; T<11000;T+=50){
        auto result = catima::calculate(p,water,T/12);
        auto result2 = catima::calculate(p,water2,T/12);
	    cout<<"T = "<<T<<" MeV, dEdx1 = "<<result.dEdxi<<", dEdx2 = "<<result2.dEdxi<<" MeV/g/cm2"<<endl;
	}

    return 0;
}
