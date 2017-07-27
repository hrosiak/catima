/**
 * this example program print out Lindhard - Soerrensen coefficients for 238U
 */

#include "catima/catima.h"
#include "catima/storage.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
       
    catima::Projectile p(238,92); // define projectile, ie 12C
    
    cout<<"projectile 238U\n";

    auto energies = catima::EnergyTable<50>(2,5); // get energy table, energies log distributed between 10^2 and 10000^5;
    for(double T:energies){
        auto ls = catima::bethek_lindhard(p(T));
	auto lsX = catima::bethek_lindhard_X(p(T));
	    cout<<"T "<<T<<", Delta LS = "<<ls<<",  X = "<<lsX<<endl;
	}

    return 0;
}
