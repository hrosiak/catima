/**
 * this example program print out Lindhard - Soerrensen coefficients
 */

#include "catima/catima.h"
#include "catima/storage.h"
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;


int main(int argc, char** argv){
    catima::Projectile p(238.0,92); // define projectile
    if(argc>2){
	double a = atof(argv[1]);
	int z = atoi(argv[2]);
	if(a>0 && z>0 && z<120){
	    p.A = a;
	    p.Z = z;
	}
    }
    cout<<"projectile: A = "<<p.A<<", Z = "<<p.Z<<"\n";
    auto energies = catima::EnergyTable<100>(1,5); // get energy table, energies log distributed between 10^1 and 10000^5;
    for(double T:energies){
        auto ls = catima::bethek_lindhard(p(T));
	auto lsX = catima::bethek_lindhard_X(p(T));
	    cout<<"T "<<T<<", Delta LS = "<<ls<<",  X = "<<lsX<<endl;
	}

    return 0;
}
