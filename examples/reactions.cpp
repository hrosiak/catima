#include "catima/catima.h"
#include "catima/reactions.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
    catima::Material target = catima::get_material(4);
    target.thickness(1.0);      // thickness in g/cm2
    catima::Projectile p(12,6); // define projectile, ie 12C
    
    double cs = 45;
    double rcsi = 870;
    double rcso = 860;
    
    cout<<"C->Be\n";
    cout<<"t(g/cm2)\t rate"<<endl;
    for(double t=0.25; t<=5;t+=0.25){
        target.thickness(t);
        double r = production_rate(45,rcsi, rcso, target);
        cout<<t<<"\t"<<r<<endl;
    }
    
    return 0;
}
