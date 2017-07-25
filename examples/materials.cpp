#include "catima/catima.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
    catima::Material graphite = catima::get_material(6);    
    catima::Material water = catima::get_material(catima::material::WATER);

    cout<<"Material info"<<endl;
    cout<<"Molar Mass = "<<graphite.M()<<endl;
    cout<<"density = "<<graphite.density()<<endl;

    cout<<"Material info"<<endl;
    cout<<"Molar Mass = "<<water.M()<<endl;
    cout<<"density = "<<water.density()<<endl;

    return 0;
}
