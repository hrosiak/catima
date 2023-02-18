#include "catima/catima.h"
#include <iostream>

using std::cout;
using std::endl;


int main(){
    catima::Material graphite = catima::get_material(6);
    catima::Material P10 = catima::get_material(catima::material::P10);

    cout<<"Material info"<<endl;
    cout<<"Molar Mass = "<<graphite.M()<<endl;
    cout<<"density = "<<graphite.density()<<endl;

    cout<<"Material info"<<endl;
    cout<<"Molar Mass = "<<P10.M()<<endl;
    cout<<"density = "<<P10.density()<<endl;

    return 0;
}
