#include "structures.h"
#include "catima/nucdata.h"


namespace catima{

bool operator==(const Projectile &a, const Projectile&b){
    if( (a.A==b.A) && (a.Z==b.Z) && (a.Q==b.Q)){
        return true;
        }
    else
        return false;
}

bool operator==(const Material &a, const Material&b){
    if(a.molar_mass != b.molar_mass)return false;
    if(a.density() != b.density())return false;
    if(a.ncomponents() != b.ncomponents())return false;
    for(int i=0;i<a.ncomponents();i++){
        if(a.stn[i] != b.stn[i])return false;
        if(a.atoms[i].A != b.atoms[i].A)return false;
        if(a.atoms[i].Z != b.atoms[i].Z)return false;
    }
    return true;
}

/*
Material::Material(const std::array<double,2> &list){
    add_element(list[0],list[1],1.0);
}
*/
Material::Material(std::initializer_list<std::array<double,3>>list,double _density):rho(_density){
    std::initializer_list<std::array<double,3>>::iterator it;
    for ( it=list.begin(); it!=list.end(); ++it){
        add_element( (*it)[0],(*it)[1],(*it)[2]);
    }
}

Material::Material(double _a, int _z, double _rho, double _th):rho(_rho),th(_th){
    add_element(_a,_z,1.0);
}

void Material::add_element(double _a, int _z, double _stn){
    double a = (_a>0)?_a:element_atomic_weight(_z);
    stn.push_back(_stn);
    atoms.push_back({a,_z});
    molar_mass += _stn*a;
}


std::pair<Target,double> Material::get_element(int i) const{
    return std::pair<Target,double>(atoms[i],stn[i]);
}

/*
void Material::calculate(){
    molar_mass = 0;
    for(int i=0;i<ncomponents();i++){
        molar_mass += stn[i]*a[i];
    }
}
*/

Layers& Layers::operator=(const Layers& other){
    
    materials.clear();
    for(auto&e : other.get_materials()){
        materials.push_back(e);
    }
    return *this;
}

void Layers::add(Material m){
    materials.push_back(m);
}

Layers operator+(const Layers &a, const Layers&b){
    Layers res;
    for(auto &e:a.materials){
        res.add(e);
    }
    
    for(auto &e:b.materials){
        res.add(e);
    }
    return res;
}

Layers operator+(const Layers &a, const Material &m){
    Layers res;
    for(auto &e:a.materials){
        res.add(e);
    }
    res.add(m);
    return res;
}


}
