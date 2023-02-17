/*
 *  Author: Andrej Prochazka
 *  Copyright(C) 2017
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONVERT_H
#define CONVERT_H
#include "catima/structures.h"
#include "catima/material_database.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "fmt/format.h"

using namespace catima;

bool mocadi_material_match(const Material &a, const Material&b){
    if(std::fabs(a.density() - b.density())> 1e-6)return false;
    if(a.ncomponents() != b.ncomponents())return false;
    for(int i=0;i<a.ncomponents();i++){
        if(a.get_element(i).stn != b.get_element(i).stn)return false;
        if(a.get_element(i).A != b.get_element(i).A)return false;
        if(a.get_element(i).Z != b.get_element(i).Z)return false;
    }
    if(a.M() != b.M())return false;
    return true;
}

int materialdb_id(const Material &m){
    for(int i=201;i<=347;i++){
        if(mocadi_material_match(get_material(i),m)){
            return i;
        }
    }
    return 0;
}

bool save_mocadi(const char* filename, const Projectile p, const Layers &layers, const Phasespace &psx={}, const Phasespace &psy={}){
    std::ofstream fw;
    fw.open(filename, std::ios::out);
    if (!fw.is_open()) { return false;}
    fw<<"epax 2\natima-1.0\noption listmode root\n";
    std::string beam = fmt::format("BEAM\n100000\n{}, 0, {}, {}\n2\n{}, {}, 0, 0, 0\n2\n{}, {}, 0, 0, 0\n1\n0,  0, 0, 0, 0\n",p.T,p.A,p.Z,psx.sigma_x,1000*psx.sigma_a, psy.sigma_x,1000*psy.sigma_a);
    fw<<beam;
    int c = 0;
    for (auto& m: layers.get_materials()){
        int z = 0;
        double a = 0.0;
        if(m.ncomponents()==1){
            z = m.get_element(0).Z;
            a = m.get_element(0).A;
            }
        else{
            z = materialdb_id(m);
        }        
        std::string mstr = fmt::format("*******\nMATTER\n{}, {}, {}\n2,{}\n0.\n0.,0.,0.\n1,1,0\n0,0,0,0\n",a,z,m.density()*1000,m.thickness_cm());
        fw<<mstr;        
        c++;
    }
    fw<<"ERWARTUNGSWERTE\nSAVE\nEND";
    fw.close();
    return true;
}

#endif