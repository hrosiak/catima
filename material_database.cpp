#include "catima/material_database.h"
#include "catima/nucdata.h"

namespace catima{

    Material get_material(int id){
        if(id>0 && id<ELEMENT_DENSITY_MAXZ){
            return Material(0,id,element_density(id),0.0);
        }
        switch(id){
            case material::PLASTIC : return Material({{0,1,10},{0,6,9}},1.032);
            case material::AIR :     return Material({{0,7,0.7810}, {0,8,0.2095},{0,18,0.0095}},0.0012);
            case material::CH2 :     return Material({{0,6,1}, {0,1,2}},0.94);
            case material::LH2 :     return Material({{0,1,1}},0.0708);
            case material::LD2 :     return Material({{2.014,1,1}},0.162);
            case material::WATER :   return Material({{0,1,2},{0,8,1}},1.0);
            case material::DIAMOND:  return Material(0,6,3.52,0.0);
            case material::GLASS :   return Material({{0,14,1},{0,8,2}},2.4 );
            case material::ALMG3:    return Material(0,13,2.67,0.0);
            case material::ARCO2_30: return Material({{0,18,7},{0,8,6},{0,6,3}}, 0.00171);
            case material::CF4 :     return Material({{0,6,1},{0,9,4}}, 0.00372 );
            case material::ISOBUTANE:return Material({{0,6,4}, {0,1,10}},0.00251);
            case material::KAPTON :  return Material({{0,1,10}, {0,6,22},{0,7,2},{0,8,5}},1.42);
            case material::MYLAR :   return Material({{0,6,5}, {0,1,4},{0,8,2}},1.38);
            case material::NAF :     return Material({{0,11,1}, {0,9,1}},2.56);
            case material::P10:      return Material({{0,18,9},{0,6,1},{0,1,4}}, 0.00166);
            case material::POLYOLEFIN:  return Material({{0,1,16},{0,6,10}}, 0.9);
            case material::CMO2:      return Material({{0,96,1},{0,8,2}}, 12.0);
            case material::SUPRASIL :   return Material({{0,14,1},{0,8,2}},2.2 );
            case material::HAVAR :   return Material({{0,27,42},{0,24,40},{0,28,13},{0,26,19},{0,74,1}},8.3);
            case material::STEEL :   return Material({{0,26,74},{0,24,18},{0,28,8}},8.0);
            default:break;
             }
        return Material();
        }

}
