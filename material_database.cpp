#include "catima/material_database.h"
#include "catima/nucdata.h"

namespace catima{

    Material get_material(int id){
        if(id>0 && id<ELEMENT_DENSITY_MAXZ){
            return Material(0,id,element_density(id),0.0);
        }
        return get_compound(static_cast<material>(id));
        }

    Material get_compound(material m){
      switch(m){
		 case material::Plastics:	 return Material({{0,1,10},{0,6,9}},1.032);
		 case material::Air:	 return Material({{0,7,0.781},{0,8,0.2095},{0,18,0.0095}},0.0012);
		 case material::CH2:	 return Material({{0,6,1},{0,1,2}},0.94);
		 case material::LH2:	 return Material({{0,1,1}},0.0708);
		 case material::LD2:	 return Material({{2.01355,1,1}},0.162);
		 case material::Water:	 return Material({{0,1,2},{0,8,1}},1);
		 case material::Diamond:	 return Material({{0,6,1}},3.52);
		 case material::Glass:	 return Material({{0,14,1},{0,8,2}},2.4);
		 case material::ALMG3:	 return Material({{0,13,1}},2.67);
		 case material::ArCO2_30:	 return Material({{0,18,7},{0,8,6},{0,6,3}},0.00171);
		 case material::CF4:	 return Material({{0,6,1},{0,9,4}},0.00372);
		 case material::Isobutane:	 return Material({{0,6,4},{0,1,10}},0.00251);
		 case material::Kapton:	 return Material({{0,1,10},{0,6,22},{0,7,2},{0,8,5}},1.42);
		 case material::Mylar:	 return Material({{0,6,5},{0,1,4},{0,8,2}},1.38);
		 case material::NaF:	 return Material({{0,11,1},{0,9,1}},2.56);
		 case material::P10:	 return Material({{0,18,9},{0,6,1},{0,1,4}},0.00166);
		 case material::PolyOlefin:	 return Material({{0,8,10},{0,1,16}},0.9);
		 case material::CmO2:	 return Material({{0,96,1},{0,8,2}},12);
		 case material::Suprasil:	 return Material({{0,14,1},{0,8,2}},2.2);
		 case material::HAVAR:	 return Material({{0,27,42},{0,24,20},{0,28,13},{0,74,3},{0,42,2},{0,26,20}},8.3);
		 case material::Steel:	 return Material({{0,26,74},{0,24,18},{0,28,8}},8);
            default:break;
             }
        return Material();
    }

}  // namespace catima
