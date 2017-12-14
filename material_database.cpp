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
		 case material::Plastics:	 return Material({{0,1,0.085},{0,6,0.915}},1.032);
		 case material::Air:	 return Material({{0,7,0.755267},{0,8,0.231781},{0,18,0.012827},{0,6,0.000124}},0.001205);
		 case material::CH2:	 return Material({{0,6,1},{0,1,2}},0.94);
		 case material::LH2:	 return Material({{0,1,1}},0.0708);
		 case material::LD2:	 return Material({{2.0141,1,1}},0.162);
		 case material::Water:	 return Material({{0,1,2},{0,8,1}},1);
		 case material::Diamond:	 return Material({{0,6,1}},3.52);
		 case material::Glass:	 return Material({{0,14,0.37722},{0,8,0.539562},{0,5,0.040064},{0,11,0.028191},{0,13,0.011644},{0,19,0.003321}},2.4);
		 case material::ALMG3:	 return Material({{0,13,0.97},{0,12,0.03}},2.67);
		 case material::ArCO2_30:	 return Material({{0,18,0.679},{0,8,0.2321},{0,6,0.0889}},0.00171);
		 case material::CF4:	 return Material({{0,6,1},{0,9,4}},0.00372);
		 case material::Isobutane:	 return Material({{0,6,4},{0,1,10}},0.00251);
		 case material::Kapton:	 return Material({{0,1,0.026362},{0,6,0.691133},{0,7,0.07327},{0,8,0.209235}},1.42);
		 case material::Mylar:	 return Material({{0,6,10},{0,1,8},{0,8,4}},1.38);
		 case material::NaF:	 return Material({{0,11,1},{0,9,1}},2.56);
		 case material::P10:	 return Material({{0,18,0.64171},{0,6,0.268247},{0,1,0.09004}},0.00166);
		 case material::Polyolefin:	 return Material({{0,8,8},{0,1,16}},0.9);
		 case material::CmO2:	 return Material({{0,96,1},{0,8,2}},12);
		 case material::Suprasil:	 return Material({{0,14,0.37722},{0,8,0.539562},{0,5,0.040064},{0,11,0.028191},{0,13,0.011644},{0,19,0.003321}},2.2);
		 case material::HAVAR:	 return Material({{0,27,0.403228},{0,24,0.169412},{0,28,0.124301},{0,74,0.089847},{0,42,0.031259},{0,26,0.181952}},8.3);
		 case material::Steel:	 return Material({{0,26,0.74621},{0,24,0.169},{0,28,0.08479}},8);
		 case material::CH4:	 return Material({{0,1,4},{0,6,1}},0.0006);
            default:break;
             }
        return Material();
    }

}  // namespace catima
