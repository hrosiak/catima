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
		 case material::CH2:	 return Material({{0,6,0.856289},{0,1,0.143711}},0.94);
		 case material::LH2:	 return Material({{0,1,1}},0.0708);
		 case material::LD2:	 return Material({{2.01355,1,1}},0.162);
		 case material::Water:	 return Material({{0,1,0.111894},{0,8,0.888106}},1);
		 case material::Diamond:	 return Material({{0,6,1}},3.52);
		 case material::Glass:	 return Material({{0,14,0.37722},{0,8,0.539562},{0,5,0.040064},{0,11,0.028191},{0,13,0.011644},{0,19,0.003321}},2.4);
		 case material::ALMG3:	 return Material({{0,13,1}},2.67);
		 case material::ArCO2_30:	 return Material({{0,18,0.679},{0,8,0.2321},{0,6,0.0889}},0.00171);
		 case material::CF4:	 return Material({{0,6,0.13636},{0,9,0.86363}},0.00372);
		 case material::Isobutane:	 return Material({{0,6,0.82758},{0,1,0.17241}},0.00251);
		 case material::Kapton:	 return Material({{0,1,0.026362},{0,6,0.691133},{0,7,0.07327},{0,8,0.209235}},1.42);
		 case material::Mylar:	 return Material({{0,6,0.625016},{0,1,0.041959},{0,8,0.333025}},1.38);
		 case material::NaF:	 return Material({{0,11,0.54761},{0,9,0.45239}},2.56);
		 case material::P10:	 return Material({{0,18,0.64171},{0,6,0.268247},{0,1,0.09004}},0.00166);
		 case material::PolyOlefin:	 return Material({{0,8,0.90909},{0,1,0.09091}},0.9);
		 case material::CmO2:	 return Material({{0,96,0.8853},{0,8,0.1147}},12);
		 case material::Suprasil:	 return Material({{0,14,0.37722},{0,8,0.539562},{0,5,0.040064},{0,11,0.028191},{0,13,0.011644},{0,19,0.003321}},2.2);
		 case material::HAVAR:	 return Material({{0,27,0.403228},{0,24,0.169412},{0,28,0.124301},{0,74,0.089847},{0,42,0.031259},{0,26,0.181952}},8.3);
		 case material::Steel:	 return Material({{0,26,0.74621},{0,24,0.169},{0,28,0.08479}},8);
		 case material::CH4:	 return Material({{0,1,0.251306},{0,6,0.748694}},0.0006);
            default:break;
             }
        return Material();
    }

}  // namespace catima
