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
		 case material::CO2:	 return Material({{0,6,1},{0,8,2}},0.001842);
		 case material::CH4:	 return Material({{0,1,4},{0,6,1}},0.0006);
		 case material::Bakelite:	 return Material({{0,1,0.057441},{0,6,0.774591},{0,8,0.167968}},1.25);
		 case material::A-150 platics:	 return Material({{0,1,0.101327},{0,6,0.7755},{0,7,0.035057},{0,8,0.0523159},{0,9,0.017422},{0,20,0.018378}},1.127);
		 case material::B-100 platics:	 return Material({{0,1,0.0654709},{0,6,0.536944},{0,7,0.0215},{0,8,0.032085},{0,9,0.167411},{0,20,0.176589}},1.45);
		 case material::Adenine:	 return Material({{0,1,5},{0,6,5},{0,7,5}},1.35);
		 case material::Ammonia:	 return Material({{0,1,3},{0,7,1}},0.000826);
		 case material::BaF2:	 return Material({{0,9,2},{0,56,1}},4.89);
		 case material::BaSO4:	 return Material({{0,8,4},{0,16,1},{0,56,1}},4.5);
		 case material::BaO:	 return Material({{0,4,1},{0,8,1}},3.01);
		 case material::BGO:	 return Material({{0,8,12},{0,32,3},{0,83,4}},7.13);
		 case material::Blood:	 return Material({{0,1,0.101866},{0,6,0.10002},{0,7,0.02964},{0,8,0.759414},{0,11,0.00185},{0,12,4e-05},{0,14,3e-05},{0,15,0.00035},{0,16,0.00185},{0,17,0.00278},{0,19,0.00163},{0,20,6e-05},{0,26,0.00046},{0,30,1e-05}},1.06);
		 case material::Bone_Compact:	 return Material({{0,1,0.063984},{0,6,0.278},{0,7,0.027},{0,8,0.410016},{0,12,0.002},{0,15,0.07},{0,16,0.002},{0,20,0.147}},1.85);
		 case material::Bone_Cortical:	 return Material({{0,1,0.047234},{0,6,0.14433},{0,7,0.04199},{0,8,0.446096},{0,12,0.0022},{0,15,0.10497},{0,16,0.00315},{0,20,0.20993},{0,30,0.0001}},1.85);
		 case material::Brain_ICRP:	 return Material({{0,1,0.110667},{0,6,0.12542},{0,7,0.01328},{0,8,0.737723},{0,11,0.00184},{0,12,0.00015},{0,15,0.00354},{0,16,0.00177},{0,17,0.00236},{0,19,0.0031},{0,20,9e-05},{0,26,5e-05},{0,30,1e-05}},1.03);
		 case material::CdTe:	 return Material({{0,48,1},{0,52,1}},6.2);
		 case material::CdWO4:	 return Material({{0,8,4},{0,48,1},{0,74,1}},7.9);
		 case material::CaCO3:	 return Material({{0,6,1},{0,8,3},{0,20,1}},2.8);
		 case material::CaF2:	 return Material({{0,9,2},{0,20,1}},3.18);
		 case material::CaWO4:	 return Material({{0,8,4},{0,20,1},{0,74,1}},6.062);
		 case material::CsF:	 return Material({{0,9,1},{0,55,1}},4.115);
		 case material::CsI:	 return Material({{0,53,1},{0,55,1}},4.51);
		 case material::Concrete:	 return Material({{0,1,0.01},{0,6,0.001},{0,8,0.529107},{0,11,0.016},{0,12,0.002},{0,13,0.033872},{0,14,0.337021},{0,19,0.013},{0,20,0.044},{0,26,0.014}},2.3);
		 case material::Eye_lens:	 return Material({{0,1,0.099269},{0,6,0.19371},{0,7,0.05327},{0,8,0.653751}},1.1);
		 case material::Lung:	 return Material({{0,1,0.101278},{0,6,0.10231},{0,7,0.02865},{0,8,0.757072},{0,11,0.00184},{0,12,0.00073},{0,15,0.0008},{0,16,0.00225},{0,17,0.00266},{0,19,0.00194},{0,20,9e-05},{0,26,0.00037},{0,30,1e-05}},1.05);
		 case material::Muscle_skeletal:	 return Material({{0,1,0.100637},{0,6,0.10783},{0,7,0.02768},{0,8,0.754773},{0,11,0.00075},{0,12,0.00019},{0,15,0.0018},{0,16,0.00241},{0,17,0.00079},{0,19,0.00302},{0,20,3e-05},{0,26,4e-05},{0,30,5e-05}},1.04);
		 case material::Muscle_strained:	 return Material({{0,1,0.101997},{0,6,0.123},{0,7,0.035},{0,8,0.729003},{0,11,0.0008},{0,12,0.0002},{0,15,0.002},{0,16,0.005},{0,19,0.003}},1.04);
		 case material::Skin:	 return Material({{0,1,0.100588},{0,6,0.22825},{0,7,0.04642},{0,8,0.619002},{0,11,7e-05},{0,12,6e-05},{0,15,0.00033},{0,16,0.00159},{0,17,0.00267},{0,19,0.00085},{0,20,0.00015},{0,26,1e-05},{0,30,1e-05}},1.1);
            default:break;
             }
        return Material();
    }

}  // namespace catima
