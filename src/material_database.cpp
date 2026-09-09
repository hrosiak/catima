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
		 case material::Air:	 return Material({{0,7,0.755267},{0,8,0.231781},{0,18,0.012827},{0,6,0.000124}},0.001205);
		 case material::CH2:	 return Material({{0,6,1},{0,1,2}},0.94);
		 case material::lH2:	 return Material({{0,1,1}},0.0708);
		 case material::lD2:	 return Material({{2.0141,1,1}},0.162);
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
		 case material::P10:	 return Material({{0,18,0.9},{0,6,0.1},{0,1,0.4}},0.00166);
		 case material::Polyolefin:	 return Material({{0,8,8},{0,1,16}},0.9);
		 case material::CmO2:	 return Material({{0,96,1},{0,8,2}},12);
		 case material::Suprasil:	 return Material({{0,14,0.37722},{0,8,0.539562},{0,5,0.040064},{0,11,0.028191},{0,13,0.011644},{0,19,0.003321}},2.2);
		 case material::HAVAR:	 return Material({{0,27,0.403228},{0,24,0.169412},{0,28,0.124301},{0,74,0.089847},{0,42,0.031259},{0,26,0.181952}},8.3);
		 case material::Steel:	 return Material({{0,26,0.74621},{0,24,0.169},{0,28,0.08479}},8);
		 case material::CO2:	 return Material({{0,6,1},{0,8,2}},0.001842);
		 case material::Methane:	 return Material({{0,1,4},{0,6,1}},0.0006);
		 case material::Methanol:	 return Material({{0,1,4},{0,6,1},{0,8,1}},0.792);
		 case material::Acetone:	 return Material({{0,1,6},{0,6,3},{0,8,1}},0.7899);
		 case material::Acetylene:	 return Material({{0,1,2},{0,6,2}},0.0010967);
		 case material::Adenine:	 return Material({{0,1,5},{0,6,5},{0,7,5}},1.35);
		 case material::Adipose_Tissue:	 return Material({{0,1,0.119477},{0,6,0.63724},{0,7,0.00797},{0,8,0.232333},{0,11,0.0005},{0,12,2e-05},{0,15,0.00016},{0,16,0.00073},{0,17,0.00119},{0,19,0.00032},{0,20,2e-05},{0,26,2e-05},{0,30,2e-05}},0.92);
		 case material::Alanine:	 return Material({{0,1,7},{0,6,3},{0,7,1},{0,8,2}},1.42);
		 case material::Bakelite:	 return Material({{0,1,0.057441},{0,6,0.774591},{0,8,0.167968}},1.25);
		 case material::AgBr:	 return Material({{0,47,1},{0,35,1}},6.473);
		 case material::AgCl:	 return Material({{0,47,1},{0,17,1}},5.56);
		 case material::AgI:	 return Material({{0,47,1},{0,53,1}},5.675);
		 case material::Al2O3:	 return Material({{0,8,3},{0,13,2}},3.97);
		 case material::Amber:	 return Material({{0,1,0.10593},{0,6,0.788974},{0,8,0.105096}},1.1);
		 case material::Ammonia:	 return Material({{0,1,3},{0,7,1}},0.000826);
		 case material::Aniline:	 return Material({{0,1,7},{0,6,6},{0,7,1}},1.0235);
		 case material::Anthracene:	 return Material({{0,1,10},{0,6,14}},1.283);
		 case material::A_150:	 return Material({{0,1,0.101327},{0,6,0.7755},{0,7,0.035057},{0,8,0.0523159},{0,9,0.017422},{0,20,0.018378}},1.127);
		 case material::B_100:	 return Material({{0,1,0.0654709},{0,6,0.536944},{0,7,0.0215},{0,8,0.032085},{0,9,0.167411},{0,20,0.176589}},1.45);
		 case material::BaF2:	 return Material({{0,9,2},{0,56,1}},4.89);
		 case material::BaSO4:	 return Material({{0,8,4},{0,16,1},{0,56,1}},4.5);
		 case material::Benzene:	 return Material({{0,1,6},{0,6,6}},0.87865);
		 case material::BeO:	 return Material({{0,4,1},{0,8,1}},3.01);
		 case material::BGO:	 return Material({{0,8,12},{0,32,3},{0,83,4}},7.13);
		 case material::Blood_ICRP:	 return Material({{0,1,0.101866},{0,6,0.10002},{0,7,0.02964},{0,8,0.759414},{0,11,0.00185},{0,12,4e-05},{0,14,3e-05},{0,15,0.00035},{0,16,0.00185},{0,17,0.00278},{0,19,0.00163},{0,20,6e-05},{0,26,0.00046},{0,30,1e-05}},1.06);
		 case material::Bone_Compact:	 return Material({{0,1,0.063984},{0,6,0.278},{0,7,0.027},{0,8,0.410016},{0,12,0.002},{0,15,0.07},{0,16,0.002},{0,20,0.147}},1.85);
		 case material::Bone_Cortical:	 return Material({{0,1,0.047234},{0,6,0.14433},{0,7,0.04199},{0,8,0.446096},{0,12,0.0022},{0,15,0.10497},{0,16,0.00315},{0,20,0.20993},{0,30,0.0001}},1.85);
		 case material::Brain_ICRP:	 return Material({{0,1,0.110667},{0,6,0.12542},{0,7,0.01328},{0,8,0.737723},{0,11,0.00184},{0,12,0.00015},{0,15,0.00354},{0,16,0.00177},{0,17,0.00236},{0,19,0.0031},{0,20,9e-05},{0,26,5e-05},{0,30,1e-05}},1.03);
		 case material::B4C:	 return Material({{0,5,4},{0,6,1}},2.52);
		 case material::BC_400:	 return Material({{0,1,10},{0,6,9}},1.032);
		 case material::nButanol:	 return Material({{0,1,10},{0,6,4},{0,8,1}},0.81);
		 case material::C_552:	 return Material({{0,1,0.02468},{0,6,0.50161},{0,8,0.004527},{0,9,0.465209},{0,14,0.003973}},1.76);
		 case material::CdTe:	 return Material({{0,48,1},{0,52,1}},6.2);
		 case material::CdWO4:	 return Material({{0,8,4},{0,48,1},{0,74,1}},7.9);
		 case material::CaCO3:	 return Material({{0,6,1},{0,8,3},{0,20,1}},2.8);
		 case material::CaF2:	 return Material({{0,9,2},{0,20,1}},3.18);
		 case material::CaO:	 return Material({{0,8,1},{0,20,1}},3.34);
		 case material::CaWO4:	 return Material({{0,8,4},{0,20,1},{0,74,1}},6.062);
		 case material::CsF:	 return Material({{0,9,1},{0,55,1}},4.115);
		 case material::CsI:	 return Material({{0,53,1},{0,55,1}},4.51);
		 case material::CCl4:	 return Material({{0,6,1},{0,17,4}},1.594);
		 case material::Tetrachloroethylene:	 return Material({{0,6,2},{0,17,4}},1.622);
		 case material::Cellophane:	 return Material({{0,1,0.062162},{0,6,0.444462},{0,8,0.493376}},1.42);
		 case material::Chlorobenzene:	 return Material({{0,1,5},{0,6,6},{0,17,1}},1.1058);
		 case material::Chloroform:	 return Material({{0,1,1},{0,6,1},{0,17,3}},1.4832);
		 case material::Cyclohexane:	 return Material({{0,1,12},{0,6,6}},0.779);
		 case material::Concrete:	 return Material({{0,1,0.01},{0,6,0.001},{0,8,0.529107},{0,11,0.016},{0,12,0.002},{0,13,0.033872},{0,14,0.337021},{0,19,0.013},{0,20,0.044},{0,26,0.014}},2.3);
		 case material::Diethyl_Ether:	 return Material({{0,1,10},{0,6,4},{0,8,1}},0.71378);
		 case material::Ethane:	 return Material({{0,1,6},{0,6,2}},0.00125324);
		 case material::Ethanol:	 return Material({{0,1,6},{0,6,2},{0,8,1}},0.7893);
		 case material::Ethylene:	 return Material({{0,1,4},{0,6,2}},0.00117497);
		 case material::Eye_lens:	 return Material({{0,1,0.099269},{0,6,0.19371},{0,7,0.05327},{0,8,0.653751}},1.1);
		 case material::Fe2O3:	 return Material({{0,8,3},{0,26,2}},5.242);
		 case material::FeO:	 return Material({{0,8,1},{0,26,1}},5.745);
		 case material::Freon_12:	 return Material({{0,6,1},{0,9,2},{0,17,2}},1.486);
		 case material::Freon_12B2:	 return Material({{0,6,1},{0,9,2},{0,35,2}},2.27);
		 case material::Freon_13:	 return Material({{0,6,1},{0,9,3},{0,17,1}},1.526);
		 case material::Freon_13B1:	 return Material({{0,6,1},{0,9,3},{0,35,1}},1.538);
		 case material::Freon_13I1:	 return Material({{0,6,1},{0,9,3},{0,53,1}},1.538);
		 case material::Gd2O2S:	 return Material({{0,8,2},{0,16,1},{0,64,2}},7.44);
		 case material::GaAs:	 return Material({{0,31,1},{0,33,1}},5.3176);
		 case material::Gel_Photo_Emulsion:	 return Material({{0,1,0.08118},{0,6,0.41606},{0,7,0.11124},{0,8,0.38064},{0,16,0.01088}},1.2914);
		 case material::Glass_Pyrex:	 return Material({{0,5,0.0400639},{0,8,0.539561},{0,11,0.0281909},{0,13,0.011644},{0,14,0.377219},{0,19,0.00332099}},2.23);
		 case material::Glass_Lead:	 return Material({{0,8,0.156453},{0,14,0.080866},{0,22,0.008092},{0,33,0.002651},{0,82,0.751938}},6.22);
		 case material::Glucose:	 return Material({{0,1,12},{0,6,6},{0,8,6}},1.54);
		 case material::Glutamine:	 return Material({{0,1,10},{0,6,5},{0,7,2},{0,8,3}},1.46);
		 case material::Glycerol:	 return Material({{0,1,8},{0,6,3},{0,8,3}},1.2613);
		 case material::Guanine:	 return Material({{0,1,5},{0,6,5},{0,7,5},{0,8,1}},1.58);
		 case material::Gypsum:	 return Material({{0,1,0.023416},{0,8,0.557572},{0,16,0.186215},{0,20,0.232797}},2.32);
		 case material::nHeptane:	 return Material({{0,1,16},{0,6,7}},0.68376);
		 case material::nHexane:	 return Material({{0,1,14},{0,6,6}},0.66603);
		 case material::KI:	 return Material({{0,19,1},{0,53,1}},3.13);
		 case material::K2O:	 return Material({{0,8,1},{0,19,2}},2.32);
		 case material::LaBr3:	 return Material({{0,57,1},{0,35,3}},5.06);
		 case material::LaOBr:	 return Material({{0,8,1},{0,35,1},{0,57,1}},6.28);
		 case material::La2O2S:	 return Material({{0,8,2},{0,16,1},{0,57,2}},5.86);
		 case material::Lung:	 return Material({{0,1,0.101278},{0,6,0.10231},{0,7,0.02865},{0,8,0.757072},{0,11,0.00184},{0,12,0.00073},{0,15,0.0008},{0,16,0.00225},{0,17,0.00266},{0,19,0.00194},{0,20,9e-05},{0,26,0.00037},{0,30,1e-05}},1.05);
		 case material::MgCO3:	 return Material({{0,12,1},{0,6,1},{0,8,3}},2.958);
		 case material::MgF2:	 return Material({{0,12,1},{0,9,2}},3.148);
		 case material::MgO:	 return Material({{0,12,1},{0,8,1}},3.6);
		 case material::MS20_Tissue:	 return Material({{0,1,0.081192},{0,6,0.583442},{0,7,0.017798},{0,8,0.186381},{0,12,0.130287},{0,17,0.0009}},1);
		 case material::Muscle_skeletal:	 return Material({{0,1,0.100637},{0,6,0.10783},{0,7,0.02768},{0,8,0.754773},{0,11,0.00075},{0,12,0.00019},{0,15,0.0018},{0,16,0.00241},{0,17,0.00079},{0,19,0.00302},{0,20,3e-05},{0,26,4e-05},{0,30,5e-05}},1.04);
		 case material::Muscle_strained:	 return Material({{0,1,0.101997},{0,6,0.123},{0,7,0.035},{0,8,0.729003},{0,11,0.0008},{0,12,0.0002},{0,15,0.002},{0,16,0.005},{0,19,0.003}},1.04);
		 case material::Muscle_sucrose:	 return Material({{0,1,0.0982341},{0,6,0.156214},{0,7,0.035451},{0,8,0.710101}},1.11);
		 case material::Muscle_no_sucrose:	 return Material({{0,1,0.101969},{0,6,0.120058},{0,7,0.035451},{0,8,0.742522}},1.07);
		 case material::Na2CO3:	 return Material({{0,6,1},{0,8,3},{0,11,1}},2.532);
		 case material::NaI:	 return Material({{0,11,1},{0,53,1}},3.667);
		 case material::NaCl:	 return Material({{0,11,1},{0,17,1}},2.165);
		 case material::Na2O:	 return Material({{0,8,1},{0,11,2}},2.27);
		 case material::NaNO3:	 return Material({{0,7,1},{0,8,3},{0,11,1}},2.261);
		 case material::Naphthalene:	 return Material({{0,1,8},{0,6,10}},1.145);
		 case material::Nitrobenzene:	 return Material({{0,1,5},{0,6,6},{0,7,1},{0,8,2}},1.199);
		 case material::N2O:	 return Material({{0,7,2},{0,8,1}},0.00183);
		 case material::Octane:	 return Material({{0,1,18},{0,6,8}},0.703);
		 case material::Paraffin:	 return Material({{0,1,0.148605},{0,6,0.851395}},0.93);
		 case material::nPentane:	 return Material({{0,1,12},{0,6,5}},0.626);
		 case material::PhotoEmulsion:	 return Material({{0,1,0.0141},{0,6,0.072261},{0,7,0.01932},{0,8,0.066101},{0,16,0.00189},{0,35,0.349103},{0,47,0.474105}},3.815);
		 case material::PuO2:	 return Material({{0,8,2},{0,94,1}},11.46);
		 case material::Polyacrylonitrile:	 return Material({{0,1,3},{0,6,3},{0,7,1}},1.184);
		 case material::Polycarbonate:	 return Material({{0,1,0.055491},{0,6,0.755751},{0,8,0.188758}},1.2);
		 case material::PMMA:	 return Material({{0,1,8},{0,6,5},{0,8,2}},1.18);
		 case material::POM:	 return Material({{0,1,2},{0,6,1},{0,8,1}},1.42);
		 case material::Polypropylene:	 return Material({{0,6,3},{0,1,6}},0.9);
		 case material::Polystyrene:	 return Material({{0,6,8},{0,1,8}},1.06);
		 case material::Propane:	 return Material({{0,1,8},{0,6,3}},0.00188);
		 case material::nPropanol:	 return Material({{0,1,8},{0,6,3},{0,8,1}},0.8035);
		 case material::PVC:	 return Material({{0,1,3},{0,6,2},{0,17,1}},1.3);
		 case material::Pyridine:	 return Material({{0,1,5},{0,6,5},{0,7,1}},0.9819);
		 case material::SiO2:	 return Material({{0,8,2},{0,14,1}},2.32);
		 case material::Skin:	 return Material({{0,1,0.100588},{0,6,0.22825},{0,7,0.04642},{0,8,0.619002},{0,11,7e-05},{0,12,6e-05},{0,15,0.00033},{0,16,0.00159},{0,17,0.00267},{0,19,0.00085},{0,20,0.00015},{0,26,1e-05},{0,30,1e-05}},1.1);
		 case material::Sucrose:	 return Material({{0,1,22},{0,6,12},{0,8,11}},1.587);
		 case material::Teflon:	 return Material({{0,6,2},{0,9,4}},2.2);
		 case material::TlCl:	 return Material({{0,17,1},{0,81,1}},7.004);
		 case material::Toluene:	 return Material({{0,1,8},{0,6,7}},0.8669);
		 case material::Trichloroethylene:	 return Material({{0,1,1},{0,6,2},{0,17,3}},1.46);
		 case material::WF6:	 return Material({{0,9,6},{0,74,1}},2.4);
		 case material::UC2:	 return Material({{0,6,2},{0,92,1}},11.28);
		 case material::UC:	 return Material({{0,6,1},{0,92,1}},13.63);
		 case material::UO2:	 return Material({{0,8,2},{0,92,1}},10.97);
		 case material::Urea:	 return Material({{0,1,0.067131},{0,6,0.2},{0,7,0.466459},{0,8,0.266411}},1.323);
		 case material::Valine:	 return Material({{0,1,11},{0,6,5},{0,7,1},{0,8,2}},1.23);
		 case material::Iodonaphthalene:	 return Material({{0,6,10},{0,1,7},{0,53,1}},1.738);
		 case material::C21H24O4:	 return Material({{0,6,21},{0,1,24},{0,8,4}},1.18);
		 case material::CoRe_Alloy:	 return Material({{0,27,17},{0,75,23},{0,24,1}},11.5);
		 case material::LLZO_electrolyte:	 return Material({{0,3,7},{0,57,3},{0,40,2},{0,8,12}},5.1);
		 case material::Nylon:	 return Material({{0,6,6},{0,1,11},{0,7,1},{0,8,1}},1.14);
		 case material::Brass:   return Material({{0,29,0.63},{0,30,0.35},{0,82,0.02}},8.73);
	         case material::RW3:  return Material({{0,1,0.0759},{0,6,0.9041},{0,6,0.008},{0,22,0.012}},1.045);
            default:break;
             }
        return Material();
    }

}  // namespace catima
