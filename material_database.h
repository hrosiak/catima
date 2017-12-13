#ifndef MATERIAL_DATABASE
#define MATERIAL_DATABASE
#include "catima/catima.h"

namespace catima{

        enum class material{
		 Plastics = 201,
		 Air = 202,
		 CH2 = 203,
		 LH2 = 204,
		 LD2 = 205,
		 Water = 206,
		 Diamond = 207,
		 Glass = 208,
		 ALMG3 = 209,
		 ArCO2_30 = 210,
		 CF4 = 211,
		 Isobutane = 212,
		 Kapton = 213,
		 Mylar = 214,
		 NaF = 215,
		 P10 = 216,
		 PolyOlefin = 217,
		 CmO2 = 218,
		 Suprasil = 219,
		 HAVAR = 220,
		 Steel = 221,
		 CH4 = 222
        };

      Material get_compound(material m);
      Material get_material(int id);
      inline Material get_material(material m){
            return get_compound(m);
      };

}

#endif
