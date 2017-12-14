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
		 Polyolefin = 217,
		 CmO2 = 218,
		 Suprasil = 219,
		 HAVAR = 220,
		 Steel = 221,
		 CO2 = 222,
		 CH4 = 223,
		 Bakelite = 224,
		 A-150 platics = 225,
		 B-100 platics = 226,
		 Adenine = 227,
		 Ammonia = 228,
		 BaF2 = 229,
		 BaSO4 = 230,
		 BaO = 231,
		 BGO = 232,
		 Blood = 233,
		 Bone_Compact = 234,
		 Bone_Cortical = 235,
		 Brain_ICRP = 236,
		 CdTe = 237,
		 CdWO4 = 238,
		 CaCO3 = 239,
		 CaF2 = 240,
		 CaWO4 = 241,
		 CsF = 242,
		 CsI = 243,
		 Concrete = 244,
		 Eye_lens = 245,
		 Lung = 246,
		 Muscle_skeletal = 247,
		 Muscle_strained = 248,
		 Skin = 249
        };

      Material get_compound(material m);
      Material get_material(int id);
      inline Material get_material(material m){
            return get_compound(m);
      };

}

#endif
