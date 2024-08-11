#ifndef MATERIAL_DATABASE
#define MATERIAL_DATABASE
#include "catima/catima.h"

namespace catima{

        enum class material{
		 Plastics = 201,
		 Air = 202,
		 CH2 = 203,
		 lH2 = 204,
		 lD2 = 205,
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
		 Methane = 223,
		 Methanol = 224,
		 Acetone = 225,
		 Acetylene = 226,
		 Adenine = 227,
		 Adipose_Tissue = 228,
		 Alanine = 229,
		 Bakelite = 230,
		 AgBr = 231,
		 AgCl = 232,
		 AgI = 233,
		 Al2O3 = 234,
		 Amber = 235,
		 Ammonia = 236,
		 Aniline = 237,
		 Anthracene = 238,
		 A_150 = 239,
		 B_100 = 240,
		 BaF2 = 241,
		 BaSO4 = 242,
		 Benzene = 243,
		 BeO = 244,
		 BGO = 245,
		 Blood_ICRP = 246,
		 Bone_Compact = 247,
		 Bone_Cortical = 248,
		 Brain_ICRP = 249,
		 B4C = 250,
		 BC_400 = 251,
		 nButanol = 252,
		 C_552 = 253,
		 CdTe = 254,
		 CdWO4 = 255,
		 CaCO3 = 256,
		 CaF2 = 257,
		 CaO = 258,
		 CaWO4 = 259,
		 CsF = 260,
		 CsI = 261,
		 CCl4 = 262,
		 Tetrachloroethylene = 263,
		 Cellophane = 264,
		 Chlorobenzene = 265,
		 Chloroform = 266,
		 Cyclohexane = 267,
		 Concrete = 268,
		 Diethyl_Ether = 269,
		 Ethane = 270,
		 Ethanol = 271,
		 Ethylene = 272,
		 Eye_lens = 273,
		 Fe2O3 = 274,
		 FeO = 275,
		 Freon_12 = 276,
		 Freon_12B2 = 277,
		 Freon_13 = 278,
		 Freon_13B1 = 279,
		 Freon_13I1 = 280,
		 Gd2O2S = 281,
		 GaAs = 282,
		 Gel_Photo_Emulsion = 283,
		 Glass_Pyrex = 284,
		 Glass_Lead = 285,
		 Glucose = 286,
		 Glutamine = 287,
		 Glycerol = 288,
		 Guanine = 289,
		 Gypsum = 290,
		 nHeptane = 291,
		 nHexane = 292,
		 KI = 293,
		 K2O = 294,
		 LaBr3 = 295,
		 LaOBr = 296,
		 La2O2S = 297,
		 Lung = 298,
		 MgCO3 = 299,
		 MgF2 = 300,
		 MgO = 301,
		 MS20_Tissue = 302,
		 Muscle_skeletal = 303,
		 Muscle_strained = 304,
		 Muscle_sucrose = 305,
		 Muscle_no_sucrose = 306,
		 Na2CO3 = 307,
		 NaI = 308,
		 NaCl = 309,
		 Na2O = 310,
		 NaNO3 = 311,
		 Naphthalene = 312,
		 Nitrobenzene = 313,
		 N2O = 314,
		 Octane = 315,
		 Paraffin = 316,
		 nPentane = 317,
		 PhotoEmulsion = 318,
		 PuO2 = 319,
		 Polyacrylonitrile = 320,
		 Polycarbonate = 321,
		 PMMA = 322,
		 POM = 323,
		 Polypropylene = 324,
		 Polystyrene = 325,
		 Propane = 326,
		 nPropanol = 327,
		 PVC = 328,
		 Pyridine = 329,
		 SiO2 = 330,
		 Skin = 331,
		 Sucrose = 332,
		 Teflon = 333,
		 TlCl = 334,
		 Toluene = 335,
		 Trichloroethylene = 336,
		 WF6 = 337,
		 UC2 = 338,
		 UC = 339,
		 UO2 = 340,
		 Urea = 341,
		 Valine = 342,
		 Iodonaphthalene = 343,
		 C21H24O4 = 344,
		 CoRe_Alloy = 345,
		 LLZO_electrolyte = 346,
		 Nylon = 347,
		 Brass = 348
        };

      Material get_compound(material m);
      Material get_material(int id);
      inline Material get_material(material m){
            return get_compound(m);
      };

}

#endif
