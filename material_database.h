#ifndef MATERIAL_DATABASE
#define MATERIAL_DATABASE
#include "catima/catima.h"

namespace catima{

    //namespace material{
        enum material{
            PLASTIC = 201,
            AIR = 202,
            CH2,
            LH2,
            LD2,
            WATER,
            DIAMOND,
            GLASS,
            ALMG3,
            ARCO2_30,
            CF4,
            ISOBUTANE,
            KAPTON,
            MYLAR,
            NAF,
            P10,
            POLYOLEFIN,
            CMO2,
            SUPRASIL,
            HAVAR
        };
    //}

    Material get_material(int);
}

#endif