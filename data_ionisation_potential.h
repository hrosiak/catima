/*
 *  Author: Andrej Prochazka
 *  Copyright(C) 2017
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IONISATION_POTENTIALS_H
#define IONISATION_POTENTIALS_H

#define IPOT_ZMAX 120
namespace catima{

const double ionisation_potentials_z[IPOT_ZMAX+1] = {
    9999999,  // 0   
    19.2, // H
    41.8, // He
    40.0,
    63.7,
    76.0,
    78.0, // C
    82.0,
    95.0,
    115.0,
    137.0,
    149.0,
    156.0,
    166.0,
    173.0,
    173.0,
    180.0,
    174.0,
    188.0,
    190.0,
    191.0, //Ca (20)
    216.0,
    233.0,
    245.0,
    257.0,
    272.0,
    286.0, //Fe 26
    297.0,
    311.0,
    322.0,
    330.0,
    334.0,
    350.0,
    347.0,
    348.0,
    343.0,
    352.0,
    363.0,
    366.0,
    379.0,
    393.0,
    417.0,
    424.0,
    428.0,
    441.0,
    449.0,
    470.0,
    470.0,
    469.0,
    488.0,
    488.0, //Sn 50
    487.0,
    485.0,
    491.0,
    482.0, // Xe 54
    488.0,
    491.0,
    501.0,
    523.0,
    535.0,
    546.0,
    560.0,
    574.0,
    580.0,
    591.0,
    614.0,
    628.0,
    650.0,
    658.0,
    674.0,
    684.0,
    694.0,
    705.0,
    718.0,
    727.0,
    736.0,
    746.0,
    757.0,
    790.0,
    790.0, // Au 79
    800.0,
    810.0,
    823.0, //Pb 82
    823.0,
    830.0,
    825.0,
    794.0,
    827.0,
    826.0,
    841.0,
    847.0,
    878.0,
    890.0, // U 92
      900.0,
      910.0,
      920.0,
      930.0,
      940.0,
      950.0,
      960.0,
     970.0,
     980.0,
     990.0,
    1000.0,
    1010.0,
    1020.0,
    1030.0,
    1040.0,
    1050.0,
    1060.0,
    1070.0,
    1080.0,
    1090.0,
    1100.0,
    1110.0,
    1120.0,
    1130.0,
    1140.0,
    1150.0,
    1160.0,
    1170.0,
    };

inline double ipot(int z){
    if(z<1 || z>IPOT_ZMAX){
        return 0.0;
    }
    return ionisation_potentials_z[z];
}

}
#endif
