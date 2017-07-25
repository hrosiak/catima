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

#ifndef nucdata_h
#define nucdata_h

#define ATOMIC_WEIGHT_MAXZ 110
#define ELEMENT_DENSITY_MAXZ 99

namespace catima{
extern double element_atomic_weights[110];
extern double element_densities[99];

inline double element_atomic_weight(int z){
    return (z>0 && z<110)?element_atomic_weights[z]:0.0;
    }

inline double element_density(int z){
    return (z>0 && z<99)?element_densities[z]:0.0;
    }

}



#endif
