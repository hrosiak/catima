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

#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include <array>
#include <initializer_list>

namespace catima{

    /**
     * Projectile class
     * Example usage:
     * \code{.cpp}
     * Projectile p(12,6);
     * p.T = 1000; // setting energy to 1000 MeV/u
     * p(1000);    // setting energy to 1000 MeV/u
     * \endcode
     */
    struct Projectile{
        double A=0;
        double Z=0;
        double Q=0;
        double T=0;
        Projectile& operator()(double e){T=e;return *this;}
        Projectile(){}
        Projectile(double a, double z, double q=0, double t=0):A(a),Z(z),Q(q),T(t){if(q==0)Q=Z;}
    };

    bool operator==(const Projectile &a, const Projectile&b);
    
    /**
      * Target class is used to store constituents of the Material class
      * its just to hold A,Z data for individual atoms.
      */
    struct Target{
        double A=0;
        int Z=0;
        double stn=1.0;
    };

    /**
     * Material
     * class to store Material information
     * This class defines the material with which projectile will interact. 
     * The class store nuclei constituents, density, thickness etc.
     */
    class Material{
        private:
            double rho=0;
            double th=0;
            double molar_mass=0;
            double i_potential=0;
            std::vector<Target>atoms;

        public:
            Material(){};
            /**
              * constructor to add 1 element into the Material, stn number of the element is set to 1.0
              * @param _a - mass number of the atom, is 0 the atomic weight of element _z is taken
              * @param _z - proton number of the atom
              * @param _rho - density of the material in g/cm3, default 0.0
              * @param _th - thickness of the material in g/cm2, default 0.0
              */
            Material(double _a, int _z, double _rho=0.0, double _th=0.0);


            /**
              * constructor to add 1 or multiple element into the Material
              * \code{.cpp}
              * Maetrial water({
                  {1,1,2},
                  {16,8,1},
                });
              * \endcode
              */
            Material(std::initializer_list<std::array<double,3>>list,double _density=0.0, double ipot = 0.0);
            
            /**
             * calculates internal variables if needed
             */
            void calculate();
            
            /**
             * add atom with mass number _a and proton number _z to the Material 
             * @param _a - mass number of the atom, is 0 the atomic weight of element _z is taken
             * @param _z - proton number of the atom
             * @param _stn - stoichiomatric number
             */
            void add_element(double _a, int _z, double _stn);
            
            /**
             * returns i-th element of the Material as a std::pair of Target and corresponding stoichiometric number
             * @param i - index of element to return
             * @return Target class
             */
            Target get_element(int i) const {return atoms[i];};

            /**
             * return weight fraction of i-th element
             * @return weight fraction
             */
            double weight_fraction(int i) const {return (atoms[i].stn<1.0)?atoms[i].stn:atoms[i].stn*atoms[i].A/M();};

            /**
              * @return number of components in Material
              */
            int ncomponents() const {return atoms.size();}

            /**
              * @return returns Molar Mass of the Material
              */
            double M() const {return molar_mass;}
            
            /**
              * @return returns density in g/cm^3
              */
            double density() const {return rho;};

            /**
              * sets density in g/cm^3
              */
            Material& density(double val){rho = val;return *this;};
            
            /**
              * @return returns thickness in g/cm^2
              */
            double thickness() const {return th;};

            /**
              * sets thickness in g/cm^2
              */
            Material& thickness(double val){th = val;return *this;};

            /**
              * set the mean ionization potential, if non elemental I should be used
              */
            Material& I(double val){i_potential = val;return *this;};
            
            /**
              * 0 if default elemental potential is used
              * @return returns ionisation potential in ev
              */
            double I() const {return i_potential;};


            friend bool operator==(const Material &a, const Material&b);
    };

    bool operator==(const Material &a, const Material&b);

    /**
      * structure to store results for calculation per Material
      */
    struct Result{
        double Ein=0.0;
        double Eout=0.0;
        double Eloss = 0.0;
        double range=0.0;
        double dEdxi=0.0;
        double dEdxo=0.0;
        double sigma_E=0.0;
        double sigma_a=0.0;
        double sigma_r=0.0;
        double tof=0.0;   
    };

    /**
      * structure to store results for calculation for multiple layers of materials, ie in catima::Layers 
      */
    struct MultiResult{
        std::vector<Result> results;
        Result total_result;
    };


    /**
     * Layers
     * class to store multiple material layers
     */

    class Layers{        
        private:
        std::vector<Material> materials;
        public:
        Layers(){};
        Layers& operator=(const Layers& other);

        /**
          * @return reference to the std::vector of stored Materials
          */
        const std::vector<Material>& get_materials()const{return materials;}

        /**
          * add Material m to the Layers
          * @param m Material 
          */
        void add(Material m);

        /**
          * @return number of stored Materials
          */
        int num()const{return materials.size();};
        
        Material& operator[](int i){return materials[i];}

        friend Layers operator+(const Layers &a, const Layers&b);
        
        friend Layers operator+(const Layers &a, const Material&b);
    };
}
#endif
