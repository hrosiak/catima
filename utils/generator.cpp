#include <iostream>
#include "catima/catima.h"
#include "catima/nucdata.h"
#include <functional>
#include <fstream>
#include <sstream>

using namespace catima;
using std::cout;
using std::endl;

//// configuration ////
constexpr double logEmin_forLS = 0.0;
constexpr double logEmax_forLS = catima::logEmax;
constexpr double num_LS_datapoints = 200;
/////////////////////////////////////////

inline double ls_energy_table( int i ) { return exp(M_LN10*(logEmin_forLS + ((double)i)*(logEmax_forLS-logEmin_forLS)/(num_LS_datapoints - 1.0))); }


int main(){
    double r;
    double r2;
    double dx;
    double x1,x2;
       
    #define LS_COEFFICENTS
    #ifdef LS_COEFFICENTS
    std::ofstream file("generated_LS_coeff.h",std::ofstream::out);
    catima::Projectile p;
    std::vector<double> ve;
    std::stringstream sv1;
    std::stringstream sv2;
    std::stringstream sx1;
    std::stringstream sx2;
    
    file<<"//LS precalculated LS coefficients \n ";
    file<<"#include \"catima/storage.h\"";
    file<<"\n\n";
    file<<"#define LS_NUM_ENERGY_POINTS "<<num_LS_datapoints<<"\n";
    file<<"#define LS_MAX_Z "<<ATOMIC_WEIGHT_MAXZ<<"\n";
    file<<"namespace catima{\n\n";
    file<<"namespace ls_coefficients{\n\n";
    file<<"// relative difference of mass number for 2nd mass point ls_coefficients_ahi \n\n";
    file<<"constexpr double a_rel_increase=0.05;\n";
    file<<"constexpr double logEmin = "<<logEmin_forLS<<";\n";
    file<<"constexpr double logEmax = "<<logEmax_forLS<<";\n";
    file<<"// energy points array in MeV/u \n";
    file<<"constexpr EnergyTable<"<<num_LS_datapoints<<"> ls_energy_table("<<logEmin_forLS<<","<<logEmax_forLS<<");\n";
     
/*
    file<<"double ls_energy_points["<<num_LS_datapoints<<"]={";
    for(int i=0;i<num_LS_datapoints;i++){
        if(i!=0)file<<",";
        file<<ls_energy_table(i);
        
    }
    file<<"};\n";
  */  
    
    sv1<<"{";
    sv2<<"{";
    sx1<<"{";
    sx2<<"{";
    double beta2;
    for(int z=1;z<ATOMIC_WEIGHT_MAXZ-1;z++){
        p.Z = z;
        if(z!=1)sv1<<",";
        if(z!=1)sv2<<",";
        if(z!=1)sx1<<",";
        if(z!=1)sx2<<",";
        sv1<<"\n\t{";
        sv2<<"\n\t{";
        sx1<<"\n\t{";
        sx2<<"\n\t{";
        for(int i=0;i<num_LS_datapoints;i++){
            p.T = ls_energy_table(i);
            p.A = element_atomic_weight((int)p.Z);
            beta2 = catima::beta_from_T(p.T);
            beta2 = beta2*beta2;
            r = bethek_lindhard(p);
            x1 = bethek_lindhard_X(p);
            
            dx = p.A*0.05;
            p.A +=dx;
            r2 = bethek_lindhard(p);
            x2 = bethek_lindhard_X(p);
            if(i!=0)sv1<<",";
            sv1<<r;
            if(i!=0)sv2<<",";
            sv2<<r2;
            if(i!=0)sx1<<",";
            sx1<<x1;
            if(i!=0)sx2<<",";
            sx2<<x2;
        }
        sv1<<"}";
        sv2<<"}";
        sx1<<"}";
        sx2<<"}";
    }
    sv1<<"\n};";
    sv2<<"\n};";
    sx1<<"\n};";
    sx2<<"\n};";
    
    file<<"\n";
    file<<"//arrays dimensions are [z][energy], z=1 starts from index=0\n";
    file<<"\n//LS coefficient for A=atomic weight\n";
    file<<"double ls_coefficients_a[]["<<num_LS_datapoints<<"]=\n";
    file<<sv1.str();
    file<<"\n";
    
    file<<"\n//LS coefficient for A=atomic weight * 1.05\n";
    file<<"double ls_coefficients_ahi[]["<<num_LS_datapoints<<"]=\n";
    file<<sv2.str();
    file<<"\n";
    
    
    file<<"\n";
    file<<"//arrays dimensions are [z][energy], z=1 starts from index=0\n";
    file<<"\n//LS X coefficient (dE straggling) for A=atomic weight\n";
    file<<"double ls_X_coefficients_a[]["<<num_LS_datapoints<<"]=\n";
    file<<sx1.str();
    file<<"\n";
    
    file<<"\n//LS X coefficient for A=atomic weight * 1.05\n";
    file<<"double ls_X_coefficients_ahi[]["<<num_LS_datapoints<<"]=\n";
    file<<sx2.str();
    file<<"\n";
    
    file<<"\n}\n"; // end of ls_coefficient namespace
    file<<"\n}\n"; // end of catima namespace
    file.close();
    
    #endif
    
    return 0;
}
