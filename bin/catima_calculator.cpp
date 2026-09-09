 #include <math.h>

#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include "catima/catima.h"
#include "catima/nucdata.h"
#include <nlohmann/json.hpp>

using namespace std;
using namespace catima;
using json = nlohmann::json;

void help(){
        std::cout<<"usage: catima_calculator config_file.json\n";
}

inline std::vector<double> linspace_vector(double a, double b, unsigned int num){
    std::vector<double> res;
    if(num>=2 && a<b){
    res.resize(num);
    double step = (b-a)/(num-1);
    for(unsigned int i=0;i<(num-1);i++){
        res[i]=a+(i*step);
        }
    res[num-1] = b;
    }
    return res;
    }

json load_json(const char *fname);
char* getCmdOption(char ** begin, char ** end, const std::string & option);
Material json_material(json &j);

int main( int argc, char * argv[] )
{
    Projectile projectile;
    Layers layers;
    std::vector<double> energies;
    Config conf;
    
    if(argc == 1 ){
        help();
        return 0;
    }
    try{
        auto j = load_json(argv[1]);
            
            // load projectile data
            if(j.count("projectile")){
                if(j["projectile"].is_array()){
                    projectile.A = j["projectile"].at(0).get<double>();
                    projectile.Z = j["projectile"].at(1).get<double>();
                }
            }
            else{
                throw std::invalid_argument("projectile field is missing");
            }   

            // load energy data
            if(j.count("energy")){
                auto e = j.at("energy");
                if(e.is_number()){
                    energies.push_back(j["energy"].get<double>());
                    }
                if(e.is_string()){
                    double _e = std::stod(j["energy"].get<std::string>());
                    energies.push_back(_e);
                    }
                if(e.is_array()){
                    for(auto &el:e){
                        if(el.is_number())
                        energies.push_back(el.get<double>());
                        }
                    }
                if(e.is_object()){
                    if(e.count("min")>0 && e.count("max")>0 && (e.count("num")>0 || e.count("step")>0)){
                        double emin = e["min"].get<double>();
                        double emax = e["max"].get<double>();
                        int num=0;
                        if(e.count("step")){
                            num = 1+(emax-emin)/e["step"].get<double>();
                        }
                        if(e.count("num")){
                            num = e["num"].get<int>();
                        }
                        if(num>2)
                            energies  =  linspace_vector(emin,emax,num);
                    }
                }
                }
                else{
                    throw std::invalid_argument("energy field is missing");
                    }
         
            if(j.count("material")){
                auto e = j.at("material");
                if(e.is_array()){
                    for(auto& entry : e){
                        if(!entry.is_object()){
                            throw std::invalid_argument("material error");
                            }
                        layers.add(json_material(entry));
                        }
                    }
                if(e.is_object()){
                    layers.add(json_material(e));
                    }
                }
            else{
                throw std::invalid_argument("material field is missing");
                }
            if(j.count("config")>0){
                auto e = j["config"];
                if(e.is_string()){
                    std::string cstr = e.get<std::string>();
                    if(cstr=="atimav1.3"){
                        conf.z_effective = z_eff_type::pierce_blann;
                        cout<<"using config: Atima v1.3\n";
                    }
                    if(cstr=="atimav1.4"){
                        conf.z_effective = z_eff_type::atima14;
                        cout<<"using config: Atima v1.4\n";
                    }
                }
            }
                
        } // end of try 
        catch(...){
            cout<<"Could not load the config file"<<"\n";
            return 0;
        }


    if(layers.num()==0){
        cout<<"no material specified\n";
        return 0;
        }

    if(energies.size()==0){
        cout<<"no energy specified\n";
        return 0;
    }
    
    cout<<"******** CAtima calculator ********\n";
    cout<<"Projectile: A = "<<projectile.A<<", Z = "<<projectile.Z<<"\n";
    cout<<"Materials:\n";
    for(unsigned int i=0;i<layers.num();i++){
        cout<<"#"<<i;
        cout<<": density = "<<layers[i].density()<<" g/cm3";
        cout<<", thickness = "<<layers[i].thickness()<<" g/cm2";
        cout<<"\n";
        }
    
    for(double e:energies){ 
        cout<<"-------- T = "<<e<<" MeV/u -------\n";
        projectile.T = e;
        auto res = calculate(projectile,layers);
        for(unsigned int i=0;i<res.results.size();i++){
            auto entry = res.results[i];
            cout<<"material #"<<i<<":\n";
            cout<<"\tEin = "<<entry.Ein<< " MeV/u\n";
            cout<<"\tEout = "<<entry.Eout<<" MeV/u\n";
            cout<<"\tsigma_E = "<<entry.sigma_E<<" MeV\n";
            cout<<"\tEloss = "<<entry.Eloss<<" MeV\n";
            cout<<"\trange = "<<entry.range<<" g/cm2\n";
            cout<<"\tsigma_r = "<<entry.sigma_r<<" g/cm2\n";
            cout<<"\tsigma_a = "<<entry.sigma_a<<" rad\n";
            cout<<"\tdEdx(Ein) = "<<entry.dEdxi<<" MeV/g/cm2\n";
            cout<<"\tTOF = "<<entry.tof<<" ns\n";
            }
        cout<<"total:\n";
        cout<<"\tEout = "<<res.total_result.Eout<<" MeV/u\n";
        cout<<"\tBeta = "<<beta_from_T(res.total_result.Eout)<<"\n";
        cout<<"\tGamma = "<<gamma_from_T(res.total_result.Eout)<<"\n";
        cout<<"\tP = "<<p_from_T(res.total_result.Eout, projectile.A)<<" MeV/c\n";
        cout<<"\tEloss = "<<res.total_result.Eloss<<" MeV\n";
        cout<<"\tsigma_E = "<<res.total_result.sigma_E<<" MeV\n";
        cout<<"\tsigma_a = "<<res.total_result.sigma_a<<" rad\n";
        cout<<"\tTOF = "<<res.total_result.tof<<" ns\n";
    }

    return 1;
}



json load_json(const char *fname){
    std::vector<std::string> res;
    std::ifstream jfile(fname,std::ifstream::in);

    if(!jfile){
        throw std::invalid_argument("Could not open config file");
    }
    
    std::string content;
    jfile.seekg(0, std::ios::end);
    content.resize(jfile.tellg());
    jfile.seekg(0, std::ios::beg);
    jfile.read(&content[0], content.size());
    jfile.close();
    
    try{
        auto j = json::parse(content);
        return j;
    }
    catch(...){
        cout<<"JSON parsing error\n";
        throw std::invalid_argument("Could not parse json file");
    }
};

Material json_material(json &j){
    if(!j.is_object()){
        throw std::invalid_argument("Wrong material definition");
        }
    try{
        double a=0;
        int z=0;
        double ipot=0.0;
        double density=0.0;
        double th=0.0;
        if(j.count("density")>0){
            density = j["density"].get<double>();
            }
        if(j.count("thickness")>0){
            th = j["thickness"].get<double>();
            }
        if(j.count("Ipot")>0){
            ipot = j["Ipot"].get<double>();
            }
        if(j.count("A")>0){
            a = j["A"].get<double>();
            }
        if(j.count("Z")>0){
            z = j["Z"].get<int>();
            }
        if(z<=0){
            cout<<"Z="<<z<<"\n";
            throw std::invalid_argument("Could not parse json file (material section)");
            }
        
        if(density<=0){
            density = element_density(z);
            if(density<=0)cout<<"Warning: material density = "<<density<<"\n";
            }
        if(th<=0){
            cout<<"Warning: material thickness = "<<th<<"\n";
            }
        
        if(z<200){
            Material m(a,z,density,th);
            if(ipot>0)m.I(ipot);
            return m;
            }
        else{
            Material m = get_material(z);
            m.thickness(th);
            return m;
            }
        
        }
    catch(...){
        cout<<"JSON parsing error: material definition\n";
        throw std::invalid_argument("Could not parse json file");
        }
    }

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return nullptr;
}

