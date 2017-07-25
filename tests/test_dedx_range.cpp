#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <utility>
#include <vector>
#include "catima/catima.h"

using std::cout;
using std::endl;

double logEmax = catima::logEmax;
double logEmin = catima::logEmin;

template<int N>
    struct EnergyTable{
	constexpr EnergyTable():values(),num(N){
	    for(auto i=0;i<N;i++){
		values[i]=exp(M_LN10*(logEmin + ((double)i)*(logEmax-logEmin)/(N - 1.0)));
		}
	    }
	double operator()(int i)const{return values[i];}
	double values[N];
	std::size_t num;
    };
    

bool expect(bool r,std::string title=""){
    if(r){
        if(title.size()>0)
            cout << "\033[1;32mpassed\033[0m "<<title<<endl;
    }
    else
    {
        cout << "\033[1;31m failed\033[0m "<<title<<endl;
    }
    return r;
}

double rdif(double v1, double v2){
    return (v1-v2)/v2;
}

void comp_dedx(catima::Projectile p, catima::Material t, double epsilon = 0.001, bool fout=false, const char* fname=""){
    double dif;
    //struct atima_results results;
    int tz = t.get_element(0).first.Z;
    double ta = t.get_element(0).first.A;
    bool res;
    cout<<"-----------------------------"<<endl;
    cout<<"projectile: A = "<<p.A<<", Z = "<<p.Z<<endl;
    cout<<"target: A = "<<ta<<", Z = "<<tz<<endl;
    EnergyTable<100> etable;
    std::ofstream f;
    if(fout){
        f.open(fname,std::ofstream::out | std::ofstream::trunc);
    }
    double v1, v2;
    for(auto _e : etable.values){
        p.T = _e;
        v1 = catima::dedx(p,_e,t);
        auto ares = catima::calculate(p(_e),t);
        v2 = ares.dEdxi;
        dif = rdif(v1,v2);
        cout<<"T = "<<_e<<endl;
        cout<<v1 << " vs "<<v2<<" -> rel. dif = "<<dif<<endl;
        if(fout){
            f<<_e<<" "<<v1<<" "<<v2<<std::endl;
        }
        res = expect(fabs(dif)<epsilon,"");    
        if(!res)exit(0);
    }
    if(fout){
        f.close();
    }
}



int main(){
    catima::Projectile p{1,1};
    catima::Projectile c{12,6};
    catima::Projectile u{238,92};
    
    catima::Material c_target(12.011,6);
    catima::Material h_target(1.00794,1);
    catima::Material pb_target(207.2,82);

    catima::Material water(
                {{1,1,2},
                {16,8,1}
                });
    
    c_target.thickness(0.5);
    c_target.density(2.0);

    h_target.density(0.001);
    h_target.thickness(10.0);

    pb_target.density(11.0);
    pb_target.thickness(0.1);

    water.density(1.0);
    water.thickness(1.0);
    
    std::vector<std::pair<double,int>> projectiles = {{1,1},{4,2},{12,6},{238,92}};

    bool fwrite = false;
    double eps = 0.01;
    comp_dedx(p,c_target,eps,fwrite,"dedx_p_c.dat");
    comp_dedx(p,h_target,eps,fwrite,"dedx_p_h.dat");
    comp_dedx(p,water,eps,fwrite,"dedx_p_water.dat");
    comp_dedx(p,pb_target,eps,fwrite,"dedx_p_pb.dat");
    
    comp_dedx(c,c_target,eps,fwrite,"dedx_c_c.dat");
    comp_dedx(c,h_target,eps,fwrite,"dedx_c_h.dat");
    comp_dedx(c,water,eps,fwrite,"dedx_c_water.dat");
    comp_dedx(c,pb_target,eps,fwrite,"dedx_c_pb.dat");
    
    comp_dedx(u,c_target,eps,fwrite,"dedx_u_c.dat");
    comp_dedx(u,h_target,eps,fwrite,"dedx_u_h.dat");
    comp_dedx(u,water,eps,fwrite,"dedx_u_water.dat");
    comp_dedx(u,pb_target,eps,fwrite,"dedx_u_pb.dat");

      
    return 0;
}
