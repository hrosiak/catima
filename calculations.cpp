#include <math.h>
#include <algorithm>
#include <cassert>
#include "catima/calculations.h"
#include "catima/constants.h"
#include "catima/data_ionisation_potential.h"
#include "catima/data_atima.h"
#include "catima/generated_LS_coeff.h"
#include "catima/nucdata.h"
#include "catima/storage.h"
#ifdef GLOBAL
extern "C"
{
    #include "global/globallib.h"
}
#endif


namespace catima{

double dedx_e(Projectile &p, const Target &t, const Config &c){
    double se = -1;
    if(p.T<=10){
        se = sezi_dedx_e(p,t);
    }
    else if(p.T>10 && p.T<30){
        double factor = 0.05 * ( p.T - 10.0 );
        se = (1-factor)*sezi_dedx_e(p,t) + factor*bethek_dedx_e(p,t,c);
    }
    else {
        se = bethek_dedx_e(p,t,c);
    }
    return se;
}


double dedx(Projectile &p, const Target &t, const Config &c){
    return dedx_e(p,t,c) + dedx_n(p,t);
}


double reduced_energy_loss_unit(const Projectile &p, const Target &t){
    double zpowers = pow(p.Z,0.23)+pow(t.Z,0.23);
    double asum = p.A + t.A;
    return 32.53*t.A*1000*p.T*p.A/(p.Z*t.Z*asum*zpowers); //projectile energy is converted from MeV/u to keV
}

double dedx_n(const Projectile &p, const Target &t){
    
    double zpowers = pow(p.Z,0.23)+pow(t.Z,0.23);
    double asum = p.A + t.A;
    double epsilon = 32.53*t.A*1000*p.T*p.A/(p.Z*t.Z*asum*zpowers); //projectile energy is converted from MeV/u to keV
    double sn=0;
    
    if(epsilon<=30){
        sn = log(1+(1.1383*epsilon))/ (2*(epsilon + 0.01321*pow(epsilon,0.21226) + 0.19593*pow(epsilon,0.5)));
    }
    else{
        sn = log(epsilon)/(2*epsilon);
    }
    sn = 100*8.4621*p.Z*t.Z*p.A*sn*Avogadro/(asum*zpowers*t.A);
    return sn;
}


double bethek_dedx_e(Projectile &p, const Target &t, const Config &c){
    assert(t.Z>0 && p.Z>0);
    assert(t.A>0 && p.A>0);
    if(p.T==0)return 0.0;
    double gamma=1.0 + p.T/atomic_mass_unit;
    double beta2=1.0-1.0/(gamma*gamma);
    assert(beta2>=0);
    double beta = sqrt(beta2);
    assert(beta>=0 && beta<1);
    //double zeta = 1.0-exp(-130.0*beta/pow(p.Z,2.0/3.0));
    //assert(zeta>=0);
    //double zp_eff = p.Z*zeta;
    double zp_eff = z_effective(p,t,c);
    
    assert(zp_eff>=0);
    double Ipot = ipot(t.Z);
    assert(Ipot>0);
    double f1 = dedx_constant*pow(zp_eff,2.0)*t.Z/(beta2*t.A);
    assert(f1>=0);
    double f2 = log(2.0*electron_mass*1000000*beta2/Ipot);
    double eta = beta*gamma;
    if(!(c.dedx&corrections::no_shell_correction) &&  eta>=0.13){ //shell corrections
        double cor = (+0.422377*pow(eta,-2)
                    +0.0304043*pow(eta,-4)
                    -0.00038106*pow(eta,-6))*1e-6*pow(Ipot,2) 
                  +(+3.858019*pow(eta,-2) 
                    -0.1667989*(pow(eta,-4))
                    +0.00157955*(pow(eta,-6)))*1.0e-9*pow(Ipot,3); 
        f2 = f2 -cor/t.Z;
    }
    f2+=2*log(gamma) -beta2;
    
    double barkas=1.0;
    if(!(c.dedx&corrections::no_barkas)){
        barkas = bethek_barkas(zp_eff,eta,t.Z);
        }
    
    double delta = bethek_density_effect(beta, t.Z);
    
    double LS = 0.0;
    if(!(c.dedx&corrections::no_lindhard)){
        //double LS = bethek_lindhard(p);
        LS = precalculated_lindhard(p);
        }
    
    double result  = (f2)*barkas + LS - delta/2.;
    result *=f1;
    
    return result;
}


double bethek_barkas(double zp_eff,double eta, double zt){
    double V2FVA[4]={0.33,0.30,0.26,0.23};
    double    VA[4]={1.,2.,3.,4.};
    double v1 = eta/(fine_structure*sqrt(zt));
    double v2fv;
    if(v1 >= 4){
        v2fv = 0.45 / sqrt(v1);
        }
    else if((v1 >= 1) && v1 < 4){//VALUES FROM THE JACKSON MC CARTHY FUNCTION //PHYS. REV. B 6 4131 P4136
        int i;
        for(i=1; i<4; i++){
	          if( VA[i] >= v1) break;
          }
	     v2fv = V2FVA[i-1]+(v1-VA[i-1])*(V2FVA[i]-V2FVA[i-1])/(VA[i]-VA[i-1]);
    }
    else{
        v2fv=0;
    }
    return 1.0+2.0 * zp_eff * v2fv /(v1*v1*sqrt(zt));
}

double bethek_density_effect(double beta, int zt){
    double gamma = 1/sqrt(1-(beta*beta));
    double x = log(beta * gamma) / 2.302585;
    int i = zt-1;
    double del = 0;
    
    if (x < density_effect::x0[i] ){
        if(density_effect::del_0[i] > 0.)del = density_effect::del_0[i] * pow(10.0,(2.*(x-density_effect::x0[i])));
        }
    else {
        del = 4.6052 * x - density_effect::c[i];
        if ( density_effect::x0[i]<= x &&  x <= density_effect::x1[i] ) del +=  density_effect::a[i] * pow((density_effect::x1[i] - x),density_effect::m[i]);
        }
    return del;
}


double bethek_lindhard(const Projectile &p){
    const double compton=3.05573356675e-3; // 1.18 fm / Compton wavelength
    double rho = exp(log(p.A)/3.0)*compton;
    double gamma=1.0 + p.T/atomic_mass_unit;
    double beta2=1.0-1.0/(gamma*gamma);
    double beta = sqrt(beta2);
    double eta = p.Z*fine_structure/beta;
    double beta_gamma_R = beta*gamma*rho;
    double sum = 0;
    int n=1;
    
    if(gamma < 10.0/rho){
        double dk[3];
        double dmk = 0;
        double dkm1 = 0;
        while(n<100){
            double k0 = n;
            int max = (n==1)?3:2;
            for(int i=0;i<max;i++){
                double k;
                if(i==0)k=k0;
                if(i==1)k=-k0 - 1.0;
                if(i==2)k=-k0;
                double l = (k>0)?k:-k-1.0;
                double signk = (k>0)?1:((k<0)?-1:0);
                double sk = sqrt(k*k-fine_structure*fine_structure*p.Z*p.Z);
                std::complex<double> cexir_n (k,-eta/gamma);
                std::complex<double> cexir_den (sk,-eta);
                std::complex<double> cexir = std::sqrt(cexir_n/cexir_den);
                std::complex<double> csketa (sk + 1.0, eta);
                std::complex<double> cpiske(0.0,(M_PI*(l-sk)/2.0) - lngamma(csketa).imag());
                std::complex<double> cedr = cexir*std::exp(cpiske);
                double H=0;
                
                //std::complex<double> ceds(0.0,0.0);
                // finite struct part
                std::complex<double> cmsketa (-sk + 1.0, eta);
                std::complex<double> cexis_den (-sk,-eta);
                std::complex<double> cexis = std::sqrt(cexir_n/cexis_den);
                std::complex<double> cpimske(0.0,(M_PI*(l+sk)/2.0) - lngamma(cmsketa).imag());
                std::complex<double> ceds = cexis*std::exp(cpimske);
                std::complex<double> cmbeta_gamma_R(0,-beta_gamma_R);
                std::complex<double> c2beta_gamma_R(0,2.0*beta_gamma_R);
                std::complex<double> c2sk_1 (2.0*sk+1,0);
                std::complex<double> cm2sk_1 (-2.0*sk+1,0);
                std::complex<double> clambda_r = cexir*std::exp(cmbeta_gamma_R)*hyperg(csketa,c2sk_1,c2beta_gamma_R);
                std::complex<double> clambda_s = cexis*std::exp(cmbeta_gamma_R)*hyperg(cmsketa,cm2sk_1,c2beta_gamma_R);
                std::complex<double> cGrGs = lngamma(cm2sk_1);
                double GrGs = clambda_r.imag()/clambda_s.imag();
                GrGs *= exp( lngamma(csketa).real() 
                             - lngamma(cmsketa).real() 
                             - lngamma(c2sk_1).real()
                             + cGrGs.real() 
                             + 2.0*sk*log(2.0*beta_gamma_R));
                if(cos(cGrGs.imag()) < 1.0)GrGs*=-1;
                if(fabs(GrGs)>1.0e-9){
                    double FrGr = sqrt((gamma-1)/(gamma+1)) * clambda_r.real()/clambda_r.imag();
                    double FsGs = sqrt((gamma-1)/(gamma+1)) * clambda_s.real()/clambda_s.imag();
                    double gz = -1.0*signk*(rho*gamma + 1.5*p.Z*fine_structure);
                    double z1 = -1.0*signk*p.Z;
                    double b0 = 1.0;
                    double a0 = (1.0 + 2.0*fabs(k))*b0/(rho-gz);
                    double a1 = 0.5*(gz+rho)*b0;
                    double an = a1;
                    double anm1 = a0;
                    double bnm1 = b0;
                    double asum = a0;
                    double bsum = b0;
                    double nn = 1.0;
                    while(fabs(anm1/asum)>1e-6 && fabs(anm1/asum)>1e-6){
                        double bn = ((rho-gz)*an + fine_structure*z1*anm1/2.0)/(2.0*nn+2.0*fabs(k)+1.0);
                        double anp1 = ((gz+rho)*bn - fine_structure*z1*bnm1/2.0)/(2.0*nn + 2.0);
                        asum += an;
                        bsum += bn;
                        nn += 1.0;
                        anm1 = an;
                        an = anp1;
                        bnm1 = bn;
                    }
                    double figi= (k>0) ? asum/bsum : bsum/asum;
                    H = (FrGr - figi)/(figi-FsGs)* GrGs;
                    
                }
                else 
                    H = 0;
                
                
                dk[i] = std::arg(cedr + H*ceds);
            }
            if(n>1)dk[2] = dmk;
        
            double sdm2 = sin(dk[2]-dk[1]);
            double term1 = k0*(k0+1.0)*sdm2*sdm2/(eta*eta*(2.0*k0 + 1.0));
            if(n>1){
                double sd2 = sin(dk[0]-dkm1);
                term1 += k0*(k0-1.0)*sd2*sd2/(eta*eta*(2.0*k0 - 1.0));
            }
            double sdd = sin(dk[0]-dk[2]);
            double term2 = k0*sdd*sdd/(eta*eta*(4.0*k0*k0 - 1.0));
            double term3 = term1 - 1.0/k0;
            sum += term2 + term3;
            n += 1;
            dmk = dk[1];
            dkm1 = dk[0];
            
        }// end of while n<100
        
    }
    else{ // ultrarelativistic limit
        sum = -log(beta_gamma_R) - 0.2;
    }   
    return sum + (0.5*beta2);
}

double bethek_lindhard_X(const Projectile &p){
    const double compton=3.05573356675e-3; // 1.18 fm / Compton wavelength
    double rho = exp(log(p.A)/3.0)*compton;
    double gamma=1.0 + p.T/atomic_mass_unit;
    double beta2=1.0-1.0/(gamma*gamma);
    double beta = sqrt(beta2);
    double eta = p.Z*fine_structure/beta;
    double beta_gamma_R = beta*gamma*rho;
    double sum = 0;
    int n=1;
    
    if(1){
        double dk[4];
        double dmk = 0;
        double dmkp1 = 0;
        double dkm1 = 0;
        double dkm2 = 0;
        
        while(n<200){
            double k0 = n;
            //int max = (n==1)?4:2;
            int max = 4;
            for(int i=0;i<max;i++){
                double k;
                if(i==0)k=k0;
                if(i==1)k=-k0 - 1.0;
                if(i==2 && n==1)k=-k0;
                if(i==3)k=-k0 - 2.0;
                double l = (k>0)?k:-k-1.0;
                double signk = (k>0)?1:((k<0)?-1:0);
                double sk = sqrt(k*k-fine_structure*fine_structure*p.Z*p.Z);
                std::complex<double> cexir_n (k,-eta/gamma);
                std::complex<double> cexir_den (sk,-eta);
                std::complex<double> cexir = std::sqrt(cexir_n/cexir_den);
                std::complex<double> csketa (sk + 1.0, eta);
                std::complex<double> cpiske(0.0,(M_PI*(l-sk)/2.0) - lngamma(csketa).imag());
                std::complex<double> cedr = cexir*std::exp(cpiske);
                double H=0;
                
                //std::complex<double> ceds(0.0,0.0);
                // finite struct part
                std::complex<double> cmsketa (-sk + 1.0, eta);
                std::complex<double> cexis_den (-sk,-eta);
                std::complex<double> cexis = std::sqrt(cexir_n/cexis_den);
                std::complex<double> cpimske(0.0,(M_PI*(l+sk)/2.0) - lngamma(cmsketa).imag());
                std::complex<double> ceds = cexis*std::exp(cpimske);
                std::complex<double> cmbeta_gamma_R(0,-beta_gamma_R);
                std::complex<double> c2beta_gamma_R(0,2.0*beta_gamma_R);
                std::complex<double> c2sk_1 (2.0*sk+1,0);
                std::complex<double> cm2sk_1 (-2.0*sk+1,0);
                std::complex<double> clambda_r = cexir*std::exp(cmbeta_gamma_R)*hyperg(csketa,c2sk_1,c2beta_gamma_R);
                std::complex<double> clambda_s = cexis*std::exp(cmbeta_gamma_R)*hyperg(cmsketa,cm2sk_1,c2beta_gamma_R);
                std::complex<double> cGrGs = lngamma(cm2sk_1);
                double GrGs = clambda_r.imag()/clambda_s.imag();
                GrGs *= exp( lngamma(csketa).real() 
                             - lngamma(cmsketa).real() 
                             - lngamma(c2sk_1).real()
                             + cGrGs.real() 
                             + 2.0*sk*log(2.0*beta_gamma_R));
                if(cos(cGrGs.imag()) < 1.0)GrGs*=-1;
                if(fabs(GrGs)>1.0e-9){
                    double FrGr = sqrt((gamma-1)/(gamma+1)) * clambda_r.real()/clambda_r.imag();
                    double FsGs = sqrt((gamma-1)/(gamma+1)) * clambda_s.real()/clambda_s.imag();
                    double gz = -1.0*signk*(rho*gamma + 1.5*p.Z*fine_structure);
                    double z1 = -1.0*signk*p.Z;
                    double b0 = 1.0;
                    double a0 = (1.0 + 2.0*fabs(k))*b0/(rho-gz);
                    double a1 = 0.5*(gz+rho)*b0;
                    double an = a1;
                    double anm1 = a0;
                    double bnm1 = b0;
                    double asum = a0;
                    double bsum = b0;
                    double nn = 1.0;
                    while(fabs(anm1/asum)>1e-6 && fabs(anm1/asum)>1e-6){
                        double bn = ((rho-gz)*an + fine_structure*z1*anm1/2.0)/(2.0*nn+2.0*fabs(k)+1.0);
                        double anp1 = ((gz+rho)*bn - fine_structure*z1*bnm1/2.0)/(2.0*nn + 2.0);
                        asum += an;
                        bsum += bn;
                        nn += 1.0;
                        anm1 = an;
                        an = anp1;
                        bnm1 = bn;
                    }
                    double figi= (k>0) ? asum/bsum : bsum/asum;
                    H = (FrGr - figi)/(figi-FsGs)* GrGs;
                    
                }
                else 
                    H = 0;
                
                
                dk[i] = std::arg(cedr + H*ceds);
            }
            if(n>1)dk[2] = dmk;
        
            double strterm1p = 0;
            double strterm1n = 0;
            double strterm2 = 0;
            double strterm3 = 0;
            double eta2 = eta*eta;
            
            double sdm2 = sin(dk[0]-dkm2);
            if(n>2){
                strterm1p = sdm2*sdm2*(k0-1)*(k0-2)/((2.0*k0 - 1.0)*(2.0*k0-3.0));
            }
            
            sdm2 = sin(dk[2]-dk[3]);
            strterm1n = sdm2*sdm2*(-k0-1)*(-k0-2)/((-2.0*k0 - 1.0)*(-2.0*k0-3.0));
            
            if(n>1){
                double sd2 = sin(dk[0]-dmkp1);
                strterm2 += (k0-1.0)*sd2*sd2/((2.0*k0 - 3.0)*(4.0*k0*k0 - 1.0));
            }
            
            double sdd = sin(dk[0]-dk[1]);
            strterm3 = sdd*sdd*(k0+1.0)*((1/(4.0*k0*k0 -1.0))+(1/(4*(k0+1.0)*(k0+1.0) - 1.0)))/(2.0*k0 + 1.0);
            
            //sum += k0*(strterm1p + strterm1n + (strterm2*2) + strterm3)/eta2;
            sum += k0*(strterm1p + strterm1n + (strterm2*2) + strterm3)/eta2;
            sum += - (2.0/k0);
            //std::cout<<n<<" "<<strterm1p<<" "<<strterm1n<<" "<<strterm2<<" "<<strterm3<<" "<<sum<<std::endl;
            n += 1;
            dmk = dk[1];
            dkm2 = dkm1;
            dkm1 = dk[0];
            dmkp1 = dk[2];
            
        }// end of while n<100
        
    }
    else{ // ultrarelativistic limit
        
    }
    return 2*bethek_lindhard(p) - sum - beta2;
    //return sum;
}



double sezi_p_se(double energy,const Target &t){
    double sp = -1;
    double e = 1000*energy; //e in keV/u
    int i = t.Z - 1;
    
    if(e<=25)e=25;
    //double sl = (proton_stopping_coef[i][0]*pow(e,proton_stopping_coef[i][1])) + (proton_stopping_coef[i][2]*pow(e,proton_stopping_coef[i][3]));
    //double sh = proton_stopping_coef[i][4]/pow(e,proton_stopping_coef[i][5]) * log( (proton_stopping_coef[i][6]/e) + (proton_stopping_coef[i][7]*e));
    double sl = (proton_stopping_coef[i][0]*catima::power(e,proton_stopping_coef[i][1])) + (proton_stopping_coef[i][2]*catima::power(e,proton_stopping_coef[i][3]));
    double sh = proton_stopping_coef[i][4]/catima::power(e,proton_stopping_coef[i][5]) * log( (proton_stopping_coef[i][6]/e) + (proton_stopping_coef[i][7]*e));
    sp = sl*sh/(sl+sh);
    e=1000*energy;
    if(e<=25){
        //sp *=(t.Z>6)?pow(e/25,0.45):pow(e/25,0.25);
        sp *=(t.Z>6)?catima::power(e/25,0.45):catima::power(e/25,0.25);
    }
    
    return 100*sp*Avogadro/t.A;
}

double sezi_dedx_e(const Projectile &p, const Target &t){
    double e=p.T*1000; // e in keV/u
    double se = 0;
    
    if(p.Z==1){
        return sezi_p_se(p.T,t);
    }
    else if(p.Z == 2){
        double a=0;
        double b=0;
        //double zeta = 0;
        
        if(e<=1)e=1;
        // He Zeff
        b = log(e);
        a = 0.2865 + b*(0.1266+ b*(-0.001429+ b*(0.02402 + b*(-0.01135 + b*0.001475))));
        double heh = 1.0 - exp(-std::min(30.,a));
        b = 7.6 - std::max(0., b);
        a = (1.0 + (0.007 + 0.00005*t.Z)*exp(- b*b ));
        heh *= a*a;
        //zeta = sqrt(heh);
        se = sezi_p_se(p.T,t)*heh*4.0; //scale proton stopping
        if(e==1)se*= sqrt(p.T*1000.0);  //vel proportional
        return se;
    }
    else{ // heavy ion
        double h1,h2,h3,h4;
        double a,q,b;
        double l1,l0,l;
        double YRmin = 0.130; //  YRmin = VR / ZP**0.67 <= 0.13 OR VR <= 1.0
        double VRmin = 1.0;
        double v=0;
        double vfermi = atima_vfermi[(int)t.Z-1];
        double yr=0;
        double zeta = 0;
        double se;
        
        v = sqrt(e/25.0)/vfermi;
        double v2=v*v;
        
        double vr = (v >= 1)? v*vfermi*(1.+ 1./(5.*v2)) : 3.0*vfermi/4.0*(1.0+v2*(2.0/3.0-v2/15.0));
        
        h1= 1./catima::power(p.Z,0.6667);
        yr = std::max(YRmin,vr*h1);
        yr = std::max(yr, VRmin*h1);
        
        //--  CALCULATE ZEFF
        a = -0.803*catima::power(yr,0.3) + 1.3167*catima::power(yr,0.6) + 0.38157*yr + 0.008983*yr*yr;
        q = std::min(1.0, std::max(0.0 , (1.0 - exp(-std::min(a, 50.0))))); //-- Q = IONIZATION LEVEL OF THE ION AT RELATIVE VELOCITY YR
        
        //-- IONIZATION LEVEL TO EFFECTIVE CHARGE
        h1 = 1./ catima::power(p.Z,0.3333);
        
        b = (std::min(0.43, std::max(0.32,0.12 + 0.025*p.Z)))*h1;
        l0 = (.8 - q * std::min(1.2,0.6 +p.Z/30.0))*h1;
        if(q < 0.2){
            l1 = 0;
            }
        else{
            if (q < std::max(0.0,0.9-0.025*p.Z)){
                l1 = b*(q-0.2)/fabs(std::max(0.0,0.9-0.025*p.Z)-0.2000001);
                }
	     	 else{
             	if(q < std::max(0.0,1.0 - 0.025*std::min(16.,p.Z))) l1 = b;
                else l1 = b*(1.0 - q)/(0.025*std::min(16.,p.Z));
                }
        }
        // calculate screening
        l = std::max(l1,l0*atima_lambda_screening[(int)p.Z-1]);
        h1 =4.0*l*vfermi/1.919;
        zeta = q + (1./(2.*(vfermi*vfermi)))*(1. - q)* log(1. + h1*h1);
         
         // ZP**3 EFFECT AS IN REF. 779?
        a = 7.6 - std::max(0.0, log(e));
        zeta = zeta*(1. + (1./(p.Z*p.Z))*(0.18 + .0015*t.Z)*exp(-a*a));
        
        h1= 1./catima::power(p.Z,0.6667);
        if (yr <= ( std::max(YRmin, VRmin*h1))){
            VRmin=std::max(VRmin, YRmin/h1);
            //--C        ..CALCULATE VELOCITY STOPPING FOR YR < YRmin
            double vmin =.5*(VRmin + sqrt(std::max(0.0,VRmin*VRmin - .8*vfermi*vfermi)));
            double eee = 25.0*vmin*vmin;
            double eval = 1;
            if((t.Z == 6) || (((t.Z == 14) || (t.Z == 32)) &&  (p.Z <= 19))) eval = 0.35;
            else eval = 0.5;
            
            h1 = zeta *p.Z;
            h4 = catima::power(e / eee,eval);
            se = sezi_p_se(eee*0.001,t) * h1*h1*h4;
            return se;
        }
        else {
            se = sezi_p_se(p.T,t)*catima::power(zeta*p.Z,2.0);
            return se;
        }
        return 0;
    }
};



double gamma_from_T(double T){
    return 1.0 + T/atomic_mass_unit;
};

double beta_from_T(double T){
    double gamma = gamma_from_T(T);
    return sqrt(1.0-1.0/(gamma*gamma));
}

double energy_straggling_firsov(double z1,double energy, double z2, double m2){
    double gamma = gamma_from_T(energy);
    double beta2=1.0-1.0/(gamma*gamma);
    double factor=4.8184E-3*pow(z1+z2,8.0/3.0)/m2;
    return factor*beta2/fine_structure/fine_structure;
    }

double angular_scattering_variance(Projectile &p, Target &t){
    if(p.T<=0)return 0.0;
    double e=p.T;
    double _p = sqrt(e*(e+2*atomic_mass_unit))*p.A;
    double beta = _p/((e+atomic_mass_unit)*p.A);
    double lr = radiation_length(t.Z,t.A);
    return 198.81 * pow(p.Z,2)/(lr*pow(_p*beta,2));
}

/// radioation lengths are taken frm Particle Data Group 2014
double radiation_length(int z, int m){
    double lr = 0;
    if(z==1){return 63.04;}
    if(z==2){return 94.32;}
    if(z==3){return 82.78;}
    if(z==4){return 65.19;}
    if(z==6){return 42.7;}
    if(z==7){return 37.99;}
    if(z==8){return 34.24;}
    if(z==9){return 32.93;}
    if(z==10){return 28.94;}
    if(z==13){return 24.01;}
    if(z==14){return 21.82;}
    if(z==17){return 19.28;}
    if(z==18){return 19.55;}
    if(z==22){return 16.16;}
    if(z==26){return 13.84;}
    if(z==29){return 12.86;}
    if(z==32){return 12.25;}
    if(z==50){return 8.82;}
    if(z==54){return 8.48;} 
    if(z==74){return 6.76;}
    if(z==78){return 6.54;}
    if(z==79){return 6.46;}
    if(z==82){return 6.37;}
    if(z==92){return 6.00;}
    
    double z2 = z*z;
    double z_13 = 1.0/pow(z,1./3.);
    double z_23 = z_13*z_13;
    double a2 = fine_structure*fine_structure*z2;
    double a4 = a2*a2;
    double a6 = a4*a2;
    lr= 716.405*m/(z2* (log(184.15*z_13) + log(1194.0*z_23)/z - -1.202*a2 + 1.0369*a4 - 1.008*a6/(1+a2) ) );
    return lr;
}


double precalculated_lindhard(const Projectile &p){
    double T = p.T;
    int z = (int)p.Z ;
    if(z>LS_MAX_Z)z=LS_MAX_Z;
    if(p.T<ls_coefficients::ls_energy_table(0))T=ls_coefficients::ls_energy_table(0);
    
    double da = (p.A - element_atomic_weight(z))/element_atomic_weight(z);
    z = z-1;
    //catima::Interpolator ls_a(ls_coefficients::ls_energy_points,ls_coefficients::ls_coefficients_a[z],LS_NUM_ENERGY_POINTS,interpolation_t::linear);
    //catima::Interpolator ls_ahi(ls_coefficients::ls_energy_points,ls_coefficients::ls_coefficients_ahi[z],LS_NUM_ENERGY_POINTS,interpolation_t::linear);
    //catima::Interpolator ls_a(ls_coefficients::ls_energy_table.values,ls_coefficients::ls_coefficients_a[z],LS_NUM_ENERGY_POINTS,interpolation_t::cspline);
    //catima::Interpolator ls_ahi(ls_coefficients::ls_energy_table.values,ls_coefficients::ls_coefficients_ahi[z],LS_NUM_ENERGY_POINTS,interpolation_t::cspline);
    double v1 = EnergyTable_interpolate(ls_coefficients::ls_energy_table,T,ls_coefficients::ls_coefficients_a[z]);
    double v2 = EnergyTable_interpolate(ls_coefficients::ls_energy_table,T,ls_coefficients::ls_coefficients_ahi[z]);

    //double dif = ls_ahi(T) - ls_a(T);
    //return ls_a(T)+(dif*da/ls_coefficients::a_rel_increase);
    double dif = v2 - v1;
    return v1+(dif*da/ls_coefficients::a_rel_increase);
    
}

double precalculated_lindhard_X(const Projectile &p){
    double T = p.T;
    int z = (int)p.Z ;
    if(z>LS_MAX_Z)z=LS_MAX_Z;
    //if(p.T<ls_coefficients::ls_energy_table(0))T=ls_coefficients::ls_energy_table(0);
    if(p.T<ls_coefficients::ls_energy_table(0))
		return 1.0;
    double da = (p.A - element_atomic_weight(z))/element_atomic_weight(z);
    z = z-1;
    
    //catima::Interpolator ls_X_a(ls_coefficients::ls_energy_table.values,ls_coefficients::ls_X_coefficients_a[z],LS_NUM_ENERGY_POINTS,interpolation_t::linear);
    //catima::Interpolator ls_X_ahi(ls_coefficients::ls_energy_table.values,ls_coefficients::ls_X_coefficients_ahi[z],LS_NUM_ENERGY_POINTS,interpolation_t::linear);
    double v1 = EnergyTable_interpolate(ls_coefficients::ls_energy_table,T,ls_coefficients::ls_X_coefficients_a[z]);
    double v2 = EnergyTable_interpolate(ls_coefficients::ls_energy_table,T,ls_coefficients::ls_X_coefficients_ahi[z]);
    
    //double dif = ls_X_ahi(T) - ls_X_a(T);
    //return ls_X_a(T)+(dif*da/ls_coefficients::a_rel_increase);
    double dif = v2 - v1;
    return v1+(dif*da/ls_coefficients::a_rel_increase);
}

double dedx_variance(Projectile &p, Target &t, const Config &c){
    double gamma = gamma_from_T(p.T);
    double cor=0;
    double beta = beta_from_T(p.T);
    double beta2 = beta*beta;
    //double zp_eff = z_effective(p,t,c);
    double zp_eff = z_eff_Pierce_Blann(p.Z, beta);
    double f = domega2dx_constant*pow(zp_eff,2)*t.Z/t.A;

    if(c.dedx_straggling == omega::atima){
        cor = 24.89 * pow(t.Z,1.2324)/(electron_mass*1e6 * beta2)*
			log( 2.0*electron_mass*1e6*beta2/(33.05*pow(t.Z,1.6364)));
	    cor = std::max(cor, 0.0 );
    }
	//double X = bethek_lindhard_X(p);
    double X = precalculated_lindhard_X(p);
    X *= gamma*gamma;
    if(p.T<30.0)
		return std::min(f*(X+cor), energy_straggling_firsov(p.Z, p.T, t.Z,t.A));
	else
		return f*(X+cor);
}

double z_effective(const Projectile &p,const Target &t, const Config &c){
    if(c.z_effective == z_eff_type::none){
        return p.Q;
    }

    double gamma=1.0 + p.T/atomic_mass_unit;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    if(c.z_effective == z_eff_type::pierce_blann){
        return z_eff_Pierce_Blann(p.Z, beta);
    }
    else if(c.z_effective == z_eff_type::anthony_landorf){
        return z_eff_Anthony_Landford(p.Z, beta, t.Z);
    }
    
    else if(c.z_effective == z_eff_type::hubert){
        return z_eff_Hubert(p.Z, p.T, t.Z);
    }
    else if(c.z_effective == z_eff_type::winger){
        return z_eff_Winger(p.Z, beta, t.Z);
    }
    else if(c.z_effective == z_eff_type::global){
        return z_eff_global(p.Z, p.T, t.Z);
    }
    else if(c.z_effective == z_eff_type::atima14){
        return z_eff_atima14(p.Z, p.T, t.Z);
    }
    else if(c.z_effective == z_eff_type::schiwietz){
        return z_eff_Schiwietz(p.Z, beta, t.Z);
    }
    else{
        return 0.0;
    }
}

double z_eff_Pierce_Blann(double z, double beta){
    return z*(1.0-exp(-0.95*fine_structure_inverted*beta/pow(z,2.0/3.0))); 
}

double z_eff_Anthony_Landford(double pz, double beta, double tz){
    double B = 1.18-tz*(7.5e-03 - 4.53e-05*tz);
    double A = 1.16-tz*(1.91e-03 - 1.26e-05*tz);
    return pz*(1.0-(A*exp(-fine_structure_inverted*B*beta/pow(pz,2.0/3.0)))); 
}

double z_eff_Hubert(double pz, double E, double tz){
    double lntz = log(tz);
    double x1,x2,x3,x4;

    if(E<2.5)
        return 0.0;
    if(tz == 4){
        x1 = 2.045 + 2.0*exp(-0.04369*pz);
        x2 = 7.0;
        x3 = 0.2643;
        x4 = 0.4171;
    }
    else if(tz==6){
        x1 = 2.584 + 1.91*exp(-0.03958*pz);
        x2 = 6.933;
        x3 = 0.2433;
        x4 = 0.3969;
    }
    else{
        x1 = (1.164 + 0.2319*exp(-0.004302*tz)) + 1.658*exp(-0.0517*pz);
        x2 = 8.144 + 0.09876*lntz;
        x3 = 0.314 + 0.01072*lntz;
        x4 = 0.5218 + 0.02521*lntz;
    }
     
    return pz*(1-x1*exp(-x2*catima::power(E,x3)*catima::power(pz,-x4)));
}

double z_eff_Winger(double pz, double beta, double tz){
	double c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13;
	double x, lnz, lnzt, a0,a1,a2,a3,a4; 

	c0 = 0.4662;
	c1 = 0.5491;
	c2 = 0.7028;
	c3 = 0.1089;
	c4 = 0.001644;
	c5 = 0.5155;
	c6 = 0.05633;
	c7 = 0.005447;
	c8 = 0.8795;
	c9 = 1.091;
	c10= 0.0008261;
	c11= 2.848;
	c12= 0.2442;
	c13= 0.00009293;
	// the following is from Helmut to correct winger formula for H and He target
    if(tz==1){
        tz = 2.6;
    }
    if(tz==2){
        tz = 2.8;
    }
	x = beta /0.012 /pow(pz,0.45);
	lnz =log(pz);
	lnzt=log(tz);

	a0 = -c0;
	a1 = -c1  * exp( c2 *lnz - c3 *lnz*lnz +c4*lnz*lnz*lnz -c5*lnzt + c6 *lnzt*lnzt);
	a2 =  c7  * exp( c8 *lnz - c9 *lnzt);
	a3 = -c10 * exp( c11*lnz - c12*lnz*lnz*lnz);
	a4 = -c13;

    return pz * (1. - exp(a0 +a1*x +a2*x*x +a3*x*x*x +a4*x*x*x*x));	
	}

double z_eff_global(double pz, double E, double tz){
    if(E>2000)
        return pz;
    else
        #ifdef GLOBAL
        return global_qmean(pz, tz, E);
        #else
        return -1;
        #endif
}

double z_eff_Schiwietz(double pz, double beta, double tz){
    double scaled_velocity;
    double c1, c2;
    double x;

    scaled_velocity = catima::power(pz,-0.543)*beta/bohr_velocity;
    c1 = 1-0.26*exp(-tz/11.0)*exp(-(tz-pz)*(tz-pz)/9);
    c2 = 1+0.030*scaled_velocity*log(tz);
    x = c1*catima::power(scaled_velocity/c2/1.54,1+(1.83/pz));
    return pz*((8.29*x) + (x*x*x*x))/((0.06/x) + 4 + (7.4*x) + (x*x*x*x) );

}

double z_eff_atima14(double pz, double T, double tz){
    double qpb;
    double qhigh,qwinger,qglobal,qmean=0;
    double c1 = 1.4;
    double c2 = 0.28;
    double c3 = 0.04;
    double beta = beta_from_T(T);
    double gamma = gamma_from_T(T);
    double emax, emin;
    qpb = z_eff_Pierce_Blann(pz,beta);
    #ifdef GLOBAL

    if(T>30.0 && T<1500.0 && pz>28){
        qglobal = z_eff_global(pz,T,tz);
        qglobal = (qglobal - qpb)*c1/catima::power(tz+1,c2)*(1-exp(-c3*T)) + qpb;
    }
    
    emax = 1500.;
    emin = 1000.;
    if(T>emax){
        qhigh = pz * (1.0-exp(-180.0*beta*catima::power(gamma,0.18)*catima::power(pz,-0.82)*catima::power(tz,0.1)));
        qmean = qhigh;
    }
    else if(T>=emin && T<=emax){
        qhigh = pz * (1.0-exp(-180.0*beta*catima::power(gamma,0.18)*catima::power(pz,-0.82)*catima::power(tz,0.1)));
        if(pz<=28){
            qwinger = z_eff_Winger(pz,beta,tz);
            qmean = ((emax-T)*qwinger + (T-emin)*qhigh)/(emax-emin);
        }
        else{
            qmean = ((emax-T)*qglobal + (T-emin)*qhigh)/(emax-emin);
        }
    }
    else{
        emax = 70.0;
        emin = 30.0;
        if(pz<=28){
            qwinger = z_eff_Winger(pz,beta,tz);
            qmean = qwinger;
        }
        else{
            if(T>=emax){
                qmean = qglobal;
                }
            else if(T<emin){
                qwinger = z_eff_Winger(pz,beta,tz);
                qmean = qwinger;
                }
            else
                {
                qwinger = z_eff_Winger(pz,beta,tz);
                qmean = ((emax-T)*qwinger + (T-emin)*qglobal)/(emax-emin);
                }
        }
    }
    #endif
    return qmean;
    }

std::complex<double> hyperg(const std::complex<double> &a,
                                                const std::complex<double> &b,
                                                const std::complex<double> &z){
        double dm = 0.0;
        std::complex<double> term(1.0, 0.0);
        std::complex<double> sumterm(1.0, 0.0);
        std::complex<double> previousterm;
        do {
            previousterm = term;
            dm += 1.0;
            std::complex<double> Cm(dm-1.0, 0.0);
            term = previousterm * ((a + Cm)/(b + Cm)) * (z/dm);
            sumterm += term;
        } while( std::abs(term) > 1.0e-6 && std::abs(previousterm) > 1.0e-6 );
        return(sumterm);
}

std::complex<double> lngamma( const std::complex<double> &z )
    {
        const static double coeff[6] = {76.18009172947146,
            -86.50532032941677,
            24.01409824083091,
            -1.231739572450155,
            0.1208650973866179e-2,
            -0.5395239384953e-5};
        double x, y;
        if(z.real() > 0) {
            x=z.real()-1.0;
            y=z.imag();
        } else {
            x=-z.real();
            y=-z.imag();
        }
        double r = sqrt((x+5.5)*(x+5.5)+y*y);
        double aterm1=y*log(r);
        double aterm2=(x+0.5)*atan2(y,(x+5.5))-y;
        double lterm1=(x+0.5)*log(r);
        double lterm2=-y*atan2(y,(x+5.5)) - (x+5.5) + 0.5*log(2.0*M_PI);
        double num=0.0;
        double denom=1.000000000190015;
        for(int j=1;j<7;j++){
            double fj=(double)j;
            double cterm=coeff[j-1]/((x+fj)*(x+fj)+y*y);
            num+=cterm;
            denom+=(x+fj)*cterm;
        }
        num*=-y;
        double aterm3=atan2(num,denom);
        double lterm3 = 0.5*log(num*num + denom*denom);
        std::complex<double> result(lterm1+lterm2+lterm3,aterm1+aterm2+aterm3);
        if(z.real() < 0){
            std::complex<double> lpi(log(M_PI), 0.0);
            result = lpi - (result + std::log(std::sin(M_PI*z)));
        }
        return(result);
    }

}
