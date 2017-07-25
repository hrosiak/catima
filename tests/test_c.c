#include <stdio.h>
#include <math.h>
#include "catima/cwrapper.h"

void expect(int r,const char *title){
    if(r){
        printf("\033[1;32mpassed\033[0m %s\n",title);
    }
    else
    {
        printf("\033[1;31mfailed\033[0m %s\n",title);
    }
};

int main()
{
    CatimaResult r;
    double dif;
    r = catima_calculate(12,6,500,12,6,1.0,2.0);

    dif = r.dEdxi - 87.78;
    expect(fabs(dif)<0.5,"dEdx");

    dif = r.Eloss - 88.05;
    expect(fabs(dif)<0.5,"Eloss");

    dif = r.sigma_a - 1.31e-3;
    expect(fabs(dif)<0.5,"sigma_a");

    dif = r.sigma_E - 0.18;
    expect(fabs(dif)<0.01,"sigma_e");

    dif = r.range - 43.7;
    expect(fabs(dif)<0.5,"range");
    
    r = catima_calculate(12,6,1000,0,206,1.0,1.0);
    dif = r.Eloss - 80.75;
    expect(fabs(dif)<1,"Eloss");

    return 1.0;
}