
namespace catima{

/**
  * return srim stopping power 
  * @param pZ - projectile Z
  * @param pZ - material Z
  * @param energy - projectile energy in MeV/u unit
  * @param use_v95 - use srim proton coefficient from version 95, otherwise version 85 will be used
  */
double srim_dedx_e(int pZ, int tZ, double energy, bool use_v95=1);

/**
  * return SRIM proton stopping power
  * @param Z - proton number of material
  * @param energy - energy per nuclein in MeV/u
  */
double p_se(int Z, double energy);
double p_se95(int Z, double energy);
}
