#ifndef PHYDRO_TEMPERATURE_DEPENDENCIES_H
#define PHYDRO_TEMPERATURE_DEPENDENCIES_H

#include <cmath>


namespace phydro{

//-----------------------------------------------------------------------
// Input:    - float, annual atm. CO2, ppm (co2)
//           - float, monthly atm. pressure, Pa (patm)
// Output:   - ca in units of Pa
// Features: Converts ca (ambient CO2) from ppm to Pa.
//-----------------------------------------------------------------------
inline float co2_to_ca(float co2, float patm ){
  float ca = ( 1.0e-6 ) * co2 * patm;         // Pa, atms. CO2
  return ca;
}



//-----------------------------------------------------------------------
// Output:   Factor fv to correct for instantaneous temperature response
//           of Vcmax for:
//
//               Vcmax(temp) = fv * Vcmax(25 deg C) 
//
//-----------------------------------------------------------------------
inline float calc_ftemp_arrhenius(float tk, float dha, float tkref = 298.15 ){

  //// Note that the following two forms are equivalent:
  //// ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
  //// ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )

  float kR   = 8.3145;     // Universal gas constant, J/mol/K

  float ftemp = exp( dha * (tk - tkref) / (tkref * kR * tk) );

  return ftemp;
}



//-----------------------------------------------------------------------
// Input:    - float, air temperature, deg C (temp)
//           - float, atmospheric pressure, Pa (patm)
// Output:   float, Pa (mmk)
// Features: Returns the temperature & pressure dependent Michaelis-Menten
//           coefficient, K (Pa).
// Ref:      Bernacchi et al. (2001), Improved temperature response 
//           functions for models of Rubisco-limited photosynthesis, 
//           Plant, Cell and Environment, 24, 253--259.
//-----------------------------------------------------------------------
inline float calc_kmm(float tc, float patm ) {
  float dhac   = 79430;      // (J/mol) Activation energy, Bernacchi et al. (2001)
  float dhao   = 36380;      // (J/mol) Activation energy, Bernacchi et al. (2001)
  float kco    = 2.09476e5;  // (ppm) O2 partial pressure, Standard Atmosphere

  //// k25 parameters are not dependent on atmospheric pressure
  float kc25 = 39.97;   // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  float ko25 = 27480;   // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  float tk = tc + 273.15;

  float kc = kc25 * calc_ftemp_arrhenius( tk, dhac );
  float ko = ko25 * calc_ftemp_arrhenius( tk, dhao );

  float po  = kco * (1e-6) * patm;         // O2 partial pressure
  float kmm = kc * (1.0 + po/ko);

  return kmm;
}




//-----------------------------------------------------------------------
// Input:    - elevation, m (elv)
// Output:   - float, atmospheric pressure at elevation 'elv', Pa (patm)
// Features: Returns the atmospheric pressure as a function of elevation
//           and standard atmosphere (1013.25 hPa)
// Depends:  - connect_sql
//           - flux_to_grid
//           - get_data_point
//           - get_msvidx
// Ref:      Allen et al. (1998)
//-----------------------------------------------------------------------
inline float calc_patm(float elv){

  // Define constants:
  float kPo = 101325;    // standard atmosphere, Pa (Allen, 1973)
  float kTo = 298.15;    // base temperature, K (Prentice, unpublished)
  float kL  = 0.0065;    // temperature lapse rate, K/m (Allen, 1973)
  float kG  = 9.80665;   // gravitational acceleration, m/s^2 (Allen, 1973)
  float kR  = 8.3145;    // universal gas constant, J/mol/K (Allen, 1973)
  float kMa = 0.028963;  // molecular weight of dry air, kg/mol (Tsilingiris, 2008)

  // Convert elevation to pressure, Pa:
  float patm = kPo* pow(1.0 - kL*elv/kTo, kG*kMa/(kR*kL));
  
  return patm;
}




//-----------------------------------------------------------------------
// Input:    float, air temperature, degrees C (tc)
// Output:   float, gamma-star, Pa (gammastar)
// Features: Returns the temperature-dependent photorespiratory 
//           compensation point, Gamma star (Pascals), based on constants 
//           derived from Bernacchi et al. (2001) study.
// Ref:      Bernacchi et al. (2001), Improved temperature response 
//           functions for models of Rubisco-limited photosynthesis, 
//           Plant, Cell and Environment, 24, 253--259.
//-----------------------------------------------------------------------
inline float calc_gammastar(float tc, float patm ) {
  float dha    = 37830;       // (J/mol) Activation energy, Bernacchi et al. (2001)
  float gs25_0 = 4.332;       // Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  float gammastar25 = gs25_0 * patm / calc_patm(0.0);

  float tk = tc + 273.15;
  float gammastar = gammastar25 * calc_ftemp_arrhenius( tk, dha );

  return gammastar;
}


//-----------------------------------------------------------------------
// arguments
// tcleaf: temperature (degrees C)
// tref: is 'to' in Nick's set it to 25 C (=298.15 K in other cals)
//
// function return variable
// fv: temperature response factor, relative to 25 deg C.
//
// Output:   Factor fv to correct for instantaneous temperature response
//           of Vcmax for:
//
//               Vcmax(temp) = fv * Vcmax(25 deg C) 
//
// Ref:      Wang Han et al. (in prep.)
//-----------------------------------------------------------------------
inline float calc_ftemp_inst_vcmax(float tcleaf, float tcgrowth, float tcref = 25.0 ){
  // loal parameters
  float Ha    = 71513;  // activation energy (J/mol)
  float Hd    = 200000; // deactivation energy (J/mol)
  float Rgas  = 8.3145; // universal gas constant (J/mol/K)
  float a_ent = 668.39; // offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
  float b_ent = 1.07;   // slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
  
  float tkref = tcref + 273.15;  // to Kelvin

  // conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C. 
  float tkleaf = tcleaf + 273.15;

  // calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
  float dent = a_ent - b_ent * tcgrowth;   // 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
  float fva = calc_ftemp_arrhenius( tkleaf, Ha, tkref );
  float fvb = (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) );
  float fv  = fva * fvb;

  return fv;
}

inline float calc_ftemp_vcmax_bernacchi(double tc){
  double dha = 65.33e3;
  double c = 26.35;
  double R = 8.31;
  return exp(c - dha/R/(tc+273.16));
}

//-----------------------------------------------------------------------
// arguments:
// tc                  // temperature (degrees C)
// function return variable:
// fr                  // temperature response factor, relative to 25 deg C.
// Output:   Factor fr to correct for instantaneous temperature response
//           of Rd (dark respiration) for:
//
//               Rd(temp) = fr * Rd(25 deg C) 
//
// Ref:      Heskel et al. (2016) used by Wang Han et al. (in prep.)
//-----------------------------------------------------------------------
inline float calc_ftemp_inst_rd(float tc){
  // loal parameters
  float apar = 0.1012;
  float bpar = 0.0005;
  float tk25  = 298.15; // 25 deg C in Kelvin

  // conversion of temperature to Kelvin
  float tk = tc + 273.15;

  float fr = exp( apar * (tc - 25.0) - bpar * (tc*tc - 25.0*25.0) );
 
  return fr; 
}


//-----------------------------------------------------------------------
// Input:    - float, air temperature (tc), degrees C
//           - float, atmospheric pressure (p), Pa
// Output:   float, density of water, kg/m^3
// Features: Calculates density of water at a given temperature and 
//           pressure using the Tumlirz Equation
// Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of 
//           pure water and sea water, Tech. Rept., Marine Physical 
//           Laboratory, San Diego, CA.
//-----------------------------------------------------------------------
inline float calc_density_h2o(float tc, float p){

  // Calculate lambda, (bar cm^3)/g:
  float my_lambda = 1788.316 + 
          21.55053*tc + 
        -0.4695911*tc*tc + 
     (3.096363e-3)*tc*tc*tc + 
    -(7.341182e-6)*tc*tc*tc*tc;

  // Calculate po, bar
  float po = 5918.499 + 
           58.05267*tc + 
         -1.1253317*tc*tc + 
     (6.6123869e-3)*tc*tc*tc + 
    -(1.4661625e-5)*tc*tc*tc*tc;

  // Calculate vinf, cm^3/g
  float vinf = 0.6980547 +
    -(7.435626e-4)*tc +
     (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
     (9.829576e-9)*tc*tc*tc*tc +
   -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
   -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
     (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
   -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc;

  // Convert pressure to bars (1 bar = 100000 Pa)
  float pbar = (1e-5)*p;
  
  // Calculate the specific volume (cm^3 g^-1):
  float v = vinf + my_lambda/(po + pbar);

  // Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  float rho = (1e3/v);

  return rho; 
}


//-----------------------------------------------------------------------
// Input:    - float, ambient temperature (tc), degrees C
//           - float, ambient pressure (p), Pa
// Return:   float, viscosity of water (mu), Pa s
// Features: Calculates viscosity of water at a given temperature and 
//           pressure.
// Depends:  density_h2o
// Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V. 
//           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New 
//           international formulation for the viscosity of H2O, J. Phys. 
//           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
//-----------------------------------------------------------------------
// [[Rcpp::export(name="cpp_calc_viscosity_h2o")]]
inline float calc_viscosity_h2o(float tc, float p) {

  // Define reference temperature, density, and pressure values:
  float tk_ast  = 647.096;    // Kelvin
  float rho_ast = 322.0;      // kg/m^3
  float mu_ast  = 1e-6;       // Pa s

  // Get the density of water, kg/m^3
  float rho = calc_density_h2o(tc, p);

  // Calculate dimensionless parameters:
  float tbar  = (tc + 273.15)/tk_ast;
  float tbarx = pow(tbar, 0.5);
  float tbar2 = tbar*tbar;
  float tbar3 = tbar*tbar*tbar;
  float rbar  = rho/rho_ast;

  // Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  float mu0 = 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3;
  mu0 = 1e2*tbarx/mu0;

  // Create Table 3, Huber et al. (2009):
  float h_array[7][6] = {
     {0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0},  // hj0
     {0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573}, // hj1
     {-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0}, // hj2
     {0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0}, // hj3
     {-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0}, // hj4
     {0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0}, // hj5
     {0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264} // hj6
  };

  // Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  float mu1 = 0.0;
  float ctbar = (1.0/tbar) - 1.0;
  // print(paste("ctbar",ctbar))
  // for i in xrange(6):
  for (int i=0; i<6; ++i){ // in 1:6){
    float coef1 = pow(ctbar, i);
    // print(paste("i, coef1", i, coef1))
    float coef2 = 0.0;
    for (int j=0; j<7; ++j){ // in 1:7){
      coef2 = coef2 + h_array[j][i] * pow(rbar - 1.0, j);
    }
    mu1 = mu1 + coef1 * coef2;    
  }
  mu1 = exp( rbar * mu1 );
  // print(paste("mu1",mu1))

  // Calculate mu_bar (Eq. 2, Huber et al., 2009)
  //   assumes mu2 = 1
  float mu_bar = mu0 * mu1;

  // Calculate mu (Eq. 1, Huber et al., 2009)
  float mu = mu_bar * mu_ast;    // Pa s

  return mu; 
}


//-----------------------------------------------------------------------
// Input:    - float, ambient temperature (tc), degrees C
// Return:   float, viscosity of water (mu), Pa s
// Features: Calculates viscosity of water at a given temperature and 
//           pressure.
// Depends:  density_h2o
// Ref:      Vogel ...
//-----------------------------------------------------------------------
inline float calc_viscosity_h2o_vogel(float tc){

	float tk = tc + 272.15;

	float a = -3.7188;
	float b = 578.919;
	float c = 137.546;

	float visc = 1e-3 * exp(a + b/(tk - c));

	return visc;
}


//----------------------------------------------------------------
// Calculates the instantaneous temperature response of the quantum
// yield efficiency based on Bernacchi et al., 2003 PCE (Equation
// and parameter values taken from Appendix B)
//----------------------------------------------------------------
inline float calc_ftemp_kphio(float tc ){
	float ftemp = 0.352 + 0.022 * tc - 3.4e-4 * tc * tc;

	return ftemp;
}


} // phydro

#endif
