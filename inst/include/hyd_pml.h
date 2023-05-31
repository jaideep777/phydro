#ifndef PHYDRO_PML_H
#define PHYDRO_PML_H

#include <stdexcept>
#include <cmath>

namespace phydro{


// Calculate saturation vapour pressure at given temperature and pressure [Pa]
// This function was adopted from the R package 'plantecophys'
// Duursma (2015) https://doi.org/10/bkmj.
// patm    Atmospheric pressure [Pa]
// TdegC   Temperature [degC]
inline double calc_esat(double TdegC, double patm = 101325) {
    // Pa in kPa
    double a = 611.21;
    double b = 17.502;
    double c = 240.97;
    double f = 1.0007 + 3.46e-8 * patm;
    double esatval = f * a * (std::exp(b * TdegC / (c + TdegC)));
    return esatval;
}


// Aerodynamic conductance [m s-1]
// To convert to mol m-2 s-1, see this: https://rdrr.io/cran/bigleaf/man/ms.to.mol.html (but not convincing)
// Refs: 
//    Eq 13 in Leuning et al (2008). https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2007WR006562
//    Eq 7 in Zhang et al (2008): https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017JD027025
//    Box 4 in https://www.fao.org/3/x0490e/x0490e06.htm 
// h_canopy    canopy height [m]
// v_wind      wind speed [m s-1]
// z_measurement    height of wind speed measurement [m]
inline double calc_g_aero(double h_canopy, double v_wind, double z_measurement){
	double k_karman = 0.41;        // von Karman's constant [-]
	double d = h_canopy*2.0/3.0;   // zero-plane displacement height [m]
	double z_om = 0.123*h_canopy;  // roughness lenghts governing transfer of water and momentum [m]
	double z_ov = 0.1*z_om;
	
	double g_aero = (k_karman*k_karman * v_wind) / ( log((z_measurement-d)/z_om)*log((z_measurement-d)/z_ov) );
	return g_aero;
}


// calculate density of dry (or moist) air [kg m-3]
// tc_air    Air temperature [degC]
// patm      Atmospheric pressure [Pa]
// vpd       Vapour pressure deficit [Pa]
inline double calc_density_air(double tc_air, double patm, double vpd, bool moist = true){
	double tk = tc_air+273.16;
	double R = 287.052874; // Specific universal gas const for dry air [J kg-1 K-1]

	// for dry air
	if (!moist) return patm / R / tk;

	// for moist air
	double vp = calc_esat(tc_air, patm) - vpd;    // vapour pressure
	double rv = 0.622 * vp/(patm - vp);   // https://glossary.ametsoc.org/wiki/Mixing_ratio
	double tv = tk * (1 + rv/0.622)/(1 + rv);  // virtual temperature. https://glossary.ametsoc.org/wiki/Virtual_temperature
	return patm / R / tv;
}


// calculate enthalpy of vapourization [J kg-1]
// tc   temperature [degC]
// Ref:
//   Eq. 8, Henderson-Sellers (1984) https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.49711046626
inline double calc_enthalpy_vap(double tc){
	double tk = tc + 273.15;
	double a = tk/(tk - 33.91);
	return 1.91846e6*a*a;
}


// calc specific heat capacity of *saturated* moist air [J kg-1 K-1]
// Ref:
//     Eq. 47, Tsilingiris (2008)
// tc   Temperature [degC]
inline double calc_cp_moist_air(double tc){
	// Adjust temperature to avoid numerical blow-up
	// Adopted temperature adjustment from SPLASH, Python version
    double my_tc = std::min(std::max(tc, 0.0), 100.0);

    // Calculate the specific heat capacity of water, J/kg/K
    // 
    // cp = 1.0045714270 +
    //      2.050632750e-3 * my_tc -
    //      1.631537093e-4 * my_tc * my_tc +
    //      6.212300300e-6 * my_tc * my_tc * my_tc -
    //      8.830478888e-8 * my_tc * my_tc * my_tc * my_tc +
    //      5.071307038e-10 * my_tc * my_tc * my_tc * my_tc * my_tc;
	double cp = (1.0045714270 +
		 my_tc * ( 2.050632750e-3 +
		 my_tc * (-1.631537093e-4 +
		 my_tc * ( 6.212300300e-6 -
		 my_tc * ( 8.830478888e-8 -
		 my_tc *   5.071307038e-10))))) * 1e3;

	return cp;
} 

// Calculate Psychrometric constant [Pa/K]
inline double calc_psychro(double tc, double patm) {
    // Constants
    const double Ma = 0.02896;  // Molecular weight dry air, kg/mol
    const double Mv = 0.018016; // Molecular weight of water vapor, kg/mol

	// calculate specific heat capacity of moist air
	double cp = calc_cp_moist_air(tc);

    // Calculate latent heat of vaporization, J/kg
    double lv = calc_enthalpy_vap(tc);

    // Calculate psychrometric constant, Pa/K
    // Eq. 8, Allen et al. (1998)
	// double psychro = cp * kMa * patm / (kMv * lv);   // J kg-1 K-1 * Pa /  J kg-1 = Pa K-1
	double psychro = cp * patm / ((Mv/Ma) * lv); // JJ: As per https://www.fao.org/3/x0490e/x0490e07.htm#psychrometric%20constant%20(g)

    return psychro;
}


// Calculates the slope of the sat pressure ~ temp curve, Pa/K
// Ref:      Eq. 13, Allen et al. (1998)
// tc        air temperature, degrees C
// output:   slope of the sat pressure temp curve, Pa/K
inline double calc_sat_slope(double tc){
	return 17.269*237.3*610.78 * exp(tc*17.269/(tc + 237.3)) / ((tc + 237.3)*(tc + 237.3));
}



inline double calc_transpiration_pml(double Iabs, double gs, double ga, double tc, double patm, double vpd){
	double rho = calc_density_air(tc, patm, vpd, true);
	double cp = calc_cp_moist_air(tc);
	double gamma = calc_psychro(tc, patm);
	double epsilon = calc_sat_slope(tc) / gamma;

	double trans = (epsilon*Iabs + (rho*cp/gamma)*ga*vpd) / (epsilon + 1 + ga/gs);
	return trans;
}


inline double calc_gs_pml(double Iabs, double Q, double ga, double tc, double patm, double vpd){
	double rho = calc_density_air(tc, patm, vpd, true);
	double cp = calc_cp_moist_air(tc);
	double gamma = calc_psychro(tc, patm);
	double epsilon = calc_sat_slope(tc) / gamma;

	double den = epsilon*Iabs + (rho*cp/gamma)*ga*vpd - (1+epsilon)*Q; 
	double gs = ga * Q / den;
	return gs;
}


inline double calc_dE_dgs_pml(double Iabs, double gs, double ga, double tc, double patm, double vpd){
	double rho = calc_density_air(tc, patm, vpd, true);
	double cp = calc_cp_moist_air(tc);
	double gamma = calc_psychro(tc, patm);
	double epsilon = calc_sat_slope(tc) / gamma;

	double num = ga * (epsilon*Iabs + (rho*cp/gamma)*ga*vpd);
	double den = epsilon*gs + gs + ga;

	return num/den/den;
}


inline double calc_dE_dgs_pml_num(double Iabs, double gs, double ga, double tc, double patm, double vpd){
	double E = calc_transpiration_pml(Iabs, gs, ga, tc, patm, vpd);
	double E_plus = calc_transpiration_pml(Iabs, gs+1e-6, ga, tc, patm, vpd);

	return (E_plus-E)/1e-6;
}


} // namespace phydro


#endif


