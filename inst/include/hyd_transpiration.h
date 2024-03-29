#ifndef PHYDRO_TRANSPIRATION_H
#define PHYDRO_TRANSPIRATION_H

#include <stdexcept>

#include "pn_integrator.h"
#include "hyd_params_classes.h"
//#include <unsupported/Eigen/SpecialFunctions>

#ifdef USE_GSL_GAMMA
#include <gsl/gsl_sf_gamma.h>
#else 
#include "incgamma.h"
#endif

namespace phydro{

// Vulnerability curve
inline double P(double psi, double psi50, double b){
	return pow(0.5, pow(psi/psi50, b));
}

inline double Pprime(double psi, double psi50, double b){
	return log(0.5)*P(psi,psi50,b)*b*pow(psi/psi50, b-1)/psi50;
}

inline double Pprimeprime(double psi, double psi50, double b){
	return log(0.5)*b*pow(psi/psi50,b-1)/psi50*Pprime(psi, psi50, b) + log(0.5)*P(psi, psi50, b)/(psi50*psi50)*b*(b-1)*pow(psi/psi50,b-2);
}


// Convert conductivity from m (m3/m2) to mol/m2/s/Mpa
inline double scale_conductivity(double K, ParEnv par_env){
	// Flow rate in m3/m2/s/Pa
	double K2 = K/par_env.viscosity_water;

	// Flow rate in mol/m2/s/Pa
	double mol_h20_per_kg_h20 = 55.5;
	double K3 = K2 * par_env.density_water * mol_h20_per_kg_h20;
	
	// Flow rate in mol/m2/s/Mpa
	double K4 = K3 * 1e6;

	return K4;
}


// integrate vulnerability curve
inline double integral_P_numerical(double dpsi, double psi_soil, double psi50, double b){
	auto p_func = [psi50, b](double psi){
		return P(psi, psi50, b);
	};

	double I = pn::Integrator().integrate(p_func, psi_soil, psi_soil-dpsi);
	return I;
}


// integrate vulnerability curve
// int P(p, p50, b) = -(p/b) * (log2)^(-1/b) * G(1/b, (x/p)^b*log2)  <--- G is unnormalized upper incomplete gamma function (GSL and gammainc impl)
//                  = -(p/b) * (log2)^(-1/b) * G(1/b) * (1 - I((x/p)^b*log2) <--- I is lower incomplete gamma integral (gammad impl)
//                  = -(p/b) * (log2)^(-1/b) * G(1/b) * (- I((pl/p)^b*log2 + I((ps/p)^b*log2) <--- I is lower incomplete gamma integral
//                  = +(p/b) * (log2)^(-1/b) * G(1/b) * (  I((pl/p)^b*log2 - I((ps/p)^b*log2) <--- I is lower incomplete gamma integral
inline double integral_P_analytical(double dpsi, double psi_soil, double psi50, double b){
	double ps = psi_soil/psi50;
	double pl = (psi_soil-dpsi)/psi50;
	double l2 = log(2);
#ifdef USE_GSL_GAMMA
	double I = -(psi50/b)*pow(l2,-1/b)*(gsl_sf_gamma_inc(1/b, l2*pow(pl,b)) - gsl_sf_gamma_inc(1/b, l2*pow(ps,b)));
#else
	double I = -(psi50/b)*pow(l2,-1/b)*(gammainc(1/b, l2*pow(pl,b)) - gammainc(1/b, l2*pow(ps,b)));
#endif
	return I;
}


inline double integral_P_approx(double dpsi, double psi_soil, double psi50, double b){
	return -P(psi_soil-dpsi/2, psi50, b)*dpsi;
}

inline double integral_P_approx2(double dpsi, double psi_soil, double psi50, double b){
	return -(P(psi_soil, psi50, b)+P(psi_soil-dpsi, psi50, b))/2 * dpsi;
}


inline double integral_P(double dpsi, double psi_soil, ParPlant par_plant){
	if      (par_plant.gs_method == GS_QNG)  return integral_P_numerical( dpsi, psi_soil, par_plant.psi50, par_plant.b);
	else if (par_plant.gs_method == GS_IGF)  return integral_P_analytical(dpsi, psi_soil, par_plant.psi50, par_plant.b);
	else if (par_plant.gs_method == GS_APX)  return integral_P_approx(    dpsi, psi_soil, par_plant.psi50, par_plant.b);
	else if (par_plant.gs_method == GS_APX2) return integral_P_approx2(   dpsi, psi_soil, par_plant.psi50, par_plant.b);
	else throw std::runtime_error("Unsupported gs_method specified");
}


inline double calc_transpiration(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double K = scale_conductivity(par_plant.conductivity, par_env);
	double E = K * -integral_P(dpsi, psi_soil, par_plant);
	return E;
}


// Calculates regulated stomatal conducatnce given the leaf water potential, 
// plant hydraulic traits, and the environment.
inline double calc_gs(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double E = calc_transpiration(dpsi, psi_soil, par_plant, par_env);
	double D = (par_env.vpd/par_env.patm);
	double gs = E/1.6/D; 
	return gs;
}


inline double calc_gsprime_analytical(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double K = scale_conductivity(par_plant.conductivity, par_env);
	double D = (par_env.vpd/par_env.patm);
	return K/1.6/D*P(psi_soil-dpsi, par_plant.psi50, par_plant.b);
}


inline double calc_gsprime_approx(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double K = scale_conductivity(par_plant.conductivity, par_env);
	double D = (par_env.vpd/par_env.patm);
	return K/1.6/D*(P(psi_soil-dpsi/2, par_plant.psi50, par_plant.b) - Pprime(psi_soil-dpsi/2, par_plant.psi50, par_plant.b)*dpsi/2);
}


inline double calc_gsprime_approx2(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double K = scale_conductivity(par_plant.conductivity, par_env);
	double D = (par_env.vpd/par_env.patm);
	return K/1.6/D* (  (P(psi_soil, par_plant.psi50, par_plant.b)+P(psi_soil-dpsi, par_plant.psi50, par_plant.b))/2 
	                  - Pprime(psi_soil-dpsi, par_plant.psi50, par_plant.b)*dpsi/2 );
}


inline double calc_gsprime(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	if      (par_plant.gs_method == GS_APX)  return calc_gsprime_approx(    dpsi, psi_soil, par_plant, par_env);
	else if (par_plant.gs_method == GS_APX2) return calc_gsprime_approx2(   dpsi, psi_soil, par_plant, par_env);
	else                                     return calc_gsprime_analytical(dpsi, psi_soil, par_plant, par_env);

}

} // phydro




#endif


