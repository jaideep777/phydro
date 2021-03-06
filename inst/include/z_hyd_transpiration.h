#ifndef PHYDRO_TRANSPIRATION_H
#define PHYDRO_TRANSPIRATION_H

#include "pn_integrator.h"
#include "hyd_params_classes.h"
//#include <unsupported/Eigen/SpecialFunctions>
#include <gsl/gsl_sf_gamma.h>

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
inline double integral_P_analytical(double dpsi, double psi_soil, double psi50, double b){
	double ps = psi_soil/psi50;
	double pl = (psi_soil-dpsi)/psi50;
	double l2 = log(2);
	//double I = -(psi50/b)*pow(l2,-1/b)*b*(Eigen::numext::igammac(1/b, l2*pow(pl,b)) - Eigen::numext::igammac(1/b, l2*pow(ps,b)));
	double I = -(psi50/b)*pow(l2,-1/b)*(gsl_sf_gamma_inc(1/b, l2*pow(pl,b)) - gsl_sf_gamma_inc(1/b, l2*pow(ps,b)));
	return I;
}


inline double integral_P_approx(double dpsi, double psi_soil, double psi50, double b){
	return P(psi_soil-dpsi/2, psi50, b)*dpsi;
}


inline double calc_transpiration(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	double K = scale_conductivity(par_plant.conductivity, par_env);
	double E = K * -integral_P_numerical(dpsi, psi_soil, par_plant.psi50, par_plant.b);
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
	return K/1.6/D*(P(psi_soil, par_plant.psi50, par_plant.b) - Pprime(psi_soil, par_plant.psi50, par_plant.b)*dpsi - Pprimeprime(psi_soil, par_plant.psi50, par_plant.b)*3*dpsi*dpsi/8);
}


inline double calc_gsprime(double dpsi, double psi_soil, ParPlant par_plant, ParEnv par_env){
	return calc_gsprime_analytical(dpsi, psi_soil, par_plant, par_env);
}

} // phydro




#endif


