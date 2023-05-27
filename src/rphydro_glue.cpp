// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//       R Interface
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifdef USINGRCPP

#include <Rcpp.h>
using namespace Rcpp;

#include "phydro.h"

namespace phydro{

inline Rcpp::List PHydroResult_to_List(const phydro::PHydroResult& res){
	return Rcpp::List::create(
	           Named("a") = res.a,
	           Named("e") = res.e,
	           Named("gs") = res.gs,
	           Named("ci") = res.ci,
	           Named("chi") = res.chi,
	           Named("vcmax") = res.vcmax,
	           Named("jmax") = res.jmax,
	           Named("dpsi") = res.dpsi,
	           Named("psi_l") = res.psi_l,
	           Named("profit") = 0,
	           Named("mc") = res.mc,
	           Named("mj") = res.mj,
	           Named("gammastar") = res.gammastar,
	           Named("kmm") = res.kmm,
	           Named("vcmax25") = res.vcmax25,
	           Named("jmax25") = res.jmax25,
	           Named("rd") = res.rd,
	           Named("isVcmaxLimited") = res.isVcmaxLimited,
	           Named("ac") = res.ac,
	           Named("aj") = res.aj
	       );
}

inline Rcpp::List rphydro_analytical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_analytical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

inline Rcpp::List rphydro_instantaneous_analytical(double vcmax25, double jmax25, double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_instantaneous_analytical(vcmax25, jmax25, tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}


#ifndef PHYDRO_ANALYTICAL_ONLY

inline Rcpp::List rphydro_numerical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_numerical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

inline Rcpp::List rphydro_instantaneous_numerical(double vcmax25, double jmax25, double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_instantaneous_numerical(vcmax25, jmax25, tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

#endif


} // phydro namespace


using namespace phydro;

RCPP_MODULE(phydro_module) {
	// Temperature dependencies
	function("calc_kmm", &calc_kmm);
	function("calc_patm", &calc_patm);
	function("calc_gammastar", &calc_gammastar);
	function("calc_density_h2o", &calc_density_h2o);
	function("calc_viscosity_h2o", &calc_viscosity_h2o);
	function("calc_ftemp_kphio", &calc_ftemp_kphio);
	// function("calc_ftemp_inst_vcmax", &calc_ftemp_inst_vcmax);
	// function("calc_ftemp_inst_jmax", &calc_ftemp_inst_jmax);
	// function("calc_ftemp_inst_rd", &calc_ftemp_inst_rd);

	// Phydro core
	function("rphydro_analytical", &rphydro_analytical);
	function("rphydro_instantaneous_analytical", &rphydro_instantaneous_analytical);

#ifndef PHYDRO_ANALYTICAL_ONLY
	function("rphydro_numerical", &rphydro_numerical);
	function("rphydro_instantaneous_numerical", &rphydro_instantaneous_numerical);
#endif

}

#endif

