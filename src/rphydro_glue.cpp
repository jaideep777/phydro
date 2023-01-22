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
	           Named("nfnct") = res.nfnct,
	           Named("niter") = res.niter,
	           Named("chi_jmax_lim") = 0,
	           Named("profit") = 0,
	           Named("mc") = res.mc,
	           Named("mj") = res.mj,
	           Named("gammastar") = res.gammastar,
	           Named("kmm") = res.kmm,
	           Named("vcmax25") = res.vcmax25
	       );
}

inline Rcpp::List rphydro_analytical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_analytical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

inline Rcpp::List rphydro_instantaneous_analytical(double vcmax, double jmax, double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_instantaneous_analytical(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}


#ifndef PHYDRO_ANALYTICAL_ONLY

inline Rcpp::List rphydro_numerical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_numerical(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

inline Rcpp::List rphydro_instantaneous_numerical(double vcmax, double jmax, double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, Rcpp::List par_plant, Rcpp::List par_cost){
	ParPlant par_plant_cpp(par_plant["conductivity"], par_plant["psi50"], par_plant["b"]);
	ParCost  par_cost_cpp(par_cost["alpha"], par_cost["gamma"]);
	return PHydroResult_to_List(phydro_instantaneous_numerical(vcmax, jmax, tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark, par_plant_cpp, par_cost_cpp));
}

#endif


} // phydro namespace


using namespace phydro;

RCPP_MODULE(phydro_module) {

	function("calc_kmm", &calc_kmm);
	function("calc_patm", &calc_patm);
	function("calc_density_h2o", &calc_density_h2o);
	function("calc_viscosity_h2o", &calc_viscosity_h2o);

	function("rphydro_analytical", &rphydro_analytical);
	function("rphydro_instantaneous_analytical", &rphydro_instantaneous_analytical);

#ifndef PHYDRO_ANALYTICAL_ONLY
	function("rphydro_numerical", &rphydro_numerical);
	function("rphydro_instantaneous_numerical", &rphydro_instantaneous_numerical);
#endif

}

#endif

