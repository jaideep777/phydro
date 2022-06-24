#ifdef USINGRCPP

#include <Rcpp.h>
using namespace Rcpp;

#include "phydro.h"

using namespace phydro;

RCPP_MODULE(phydro_module) {

	function("rphydro_numerical", &rphydro_numerical);
	function("rphydro_analytical", &rphydro_analytical);
	function("rphydro_instantaneous_numerical", &rphydro_instantaneous_numerical);
	function("rphydro_instantaneous_analytical", &rphydro_instantaneous_analytical);

	function("calc_kmm", &calc_kmm);
	function("calc_patm", &calc_patm);
	function("calc_density_h2o", &calc_density_h2o);
	function("calc_viscosity_h2o", &calc_viscosity_h2o);

}

#endif

