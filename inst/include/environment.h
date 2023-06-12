#ifndef PHYDRO_ENVIRONMENT_H
#define PHYDRO_ENVIRONMENT_H

#include "temperature_dependencies_physical.h"

namespace phydro{

enum GsMethod{GS_IGF, GS_QNG, GS_APX, GS_APX2};
enum ETMethod{ET_DIFFUSION, ET_PM};

class ParEnv{
	public:
	double tc;     // Temperature [degC]
	double patm;   // Atmospheric pressure [Pa]
	double vpd;    // VPD [Pa]
	double Rn;     // net radiation [w m-2]
	double v_wind; // Wind speed [m s-1]

	double viscosity_water;  // [Pa s]
	double density_water;    // [kg m-3]

	GsMethod gs_method = GS_IGF;
	ETMethod et_method = ET_DIFFUSION;

	ParEnv(double _tc, double _patm, double _vpd, double _Rn){
		tc = _tc;
		vpd = _vpd;
		patm = _patm;
		Rn = _Rn;
		v_wind = 3; // global average value
		viscosity_water = calc_viscosity_h2o(tc, patm);
		density_water = calc_density_h2o(tc, patm);
	}
};


} // namespace phydro

#endif


