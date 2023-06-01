#ifndef PHYDRO_PARAMS_CLASSES_H
#define PHYDRO_PARAMS_CLASSES_H

#include <unordered_map>

#include "temperature_dependencies.h"

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


class ParPlant{
	public:
	double conductivity;
	double psi50;
	double b;
	
	ParPlant(double _conductivity, double _psi50, double _b){
		conductivity = _conductivity;
		psi50 = _psi50;
		b = _b;
	}
};


class ParCost{
	public:
	double alpha;
	double gamma;

	ParCost(double _a, double _g){
		alpha = _a;
		gamma = _g;
	}
};


class ParPhotosynth{
	public:
	double kmm;
	double gammastar;
	double phi0;
	double ca;
	double delta;  // TODO: Replace name with brd / rdark

	FtempVcmaxJmaxMethod ftemp_vj_method;
	FtempRdMethod        ftemp_rd_method;
	FtempBrMethod        ftemp_br_method;

	double Iabs;
	double patm;

	double fT_vcmax;
	double fT_jmax;
	double fT_rd;

	ParPhotosynth(double _tc, double _patm, double _kphio, double _co2, double _ppfd, double _fapar, double _rdark25,
				  FtempVcmaxJmaxMethod _ftemp_vj_method = FV_kumarathunge19, 
				  FtempRdMethod        _ftemp_rd_method = FR_heskel16, 
				  FtempBrMethod        _ftemp_br_method = FB_atkin15){
		
		ftemp_vj_method = _ftemp_vj_method;
		ftemp_rd_method = _ftemp_rd_method;
		ftemp_br_method = _ftemp_br_method;

		fT_vcmax = calc_ftemp_inst_vcmax(_tc, _tc, 25.0, ftemp_vj_method);
		fT_jmax  = calc_ftemp_inst_jmax(_tc, _tc, _tc, 25, ftemp_vj_method);
		fT_rd    = calc_ftemp_inst_rd(_tc, _ftemp_rd_method);

		kmm = calc_kmm(_tc, _patm);
		gammastar = calc_gammastar(_tc, _patm);
		phi0 = _kphio*calc_ftemp_kphio(_tc);
		Iabs = _ppfd * _fapar;
		ca = _co2 * _patm * 1e-6;
		patm = _patm;
		delta = _rdark25 * fT_rd / fT_vcmax;

	}

};

class ParControl{
	public:
	GsMethod             gs_method = GS_IGF;
	ETMethod             et_method = ET_DIFFUSION;
	FtempVcmaxJmaxMethod ftemp_vj_method = FV_kumarathunge19;
	FtempRdMethod        ftemp_rd_method = FR_heskel16; 
	FtempBrMethod        ftemp_br_method = FB_atkin15;
};


} // phydro

#endif
