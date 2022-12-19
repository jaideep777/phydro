#ifndef PHYDRO_PARAMS_CLASSES_H
#define PHYDRO_PARAMS_CLASSES_H

#include "temperature_dependencies.h"

namespace phydro{

class ParEnv{
	public:
	double tc;
	double patm;
	double vpd;
	double viscosity_water;
	double density_water;

	ParEnv(double _tc, double _patm, double _vpd){
		tc = _tc;
		vpd = _vpd;
		patm = _patm;
		viscosity_water = calc_viscosity_h2o(tc, patm);
		density_water = calc_density_h2o(tc, patm);
	}
};


enum GsMethod{GS_IGF, GS_QNG, GS_APX, GS_APX2};

class ParPlant{
	public:
	double conductivity;
	double psi50;
	double b;
	
	GsMethod gs_method = GS_IGF;

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

	double Iabs;
	double patm;

	ParPhotosynth(double _tc, double _patm, double _kphio, double _co2, double _ppfd, double _fapar, double _rdark){
		kmm = calc_kmm(_tc, _patm);
		gammastar = calc_gammastar(_tc, _patm);
		phi0 = _kphio*calc_ftemp_kphio(_tc);
	   	Iabs = _ppfd * _fapar;
		ca = _co2 * _patm * 1e-6;
		patm = _patm;
		delta = _rdark;
	}

};

} // phydro

#endif
