#ifndef PHYDRO_PHOTOSYNTHESIS_H
#define PHYDRO_PHOTOSYNTHESIS_H

#include "hyd_params_classes.h"

namespace phydro{

inline double QUADM(double A, double B, double C){
	return (-B - sqrt(B*B - 4*A*C)) / (2*A);
}

inline double QUADP(double A, double B, double C){
	return (-B + sqrt(B*B - 4*A*C)) / (2*A);
}

struct ACi{
	double a;
	double ci;
};


inline ACi calc_assim_rubisco_limited(double _gs, double vcmax, ParPhotosynth par_photosynth){
	double ca = par_photosynth.ca;            // ca is in Pa
	double gs = _gs * 1e6/par_photosynth.patm; // convert to umol/m2/s/Pa

	double d = par_photosynth.delta;

	double A = -1.0 * gs;
	double B = gs * ca - gs * par_photosynth.kmm - vcmax*(1-d);
	double C = gs * ca * par_photosynth.kmm + vcmax * (par_photosynth.gammastar + par_photosynth.kmm*d);

	ACi res;
	res.ci = QUADM(A,B,C);
	res.a  = gs*(ca-res.ci);

	return res;

}


inline ACi calc_assim_light_limited(double _gs, double jmax, ParPhotosynth par_photosynth){
	double ca = par_photosynth.ca;             // ca is in Pa
	double gs = _gs * 1e6/par_photosynth.patm;  // convert to umol/m2/s/Pa

	gs += 1e-12;

	double d = par_photosynth.delta;

	double phi0iabs = par_photosynth.phi0 * par_photosynth.Iabs;
	double jj = 4*phi0iabs/jmax;
	double jlim = phi0iabs / sqrt(1+ jj*jj);

	double A = -1.0 * gs;
	double B = gs * ca - gs * 2 * par_photosynth.gammastar - jlim*(1-d);
	double C = gs * ca * 2*par_photosynth.gammastar + jlim * (par_photosynth.gammastar + d*par_photosynth.kmm);

	ACi res;
	res.ci = QUADM(A,B,C);
	res.a  = gs*(ca-res.ci);

	return res;

}

inline ACi calc_assimilation_limiting(double vcmax, double jmax, double gs, ParPhotosynth par_photosynth){
	auto Ac = calc_assim_rubisco_limited(gs, vcmax, par_photosynth);
	auto Aj = calc_assim_light_limited(gs, jmax, par_photosynth);

	if (Ac.ci > Aj.ci ) return Ac; 
	else				return Aj;
}

} // phydro

#endif
