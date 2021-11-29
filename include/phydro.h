#ifndef PHYDRO_PHYDRO_H
#define PHYDRO_PHYDRO_H


#include "hyd_analytical_solver.h"

namespace phydro{

struct PHydroResult{
	double a;
	double e;
	double gs;
	double ci;
	double chi;
	double vcmax;
	double jmax;
	double dpsi;
	double psi_l;
	double nfnct;
	double niter;
};

inline PHydroResult phydro_analytical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1)){
	
	double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark);
	ParEnv        par_env(tc, pa, vpd);

	auto     bounds = calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost);
	auto   dpsi_opt = pn::zero(bounds.Iabs_bound*0.001, bounds.Iabs_bound*0.999, [&](double dpsi){return dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost).dPdx;}, 1e-6);
	double        x = calc_x_from_dpsi(dpsi_opt.root, psi_soil, par_plant, par_env, par_photosynth, par_cost);  	
	double       gs = calc_gs(dpsi_opt.root, psi_soil, par_plant, par_env);
	double        J = calc_J(gs, x, par_photosynth);	
	double     jmax = calc_jmax_from_J(J, par_photosynth); 
	double    vcmax = (J/4.0)*(x*par_photosynth.ca + par_photosynth.kmm)/(x*par_photosynth.ca + 2*par_photosynth.gammastar);
	double        a = gs*(par_photosynth.ca/par_photosynth.patm*1e6)*(1-x);

	PHydroResult res;
	res.a = a;
	res.e = 1.6*gs*vpd/par_env.patm;
	res.ci = x*par_photosynth.ca;
	res.gs = gs;
	res.chi = x;
	res.vcmax = vcmax;
	res.jmax = jmax;
	res.dpsi = dpsi_opt.root;
	res.psi_l = psi_soil - dpsi_opt.root;
	res.nfnct = dpsi_opt.nfnct;

	return res;

}

} // phydro

#include "hyd_numerical_solver.h"


namespace phydro{

inline PHydroResult phydro_numerical(double tc, double ppfd, double vpd, double co2, double elv, double fapar, double kphio, double psi_soil, double rdark, ParPlant par_plant, ParCost par_cost = ParCost(0.1,1)){
	
	double pa = calc_patm(elv);

	ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark);
	ParEnv        par_env(tc, pa, vpd);

	auto     opt = optimize_midterm_multi(psi_soil, par_cost, par_photosynth, par_plant, par_env);
	double    gs = calc_gs(opt.dpsi, psi_soil, par_plant, par_env);
	auto      aj = calc_assim_light_limited(gs, opt.jmax, par_photosynth); 	
	double vcmax = vcmax_coordinated_numerical(aj.a, aj.ci, par_photosynth);

	PHydroResult res;
	res.a = aj.a;
	res.e = 1.6*gs*vpd/par_env.patm;
	res.ci = aj.ci;
	res.gs = gs;
	res.chi = aj.ci/par_photosynth.ca;
	res.vcmax = vcmax;
	res.jmax = opt.jmax;
	res.dpsi = opt.dpsi;
	res.psi_l = psi_soil - opt.dpsi;

	return res;

}

} // phydro

#endif
