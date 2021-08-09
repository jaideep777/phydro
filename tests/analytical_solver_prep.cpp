#include <iostream>
#include <iomanip>
#include <fstream>
#include "hyd_analytical_solver.h"
#include <chrono>
#include <vector>
using namespace std;


double err(double x, double ref){
	return std::min(abs(x-ref), abs(x/ref-1));
}

int check(double x, double ref, double err=1e-5){
	//cout << "err: " << abs(x-ref) << " " << abs(x/ref-1)<< "\n";
	//cout << "comp: " << x << " " << ref << "|" << abs(x-ref) << "\n";
	if (abs(x-ref) < err || abs(x/ref-1) < err) return 0;
	else return 1;
}


vector<double> seq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(start + double(i)*(end-start)/(length-1));
	return x;
}

vector<double> lseq(double start, double end, int length){
	vector<double> x;
	for (int i=0; i<length; ++i) x.push_back(exp(log(start) + double(i)*(log(end)-log(start))/(length-1)));
	return x;
}

int main(){

	double psi_soil = -2;

	double kphio = 0.087;        // quantum yield efficiency
	// c_molmass <- 12.0107 # molar mass, g / mol

	// Define environmental conditions
	double tc = 24.591837;             // temperature, deg C
	double ppfd = 300;          // umol/m2/s
	double vpd  = 1000;         // Pa
	double co2  = 400;          // ppm
	double elv  = 0;            // m.a.s.l.
	double fapar = 0.7;         // fractioni
	double rdark = 0.0;

	double pa = phydro::calc_patm(elv);

	phydro::ParCost       par_cost(0.1, 0.1);
	phydro::ParPlant      par_plant(3e-17, -2, 2);
	phydro::ParPhotosynth par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark);
	phydro::ParEnv        par_env(tc, pa, vpd);

	par_plant.gs_method = phydro::GS_IGF;

    // if (check(calc_Aj_max(0.1, 0.75, par_photosynth), 14.86594) == 1) return 1;
	
	//if (check(calc_jmax_from_Ajmax(11, par_photosynth), 91.18512) == 1) return 1;
	
	//if (check(calc_djmax_dAjmax(11, par_photosynth), 35.60196, 1e-5) == 1) return 1;
	
	cout << setprecision(10) << "dAjmax_dchi = " << calc_dAjmax_dchi(0.1, 0.75, par_photosynth) << "\n";
	//if (fabs(calc_dAjmax_dchi(0.1, 0.75, par_photosynth) - (-67.00399)) > 5e-5) return 1;
	
	cout << setprecision(10) << "x from dpsi = " << calc_x_from_dpsi(1.23, psi_soil, par_plant, par_env, par_photosynth, par_cost) << "\n";
	//if (fabs(calc_x_from_dpsi(1.23, psi_soil, par_plant, par_env, par_photosynth, par_cost) - 0.736516) > 1e-5) return 1;
	
	cout << setprecision(10) << "delta from dpsi = " << calc_delta_from_dpsi(1.23, psi_soil, par_plant, par_env, par_photosynth, par_cost) << "\n";
	//if (fabs(calc_delta_from_dpsi(1.23, psi_soil, par_plant, par_env, par_photosynth, par_cost) - 9.038855) > 1e-5) return 1;
	
	auto dfdx = dFdx(1.23, psi_soil, par_plant, par_env, par_photosynth, par_cost);
	cout << "dPdx = " << dfdx.dPdx << "\najmax = " << dfdx.ajmax << "\ndjmax_dajmax = " << dfdx.djmax_dajmax << "\ndajmax_dchi = " << dfdx.dajmax_dchi << "\n";
	//if (check(dfdx.dPdx, 4.094848) == 1) return 1;
	//if (check(dfdx.ajmax, 7.146608) == 1) return 1;
	//if (check(dfdx.djmax_dajmax, 7.193898) == 1) return 1;
	//if (check(dfdx.dajmax_dchi, -30.87851) == 1) return 1;

	auto bounds = calc_dpsi_bound(psi_soil, par_plant, par_env, par_photosynth, par_cost);	
	cout << "exact = " << bounds.exact << "\napprox_O2 = " << bounds.approx_O2 << "\nIabs_bound = " << bounds.Iabs_bound << "\n";
	//if (check(bounds.exact, 3.103414) == 1) return 1;
	//if (check(bounds.approx_O2, 1.706554) == 1) return 1;
	//if (check(bounds.Iabs_bound, 1.711369) == 1) return 1;

	auto dpsi_opt = pn::zero(bounds.Iabs_bound*0.001, bounds.Iabs_bound*0.999, [&](double dpsi){return dFdx(dpsi, psi_soil, par_plant, par_env, par_photosynth, par_cost).dPdx;}, 1e-6);
	cout << "optimal dpsi = " << dpsi_opt.root << "\n";
	cout << "func evals   = " << dpsi_opt.nfnct << "\n";
	if (check(dpsi_opt.root, 1.119076) == 1) return 1;

	return 0;
}





