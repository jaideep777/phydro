#include "phydro.h"
using namespace std;
using namespace phydro;

int main(){
	
	double tc = 20;
	double vpd = 810.6;
	double co2 = 400;
	double ppfd = 1200;
	double pa = calc_patm(0);
	double kphio = 0.087;
	double fapar = 0.99;
	double rdark = 0.02;
	double psi_soil = -0.4137931;
	
	ParCost        par_cost(0.118514, 1.227068);
	ParPlant       par_plant(7.457324e-17, -1.039539, 1);
	ParPhotosynth  par_photosynth(tc, pa, kphio, co2, ppfd, fapar, rdark);
	ParEnv         par_env(tc, pa, vpd, ppfd/2);

	double jmax = 117.0184518;
	double vcmax = 55.6279401;

	PHydro_Profit_Inst P(vcmax, jmax, psi_soil, par_cost, par_photosynth, par_plant, par_env);
	for (int i=0; i<50; ++i){
		double dpsi = 2.0/49.0*i;
		VectorXd x(1);
		x << dpsi;
		cout << dpsi << " " << P.value(x) << "\n";
	}

	double dpsi_opt_num = optimize_shortterm_multi(vcmax, jmax, psi_soil, par_cost, par_photosynth, par_plant, par_env);
	cout << "Opt = " << dpsi_opt_num << "\n";

	auto dpsi_opt = pn::zero(0, 20, [&](double dpsi){return calc_dP_ddpsi(dpsi, vcmax, jmax, psi_soil, par_plant, par_env, par_photosynth, par_cost);}, 1e-6);
	cout << "Opt analytical = " << dpsi_opt.root << "\n";

	return 0;
}

