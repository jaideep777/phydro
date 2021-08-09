#include <iostream>
#include "hyd_transpiration.h"
#include <chrono>
using namespace std;

int main(){

	double tc=25;
	double p = phydro::calc_patm(0);
	double vpd = 1000;
	
	double g;
	
	auto t1 = std::chrono::high_resolution_clock::now();
	for (int i=0; i<1000; ++i){
		double psi_s = -6.0 + i*(6.0)/(1000-1);
		g = phydro::calc_gs(1, psi_s, phydro::ParPlant(3e-17, -2, 2), phydro::ParEnv(tc, p, vpd));
	}
	auto t2 = std::chrono::high_resolution_clock::now();

	cout << g << endl;
	cout << "Time required: " << (std::chrono::duration<double, std::milli> (t2 - t1)).count() << " ms\n";


	if (abs(g-0.1116382) < 1e-5) return 0;
	else return 1;
}


