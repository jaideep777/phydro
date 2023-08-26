#ifndef PHYDRO_PARAMS_CLASSES_H
#define PHYDRO_PARAMS_CLASSES_H

namespace phydro{

class ParPlant{
	public:
	double conductivity;
	double psi50;
	double b;

	double tchome = 25;

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


} // phydro

#endif
