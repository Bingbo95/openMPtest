#pragma once
#include <vector>
#include <string>

namespace Hydrate_Related_Parameters
{
	struct InhibitorParameters {
		int inhibitor_flag;
		double Max_TShift;
		double Y_atMax_TShift;
		double InhibMW;
		double InhibDens;
		double InhibEnthSol;
		double InhibCpCoeff;
	};
	extern int NCom;
	extern int hydrN;
	extern int N_ThC;
	extern int N_SpH;
	extern int N_Rho;
	extern int EquationOption;

	extern double moleF;
	extern std::vector<double> An, Bn, Dn;
	extern std::string nameG;
	extern InhibitorParameters inhibp;
}