#pragma once
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>

namespace Hydrate_Related_Parameters
{
	/*struct InhibitorParameters {
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
	extern InhibitorParameters inhibp;*/
	//struct InhibitorParameters {
	//	int inhibitor_flag;
	//	std::string Inhibitor_Name;
	//	double Max_TShift;
	//	double Y_atMax_TShift;
	//	double InhibMW;
	//	double InhibDens;
	//	double InhibEnthSol;
	//	double InhibCpCoeff[3];
	//};
	//struct HydrateGas {
	//	std::string NameG;
	//	int HydrNum;
	//	double MoleFracG;
	//};
	//extern std::vector<HydrateGas> hydrGasP;
	//extern int NumHydrCom;
	//extern int N_ThC;
	//extern int N_SpH;
	//extern int N_Rho;
	//extern int EquationOption;
	//extern int SandFlow_flag;
	//extern std::vector<double> p_ThC, p_SpH, p_Rho;
	//extern InhibitorParameters inhibp;
	struct InhibitorParameters {
		int Inhibitor_flag;
		std::string Inhibitor_Name;
		double Max_TShift;
		double Y_atMax_TShift;
		double InhibMW;
		double InhibDens;
		double InhibEnthSol;
		double InhibCpCoeff[3];
	};
	struct HydrateGas {
		std::string NameG;
		double HydrNum;
		double MoleFracG;
		double HydrMW;//
		double MassFracGas;//
		double MassFracH2O;//
	};
	//extern std::vector<std::vector<double>> xhydrsolid;
	extern std::vector<HydrateGas> hydrGasP;
	extern int NumHydrCom;
	extern int N_ThC;
	extern int N_SpH;
	extern int N_Rho;
	extern int EquationOption;
	extern int eqTableLength;
	extern int eqSectionNum;
	//extern int SandFlow_flag;//????????????--> MODEL block
	extern double P_quad_hyd;//
	extern double MW_Hydrate;// All of the gas hydrate, not only methane hydrate.
	extern double InitialHydrateMass;//???????
	extern double RefDepartureEnth;//???????
	extern std::vector<double> p_ThC, p_SpH, p_Rho;
	extern std::vector<double> tempHi, tempLow;
	extern std::vector<double> eqPressure, eqTemper;
	extern std::vector<std::vector<double>> acof;
	extern InhibitorParameters inhibp;

	// Bingham Fluid?????????????
	//extern bool pressure_gradiant_eff;
	//extern double EffectDP_fac1;
	//extern double EffectDP_fac2;

	extern const double MW_CH4;//
	extern const double Tc_qup_hyd;//
	extern const double Tk_qup_hyd;//
	extern std::string gas_viscosity_option;		// Selection to compute the viscosity of gas.
	extern std::string gas_cubic_EOS;				// The cubic equation of state.
	extern bool Lee_Kessler_method;		// Flag indicating whether the Lee-Kessker method is used to compute the departure enthalpy.
	extern bool noInhib;
}