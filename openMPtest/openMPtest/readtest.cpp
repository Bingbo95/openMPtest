#include <iostream>;
#include <fstream>;
#include <string>;		// for getline function
#include <vector>;
#include <algorithm>;
#include <iomanip>; // for setw and precision
#include <boost/algorithm/string.hpp>
#include "readtest.h"


using namespace std;

ifstream MainInputFile;
ofstream LogFile;
std::string LineText;
int Line_Num;

std::string readLineText(int length = 80);
std::string readLineText(int lineLength)
{
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//if (myid == iMaster) getline(MainInputFile, LineText);
	getline(MainInputFile, LineText);
	//pc->BroadcastStr(LineText, iMaster);
	int nl = LineText.length() - 1;
	if (nl >= 0)
		if (int(LineText.at(nl)) == 13)  LineText.replace(nl, 1, 1, ' ');
	nl++;
	Line_Num += 1;
	if (nl > lineLength) { return LineText; }
	else
	{
		LineText.append(lineLength - nl, ' ');
		return LineText;
	}
}

int readSubStr(std::string& strSub, int& nIdx, std::vector<std::string> AryDescribles)
{
	//		using namespace General_Control_Parameters;
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	//		using namespace Diffusion_Parameters;
	//cout << AryDescribles[1] << endl;
	bool hasKeyWord = true;
	int orinIdx;
	if (nIdx < -2) {
		hasKeyWord = false;
		orinIdx = -nIdx - 10;
		nIdx = -1;
	}
	LineText = readLineText();	// 从输入文件中读取数据。
	if (!LineText.empty())	// 如果LineText非空，则第一句用来将非空字符前的空格去掉，第二局是将其后面的空格去掉。
	{
		LineText.erase(0, LineText.find_first_not_of(" "));
		LineText.erase(LineText.find_last_not_of(" ") + 1);
	}

	if (LineText.length() == 0 || LineText.compare("/") == 0) return -1;
	//cout << LineText.length() << endl;
	//cout << LineText.substr(0, 2) << endl;
	if (LineText.length() < 2 || LineText.substr(0, 2).compare("//") == 0) return 1;

	std::string strTmp = LineText;
	transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);

	for (int i = 0; i < AryDescribles.size(); i++)
	{
		//cout << (strTmp.find(AryDescribles[i]) == -1) << endl;

		if (strTmp.find(AryDescribles[i]) == -1) continue;
		nIdx = i;
		break;
	}
	if (nIdx == -1) {
		if (hasKeyWord) {
			//if (myid == iMaster)
			//	LogFile << "Warning: No keyword was found at text line(" << Line_Num << "): '" << LineText << "'. It was negleted!" << std::endl;
			return 1;
		}
		else {
			nIdx = orinIdx;
			strSub = strTmp;
			if (strSub == "ENDTABLE") return 1;
			return 0;
		}
	}

	int nSub = strTmp.rfind(AryDescribles[nIdx]);
	strTmp = LineText.substr(nSub);
	nSub = strTmp.rfind(':');

	if (nSub < 0) {
		LogFile << "Warning: No data indicator ':' was found at text line(" << Line_Num << "): '" << LineText << "'.\n";
		nSub = AryDescribles[nIdx].length() + 2;
		LogFile << "         The input data is : " << strTmp.substr(nSub + 1) << ", pleease make sure is it correct!" << std::endl;
	}

	strSub = strTmp.substr(nSub + 1);
	if (!strSub.empty())
	{
		strSub.erase(0, strSub.find_first_not_of(" "));
		nSub = strSub.find("//");
		if (nSub != -1) strSub.erase(nSub);
		strSub.erase(strSub.find_last_not_of(" ") + 1);
	}

	return 0;
}

// Test of SOLVE module
void readSOLVE()
{
	//using namespace Fluid_Proterties;
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Flow;

	//if (myid == iMaster) {
	//	SOLVE_read = true;
	LogFile << "--->>> Start reading SOLVER INFO" << std::endl;
	//}

	int AMGCL_Solver = 0, Max_NumIterations = 0;
	double MS_convergence_crit = 0.0, AIM_satrg = 0.01;
	char solverSelection = {};
	std::string PreConditioner = "";
	std::string MatrixSolver = "";

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("SLibrary");
	AryDescribles.push_back("Pre_C");
	AryDescribles.push_back("LSolver");
	AryDescribles.push_back("Maxit");
	AryDescribles.push_back("c_crit");
	AryDescribles.push_back("AIM_sa");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		switch (nIdx)
		{
		case 0:
		{
			char szchar = 'A';
			if (strSub.length() > 0) szchar = strSub[0];
			if (szchar == 'P' || szchar == 'p') szchar = 'P';
			else if (szchar == 'T' || szchar == 't') szchar = 'T';
			else if (szchar == 'F' || szchar == 'f') szchar = 'F';
			else szchar = 'A';
			solverSelection = szchar;
		}
		break;
		case 1:
		{
			if (strSub.length() > 0) PreConditioner = strSub;
		}
		break;
		case 2:
		{
			if (strSub.length() > 0) MatrixSolver = strSub;
		}
		break;
		case 3:
		{
			if (strSub.length() > 0) Max_NumIterations = atoi(strSub.c_str());

		}
		break;
		case 4:
		{
			if (strSub.length() > 0) MS_convergence_crit = atof(strSub.c_str());;
		}
		break;
		case 5:
		{
			if (strSub.length() > 0) AIM_satrg = atof(strSub.c_str());;
		}
		break;
		}

	}

	LogFile << "SLibrary = " << solverSelection << endl;
	LogFile << "Pre_C = " << PreConditioner << endl;
	LogFile << "LSolver = " << MatrixSolver << endl;
	LogFile << "Maxit = " << Max_NumIterations << endl;
	LogFile << "SLibrary = " << MS_convergence_crit << endl;
	LogFile << "SLibrary = " << AIM_satrg << endl;
	//if (Solver_Parameters::Max_NumIterations < 30) Solver_Parameters::Max_NumIterations = 300;
	//if (Solver_Parameters::MS_convergence_crit <= 0.0) Solver_Parameters::MS_convergence_crit = 1.0e-6;
	//if (Solver_Parameters::solverSelection == 'A') {
	//	if (Solver_Parameters::PreConditioner == "CPR") Solver_Parameters::AMGCL_Solver = 1;
	//	else if (Solver_Parameters::PreConditioner == "CPR_DRS")  Solver_Parameters::AMGCL_Solver = 2;
	//	else if (Solver_Parameters::PreConditioner == "AMG")  Solver_Parameters::AMGCL_Solver = 5;
	//	else if (Solver_Parameters::PreConditioner == "ILU")  Solver_Parameters::AMGCL_Solver = 8;
	//	else if (Solver_Parameters::PreConditioner == "SCHUR")  Solver_Parameters::AMGCL_Solver = 9;
	//	else Solver_Parameters::AMGCL_Solver = 1;
	//}
	//if (General_Control_Parameters::AIM_satrg <= 0.0 || General_Control_Parameters::AIM_satrg > 0.9)
	//	General_Control_Parameters::AIM_satrg = 0.01;

	//if (myid == iMaster) {
	//	LogFile << std::endl << "--->>> End reading SOLVER INFO" << std::endl << std::endl;
	//	LogFile << "===============================================================" << std::endl << std::endl;
	//}
}

// Test of TIMES module
void readTIMES()
{
	//using namespace General_Control_Parameters;
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;

	int  Max_NumPrintTimes = 0;
	double  PrintTimeIncrement = 0.0;
	int timeUnit = 0;

	std::vector<std::string> AryDescribles;
	std::vector<double> OutputTimes;
	AryDescribles.push_back("OutputTime");
	AryDescribles.push_back("MaxPrintTimes");
	AryDescribles.push_back("OPTimeIncrement");
	AryDescribles.push_back("TimeUnit");
	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		switch (nIdx)
		{
		case 0:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				OutputTimes.push_back(atof(p));
				p = strtok(NULL, ",");
			}
		}
		break;
		case 1:
		{
			Max_NumPrintTimes = atoi(strSub.c_str());
		}
		break;
		case 2:
		{
			PrintTimeIncrement = atof(strSub.c_str());
		}
		break;
		case 3:
		{
			timeUnit = atoi(strSub.c_str());
		}
		break;
		}
	}
	int NumPrintTimes = OutputTimes.size();
	if (NumPrintTimes < Max_NumPrintTimes) {
		for (int i = NumPrintTimes; i < Max_NumPrintTimes; i++) {
			double rtem;
			if (i == 0) rtem = 0.0;
			else rtem = OutputTimes[i - 1];
			OutputTimes.push_back(rtem + PrintTimeIncrement);
		}
	}

	NumPrintTimes = OutputTimes.size();
	for (int i = 0; i < NumPrintTimes; i++) {
		if (timeUnit == 1) OutputTimes[i] = OutputTimes[i] * 3600.0;   //use hour as  input
		if (timeUnit == 2) OutputTimes[i] = OutputTimes[i] * 86400.0;  //use day as input
		if (timeUnit == 3) OutputTimes[i] = OutputTimes[i] * 86400.0 * 365.24;  //use year as input
	}
}

// Test of FLUID
void readFLUID() {
	//using namespace Fluid_Proterties;
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Flow;

	//if (myid == iMaster) {
	//	FLUID_read = true;
	LogFile << "--->>> Start reading FLUID" << std::endl << std::endl;
	//}
	struct Fluid {
		std::string FluidName;
		double Density;              // Fluid density(kg / m3)
		double Beta;                 // Fluid compressibility factor,  beta. rho_air= beta*p_air (kg/m**3/Pa)
		double Visco;                // Fluid viscosity (Pa.s)
		double ThermalExansionC;     // Thermal expansion coefficient(1/ C^o)
		double BinghamG;             // 2 parameter G   (Pa/m)
		double PowerLawN;            // Power-law n;
		double PowerLawH;            // Power-law H;
		double MinViscoPLFluid;      // Min.viscosity for power - law fluid, default = 1.e-5 (Pa.s)
		double ThermalCondc;         // Fluid thermal conductivity (W / m / C)
		short FluidType;
		double ref_P, ref_T;
	};
	std::vector<Fluid> Fluids;
	std::vector<int> AryIdx;
	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("FluidName");
	AryDescribles.push_back("FluidDensity");
	AryDescribles.push_back("FluidCom");
	AryDescribles.push_back("FluidVis");
	AryDescribles.push_back("ThermalExansionC");
	AryDescribles.push_back("FluidBingG");
	AryDescribles.push_back("FluidPowN");
	AryDescribles.push_back("FluidPowH");
	AryDescribles.push_back("FluidMinVis");
	AryDescribles.push_back("FluidTpye");
	AryDescribles.push_back("referT");
	AryDescribles.push_back("referP");
	AryDescribles.push_back("ThermalCondc");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}
	struct Fluid Fluid_Tmp;
	memset(&Fluid_Tmp, 0, sizeof(Fluid_Tmp));
	int nIdx = -1;
	//cout << sizeof(Fluid_Tmp);
	std::string strSub = "";
	Fluid_Tmp.FluidName = " ";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		if (std::find(AryIdx.begin(), AryIdx.end(), nIdx) != AryIdx.end())
		{
			if (nIdx != 0) {
				LogFile << "Warning: Parameter was redefined at text line(" << Line_Num << "): '" << LineText << "'. It was negleted!" << std::endl;
				continue;
			}
			AryIdx.clear();
			//if (Fluid_Tmp.FluidType < 0 || Fluid_Tmp.FluidType>2) Fluid_Tmp.FluidType = 0;
			//Fluids.push_back(Fluid_Tmp);
			memset(&Fluid_Tmp, 0, sizeof(Fluid_Tmp));
		}
		if (std::find(AryIdx.begin(), AryIdx.end(), nIdx) == AryIdx.end()) AryIdx.push_back(nIdx);

		switch (nIdx)
		{
		case 0:
		{
			Fluid_Tmp.FluidName = strSub.c_str();	// Read the input fulidname.
		}
		break;
		case 1:
		{
			if (strSub.length() > 0) Fluid_Tmp.Density = atof(strSub.c_str());
		}
		break;
		case 2:
		{
			if (strSub.length() > 0) Fluid_Tmp.Beta = atof(strSub.c_str());
		}
		break;
		case 3:
		{
			if (strSub.length() > 0) Fluid_Tmp.Visco = atof(strSub.c_str());
		}
		break;
		case 4:
		{
			if (strSub.length() > 0) Fluid_Tmp.ThermalExansionC = atof(strSub.c_str());
		}
		break;
		case 5:
		{
			if (strSub.length() > 0) Fluid_Tmp.BinghamG = atof(strSub.c_str());
		}
		break;
		case 6:
		{
			if (strSub.length() > 0) Fluid_Tmp.PowerLawN = atof(strSub.c_str());
		}
		break;
		case 7:
		{
			if (strSub.length() > 0) Fluid_Tmp.PowerLawH = atof(strSub.c_str());
		}
		break;
		case 8:
		{
			if (strSub.length() > 0) Fluid_Tmp.MinViscoPLFluid = atof(strSub.c_str());
		}
		break;
		case 9:
		{
			if (strSub.length() > 0) Fluid_Tmp.FluidType = atoi(strSub.c_str());
		}
		break;
		case 10:
		{
			if (strSub.length() > 0) Fluid_Tmp.ref_T = atof(strSub.c_str());
		}
		break;
		case 11:
		{
			if (strSub.length() > 0) Fluid_Tmp.ref_P = atof(strSub.c_str());
		}
		break;
		case 12:
		{
			if (strSub.length() > 0) Fluid_Tmp.ThermalCondc = atof(strSub.c_str());
		}
		break;
		}
	}
	{
		AryIdx.clear();
		if (Fluid_Tmp.FluidType < 0 || Fluid_Tmp.FluidType>2) Fluid_Tmp.FluidType = 0;
		Fluids.push_back(Fluid_Tmp);
		memset(&Fluid_Tmp, 0, sizeof(Fluid_Tmp));
	}
	//

	if (Fluids.size() == 1) {
		Fluids.push_back(Fluid_Tmp);
		Fluids.push_back(Fluid_Tmp);
		Fluids[1].Density = 0.0;
		Fluids[1].ThermalCondc = 0.5918;   //treated as water
		Fluids[1].FluidName = "Water";
		cout << Fluids[0].FluidName << endl;
		Fluids[1].FluidType = 0;
		Fluids[2].Beta = Fluids[1].Beta;
		Fluids[2].Density = Fluids[1].Density;
		Fluids[2].Visco = Fluids[1].Visco;
		Fluids[2].ThermalCondc = 0.5918;
		Fluids[2].FluidName = "Aqu2";
		Fluids[2].FluidType = 0;
	}
	else if (Fluids.size() == 2) {
		Fluids.push_back(Fluid_Tmp);
		Fluids[2].Beta = Fluids[1].Beta;
		Fluids[2].Density = Fluids[1].Density;
		Fluids[2].Visco = Fluids[1].Visco;
		Fluids[2].ThermalCondc = 0.5918;
		Fluids[2].FluidName = "Aqu2";
		Fluids[2].FluidType = 0;
	}
	if (Fluids[2].FluidName == "") Fluids[2].FluidName = "Aqu2";
	if (Fluids[0].FluidName == "") Fluids[0].FluidName = "Gas";
	if (Fluids[1].FluidName == "") Fluids[1].FluidName = "Water";

	if (Fluids[1].PowerLawN > 0.0 && Fluids[1].PowerLawH > 0.0) Fluids[1].FluidType = 1;
	if (Fluids[1].Visco > 0.0 && Fluids[1].BinghamG > 0.0) Fluids[1].FluidType = 2;
	if (Fluids[2].PowerLawN > 0.0 && Fluids[2].PowerLawH > 0.0) Fluids[2].FluidType = 1;
	if (Fluids[2].Visco > 0.0 && Fluids[2].BinghamG > 0.0) Fluids[2].FluidType = 2;

	if (Fluids[1].FluidType == 2 && Fluids[1].BinghamG == 0.0) {
		LogFile << "Warning:: The fluid '" << Fluids[1].FluidName << "' is specified as 2 fluid,but BinghanG=0.0, which must be >0.0 for 2 fluid.\n";
		LogFile << "    '" << Fluids[1].FluidName << "' will be simulated as 0 fluid" << std::endl;
	}
	if (Fluids[2].FluidType == 2 && Fluids[2].BinghamG == 0.0) {
		LogFile << "Warning:: The fluid '" << Fluids[2].FluidName << "' is specified as 2 fluid,but BinghanG=0.0, which must be >0.0 for 2 fluid.\n";
		LogFile << "    '" << Fluids[2].FluidName << "' will be simulated as 0 fluid" << std::endl;
	}

	if (Fluids[1].FluidType == 1) Fluids[1].Visco = 1.0;
	if (Fluids[2].FluidType == 1) Fluids[2].Visco = 1.0;

	if (Fluids[1].FluidType == 1)
	{
		if (Fluids[1].MinViscoPLFluid == 0.0) Fluids[1].MinViscoPLFluid = 1.0e-5;
	}
	if (Fluids[2].FluidType == 1) {
		if (Fluids[2].MinViscoPLFluid == 0.0) Fluids[2].MinViscoPLFluid = 1.0e-5;
	}

	if (Fluids[0].ThermalCondc == 0.0) Fluids[0].ThermalCondc = 0.0262;   //treated as air
	if (Fluids[1].ThermalCondc == 0.0) Fluids[1].ThermalCondc = 0.5918;   //treated as water
	if (Fluids[2].ThermalCondc == 0.0) Fluids[2].ThermalCondc = 0.5918;

	//CompuDefaultFluidProperties(ref_pressure, ref_temperature);
	//if (Fluids[1].FluidType == 1 || Fluids[2].FluidType == 1) IniNonNewtonianPLViscosity();

	//if (myid == iMaster) {
	LogFile << std::endl << "--->>> End reading FLUID" << std::endl << std::endl;
	LogFile << "===============================================================" << std::endl << std::endl;
	//}
}

// Test of PVTDA
void readPVTDA() {
	//using namespace PVT_Data;
	//using namespace Geologic_Media_Properties;
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//if (myid == iMaster) {
	//	PVTDA_read = true;
	LogFile << "--->>> Start reading PVTDA" << std::endl << std::endl;
	//}
	struct PVT_Table_Contents {
		double P_Oil;          // pressure
		double B_Oil;          // oil formation volume factor at P=P_Oil
		double B_Gas;          // gas formation volume factor at P=P_Oil
		double Rs_GO;          // solution gas-oil ratio at P=p_Oil
		double Rs_GW;          // solution gas-water ratio at P=P_water
		double Vis_O;          // oil viscosity at P=P_Oil
		double Vis_G;          // gas viscosity at P=P_Gas

	};

	struct PVT {
		std::vector <PVT_Table_Contents> tb;
		double Pb_Ori;         // original bubble point pressure (Pa).
		double B_Water_Ori;    // water formation volume factor at Po=Pb_ori (default=1)
		short int Num_PVT;    // number of entries of PVT Table
		short int Num_PVT3;   // number of entries of PVT Table with P_Oil <= Pb_ori
		double B_Oil_Slp;     // slope of P_Oil vs. B_Oil for P_Oil > Pb_Ori
		double Rsw_Slp;
		double Rs_Slp;        // slope of P_Oil vs. Rs_GO for P_Oil > Pb_Ori
		double Viso_Slp;
		double temperature;
	};
	std::vector<PVT> PVTs;
	int Num_Temperature;
	int NumPhases = 3;
	struct PVT_Table_Contents tb_tmp;
	struct PVT PVTs_tmp;
	int nt_pvt = 0;
	std::vector<std::string> vctAry;

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("PVTatT");
	AryDescribles.push_back("PB_I");
	AryDescribles.push_back("BW_I");
	AryDescribles.push_back("PVT_TABLE");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	memset(&PVTs_tmp, 0, sizeof(PVTs_tmp));
	bool bChange = false;
	int nIdx = -1;
	std::string strSub = "";

	while (true)
	{
		if (nIdx == 3) nIdx = -nIdx - 10;
		else nIdx = -1;

		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		switch (nIdx)
		{
		case 0:
		{
			if (bChange) PVTs.push_back(PVTs_tmp);
			PVTs_tmp.tb.clear();
			memset(&PVTs_tmp, 0, sizeof(PVTs_tmp));
			bChange = false;
			PVTs_tmp.temperature = atof(strSub.c_str());
		}
		break;
		case 1:
		{
			PVTs_tmp.Pb_Ori = atof(strSub.c_str());
			bChange = true;
		}
		break;
		case 2:
		{
			PVTs_tmp.B_Water_Ori = atof(strSub.c_str());
			bChange = true;
		}
		break;
		case 3:
		{
			vctAry.clear();
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				vctAry.push_back(p);
				p = strtok(NULL, ",");
			}
			if (vctAry.size() > 0) tb_tmp.P_Oil = atof(vctAry[0].c_str());
			if (vctAry.size() > 1) tb_tmp.B_Oil = atof(vctAry[1].c_str());
			if (vctAry.size() > 2) tb_tmp.B_Gas = atof(vctAry[2].c_str());
			if (vctAry.size() > 3) tb_tmp.Rs_GO = atof(vctAry[3].c_str());
			if (vctAry.size() > 4) tb_tmp.Vis_O = atof(vctAry[4].c_str());
			if (vctAry.size() > 5) tb_tmp.Vis_G = atof(vctAry[5].c_str());
			if (vctAry.size() > 6) tb_tmp.Rs_GW = atof(vctAry[6].c_str());

			PVTs_tmp.tb.push_back(tb_tmp);
			bChange = true;
		}
		break;
		}
	}
	if (bChange) PVTs.push_back(PVTs_tmp);
	nt_pvt = PVTs.size();
	//
	for (int j = 0; j < nt_pvt; j++) {
		PVTs_tmp = PVTs[j];
		PVTs_tmp.Num_PVT = PVTs_tmp.tb.size();
		if (PVTs_tmp.Num_PVT > 0) {
			PVTs_tmp.Num_PVT3 = 1;
			for (int k = 0; k < PVTs_tmp.tb.size(); k++) {
				if (k > 0) {
					if (PVTs_tmp.tb[k - 1].P_Oil <= PVTs_tmp.Pb_Ori &&
						PVTs_tmp.tb[k].P_Oil > PVTs_tmp.Pb_Ori) PVTs_tmp.Num_PVT3 = k;
				}
			}
			//
			if (NumPhases == 3) {
				int iTem = PVTs_tmp.Num_PVT - 1;
				if (PVTs_tmp.tb[0].P_Oil >= PVTs_tmp.Pb_Ori || PVTs_tmp.tb[iTem].P_Oil <= PVTs_tmp.Pb_Ori) {
					LogFile << "Warening: for three phase flow, the input pressures in PVT table should cover the range for both larger and less\n"
						<< "   than its initial bubble point pressure. Please check the PVT table at temperature: "
						<< PVTs_tmp.temperature << " oC." << std::endl;
				}
			}

			int n_slop = PVTs_tmp.Num_PVT - PVTs_tmp.Num_PVT3;
			double boslp_s = 0.0;
			double rsslp_s = 0.0;
			double visoslp_s = 0.0;
			for (int k = 0; k < n_slop; k++) {
				int k_npvt3 = PVTs_tmp.Num_PVT3 + k;
				double boslp_i = (PVTs_tmp.tb[k_npvt3].B_Oil - PVTs_tmp.tb[k_npvt3 - 1].B_Oil) /
					(PVTs_tmp.tb[k_npvt3].P_Oil - PVTs_tmp.tb[k_npvt3 - 1].P_Oil);
				double rsslp_i = (PVTs_tmp.tb[k_npvt3].Rs_GO - PVTs_tmp.tb[k_npvt3 - 1].Rs_GO) /
					(PVTs_tmp.tb[k_npvt3].P_Oil - PVTs_tmp.tb[k_npvt3 - 1].P_Oil);
				double visoslp_i = (PVTs_tmp.tb[k_npvt3].Vis_O - PVTs_tmp.tb[k_npvt3 - 1].Vis_O) /
					(PVTs_tmp.tb[k_npvt3].P_Oil - PVTs_tmp.tb[k_npvt3 - 1].P_Oil);
				boslp_s = boslp_s + boslp_i;
				rsslp_s = rsslp_s + rsslp_i;
				visoslp_s = visoslp_s + visoslp_i;
			}
			if (n_slop >= 1) {
				PVTs_tmp.B_Oil_Slp = boslp_s / n_slop;
				PVTs_tmp.Rs_Slp = rsslp_s / n_slop;
				PVTs_tmp.Viso_Slp = visoslp_s / n_slop;
			}
		}
		if (PVTs_tmp.B_Water_Ori < 0) { PVTs_tmp.B_Water_Ori = 1.00; }
		PVTs[j] = PVTs_tmp;
	}
	Num_Temperature = nt_pvt;

	for (int i = 0; i < nt_pvt; i++)
	{
		for (int j = i + 1; j < nt_pvt; j++)
		{
			if (PVTs[i].temperature <= PVTs[j].temperature) continue;
			PVTs_tmp = PVTs[i];
			PVTs[i] = PVTs[j];
			PVTs[j] = PVTs_tmp;
		}
	}

	//if (myid == iMaster) {
	LogFile << "Have read " << nt_pvt << " set PVT data " << std::endl;
	LogFile << std::endl << "--->>> End reading PVTDA" << std::endl << std::endl;
	LogFile << "===============================================================" << std::endl << std::endl;
	//}
}

// Test of FOFT, COFT, GOFT

void readFOFT() {
	//using namespace Time_Series_Parameters;
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	//using namespace General_Control_Parameters;
	struct Observation_Elem {
		int num;
		std::string name;
	};
	struct Observation_Conx {
		int num;
		std::string Ele1Name, Ele2Name;
	};
	struct Observation_SS {
		int num;
		std::string name;
	};
	std::vector<Observation_Elem> ObsElem;   //The observation element list
	std::vector<Observation_Conx> ObsConx;   //The observation connection list
	std::vector<Observation_SS> ObsSS;       //The observation SS list
	int NumObsElem;                          //The # of elements in the observation element list
	int NumObsConx;                          //The # connections in the observation connection list
	int NumObsSS;                            //The # of sources / sinks in the observation SS list
	short ObsElemVar[11], ObsConxVar[11], ObsSSVar[11];
	bool FirstVisitTS = true;
	//if (myid == iMaster) {
	//	ELEME_read = true;
	LogFile << std::endl << "--->>> Start reading FOFTS" << std::endl;
	//}

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("OBElements");
	AryDescribles.push_back("ObsElemVar");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}
	struct Observation_Elem ObsElem_Tmp;

	for (int j = 0; j <= 10; j++) ObsElemVar[j] = 0;

	NumObsElem = 0;
	int i_Num = 1;
	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		switch (nIdx)
		{
		case 0:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				ObsElem_Tmp.num = i_Num; i_Num++;
				strTmp = p;
				strTmp.erase(0, strTmp.find_first_not_of(" "));
				strTmp.erase(strTmp.find_last_not_of(" ") + 1);
				ObsElem_Tmp.name = strTmp;
				ObsElem.push_back(ObsElem_Tmp);
				p = strtok(NULL, ",");
			}
		}
		break;
		case 1:
		{
			int j = 0;
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				ObsElemVar[j] = atoi(p);
				j++;

				if (j >= 10) break;
				p = strtok(NULL, ",");
			}
		}
		break;
		}
	}
	NumObsElem = i_Num - 1;
	if (NumObsElem > 0) {
		FirstVisitTS = true;
	}

	//if (myid == iMaster) {
	LogFile << std::endl << "Have read " << NumObsElem << " grid(s) for time series output." << std::endl;
	LogFile << std::endl << "--->>> End reading FOFTS" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}
}

void readCOFT() {
	//using namespace Time_Series_Parameters;
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	struct Observation_Elem {
		int num;
		std::string name;
	};
	struct Observation_Conx {
		int num;
		std::string Ele1Name, Ele2Name;
	};
	struct Observation_SS {
		int num;
		std::string name;
	};
	std::vector<Observation_Elem> ObsElem;   //The observation element list
	std::vector<Observation_Conx> ObsConx;   //The observation connection list
	std::vector<Observation_SS> ObsSS;       //The observation SS list
	int NumObsElem;                          //The # of elements in the observation element list
	int NumObsConx;                          //The # connections in the observation connection list
	int NumObsSS;                            //The # of sources / sinks in the observation SS list
	short ObsElemVar[11], ObsConxVar[11], ObsSSVar[11];
	bool FirstVisitTS = true;
	int ElemNameLength = 5;
	//if (myid == iMaster) {
	//	Conx__read = true;
	LogFile << std::endl << "--->>> Start reading COFTS" << std::endl;
	//}

	struct Observation_Conx ObsConx_Tmp;
	NumObsConx = 0;
	int i_Num = 1;

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("OBConnections");
	AryDescribles.push_back("ObsConxVar");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	for (int j = 0; j <= 10; j++) ObsConxVar[j] = 0;
	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		switch (nIdx)
		{
		case 0:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				ObsConx_Tmp.num = i_Num; i_Num++;
				strTmp = p;
				strTmp.erase(0, strTmp.find_first_not_of(" "));
				strTmp.erase(strTmp.find_last_not_of(" ") + 1);
				ObsConx_Tmp.Ele1Name = strTmp.substr(0, ElemNameLength);
				ObsConx_Tmp.Ele2Name = strTmp.substr(ElemNameLength, ElemNameLength);
				ObsConx.push_back(ObsConx_Tmp);
				p = strtok(NULL, ",");
			}
		}
		break;
		case 1:
		{
			int j = 0;
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				ObsConxVar[j] = atoi(p);
				j++;

				if (j >= 10) break;
				p = strtok(NULL, ",");
			}
		}
		break;
		}
	}
	NumObsConx = i_Num - 1;
	//if (myid == iMaster) {
	LogFile << std::endl << "Have read " << NumObsConx << " connection(s) for time series output." << std::endl;
	LogFile << std::endl << "--->>> End reading COFTS" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}
}

void readGOFT() {
	//using namespace Time_Series_Parameters;
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	struct Observation_Elem {
		int num;
		std::string name;
	};
	struct Observation_Conx {
		int num;
		std::string Ele1Name, Ele2Name;
	};
	struct Observation_SS {
		int num;
		std::string name;
	};
	std::vector<Observation_Elem> ObsElem;   //The observation element list
	std::vector<Observation_Conx> ObsConx;   //The observation connection list
	std::vector<Observation_SS> ObsSS;       //The observation SS list
	int NumObsElem;                          //The # of elements in the observation element list
	int NumObsConx;                          //The # connections in the observation connection list
	int NumObsSS;                            //The # of sources / sinks in the observation SS list
	short ObsElemVar[11], ObsConxVar[11], ObsSSVar[11];
	bool FirstVisitTS = true;
	//if (myid == iMaster) {
	//	SS_Ti_read = true;
	LogFile << std::endl << "--->>> Start reading GOFTS" << std::endl;
	//}

	std::string strTmp;
	int i_Num = 1;
	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("OBSS");
	struct Observation_SS NumObsSS_Tmp;
	NumObsSS = 0;
	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		switch (nIdx)
		{
		case 0:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				NumObsSS_Tmp.num = i_Num; i_Num++;
				strTmp = p;
				strTmp.erase(0, strTmp.find_first_not_of(" "));
				strTmp.erase(strTmp.find_last_not_of(" ") + 1);
				NumObsSS_Tmp.name = strTmp;
				ObsSS.push_back(NumObsSS_Tmp);
				p = strtok(NULL, ",");
			}
		}
		break;
		}
	}
	NumObsSS = i_Num - 1;
	//for (int i = 0; i < NumObsSS; i++) std::cout << ObsSS[i].name << std::endl;
	//if (myid == iMaster) {
	LogFile << std::endl << "Have read " << NumObsSS << " grid(s)" << std::endl;
	LogFile << std::endl << "--->>> End reading GOFTS" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}
}

std::string trim(std::string& s)
{
	if (s.empty()) { return s; }
	s.erase(0, s.find_first_not_of(" "));
	s.erase(s.find_last_not_of(" ") + 1);
	return s;
}

// Test of MESH_ELEME reading
void readMESH_Eleme() //Seperate the keywords eleme and conne
{
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	//if (myid != iMaster) return;
	std::ofstream saveMesh;
	saveMesh.open("MESH");
	//ELEME_read = true;
	LogFile << "--->>> Start reading ELEME" << std::endl << std::endl;
	saveMesh << "ELEME" << std::endl;
	getline(MainInputFile, LineText);
	Line_Num++;
	std::vector <std::string> fields;
	char szOutData[255]; double fTmp;
	while (trim(LineText).length() > 0) {
		LineText = trim(LineText);
		if (LineText.length() == 0 || LineText.compare("ENDFI") == 0 || LineText.compare("/") == 0 || LineText.compare("END_SECTION") == 0) break;

		fields.clear();
		boost::split(fields, LineText, boost::is_any_of(","));

		if (fields.size() <= 6)
		{
			while (fields.size() < 6) fields.push_back("");
			sprintf(szOutData, "%-15s%-5s", fields[0].c_str(), fields[1].c_str());
			saveMesh << szOutData;
			for (int i = 2; i < 6; i++)
			{
				fTmp = atof(fields[i].c_str());
				if (fTmp < 0)
				{
					sprintf(szOutData, "%10.3e", fTmp);
				}
				else
				{
					sprintf(szOutData, "%10.4e", fTmp);
				}
				saveMesh << szOutData;
			}
			saveMesh << std::endl;
		}
		else
		{
			while (fields.size() < 9) fields.push_back("");
			sprintf(szOutData, "%-15s%-5s", fields[0].c_str(), fields[1].c_str());
			saveMesh << szOutData;
			for (int i = 2; i < 9; i++)
			{
				fTmp = atof(fields[i].c_str());
				if (fTmp < 0)
				{
					sprintf(szOutData, "%10.3e", fTmp);
				}
				else
				{
					sprintf(szOutData, "%10.4e", fTmp);
				}
				saveMesh << szOutData;
			}
			saveMesh << std::endl;
		}
		getline(MainInputFile, LineText);
		Line_Num++;
	}
	saveMesh << std::endl;
	saveMesh.close();

	LogFile << std::endl << "--->>> End reading ELEME" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl << std::endl;
}

// Test of CONNE reading
void readMESH_Conne() {
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;

	//if (myid != iMaster) return;
	int ElemNameLength = 5;
	std::ofstream saveMesh;
	saveMesh.open("MESH", std::ios::app);
	//CONNE_read = true;
	LogFile << "--->>> Start reading CONNE" << std::endl << std::endl;
	saveMesh << "CONNE" << std::endl;
	getline(MainInputFile, LineText);
	Line_Num++;
	std::vector <std::string> fields;
	char szOutData[255]; double fTmp;
	while (trim(LineText).length() > 0) {
		LineText = trim(LineText);
		if (LineText.length() == 0 || LineText.compare("/") == 0 || LineText.compare("ENDFI") == 0 || LineText.compare("END_SECTION") == 0) break;

		fields.clear();
		boost::split(fields, LineText, boost::is_any_of(","));

		while (fields.size() < 8) fields.push_back("");
		saveMesh << std::setw(ElemNameLength) << fields[0];
		saveMesh << std::setw(ElemNameLength) << fields[1];
		saveMesh << std::setw(28 - 2 * ElemNameLength) << "";
		saveMesh << std::setw(2) << atoi(fields[2].c_str());
		//cout << fields[2] << endl;
		//cout << fields[2].c_str() << endl;
		for (int i = 3; i < 8; i++)
		{
			fTmp = atof(fields[i].c_str());
			if (fTmp < 0)
			{
				sprintf(szOutData, "%10.3e", fTmp);
			}
			else
			{
				sprintf(szOutData, "%10.4e", fTmp);
			}
			saveMesh << szOutData;
		}
		saveMesh << std::endl;
		getline(MainInputFile, LineText);
		Line_Num++;
	}
	saveMesh << std::endl;
	saveMesh.close();

	LogFile << std::endl << "--->>> End reading CONNE" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl << std::endl;
}

// Test of INCON reading
void readINCON()
{
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	std::string str_tmp = " ";
	short nPrimary = 6;
	//if (myid != iMaster) return;
	//INCON_read = true;
	LogFile << "--->>> Start reading INCON" << std::endl;
	std::ofstream incon_out;
	incon_out.open("INCON");	// INCON input
	incon_out << "INCON" << std::endl;

	getline(MainInputFile, LineText);
	Line_Num += 1;
	std::string mainInputFileName;
	std::string strDescrible = "ELEMENTNAME";
	std::string strName = "";
	std::string strOut = "";
	double xn[10];
	while (LineText.length() > 0)
	{
		if (!LineText.empty())
		{
			LineText.erase(0, LineText.find_first_not_of(" "));
			LineText.erase(LineText.find_last_not_of(" ") + 1);
		}

		if (LineText.length() == 0 || LineText.compare("/") == 0) break;
		if (LineText.length() < 2 || LineText.substr(0, 2).compare("//") == 0)
		{
			getline(MainInputFile, LineText);
			Line_Num += 1;
			continue;
		}

		std::string strTmp = LineText;
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);

		int nPos = strTmp.find(strDescrible);
		if (nPos != -1)
		{
			LineText = LineText.substr(nPos + strDescrible.size());
			nPos = LineText.rfind(':');
			strTmp = LineText.substr(nPos + 1);
		}
		else
		{
			strTmp = LineText;
		}
		int nSub;
		if (!strTmp.empty())
		{
			strTmp.erase(0, strTmp.find_first_not_of(" "));
			nSub = strTmp.find("//");
			if (nSub != -1) strTmp.erase(nSub);
			strTmp.erase(strTmp.find_last_not_of(" ") + 1);
		}

		getline(MainInputFile, LineText);
		Line_Num += 1;
		nSub = LineText.find("//");
		if (nSub != -1) LineText.erase(nSub);
		std::size_t found = LineText.find_first_not_of("0123456789+-Ee., ");
		if (found != std::string::npos) {
			LogFile << "Error in reading file '" << mainInputFileName << "' at the text line: " << Line_Num;
			LogFile << " for reading expected double variables at: '" << LineText << "'" << std::endl;
			std::cout << "Abnormal stop, please see the log file for details." << std::endl;
			exit(1);
		}

		strName = strTmp;
		char* p = strtok(const_cast<char*>(LineText.c_str()), ",");
		int nTmp = 0;

		for (int i = 0; i < 10; i++) xn[i] = -999.666;
		while (p != NULL)
		{
			xn[nTmp++] = atof(p);
			p = strtok(NULL, ",");
		}
		nPrimary = nTmp;
		strOut = strName;
		if (strName.length() == 0) { getline(MainInputFile, LineText); Line_Num += 1; continue; }
		for (int i = strName.length(); i < 15; i++) strOut += " ";
		incon_out << strOut;

		for (int i = 0; i < nPrimary - 1; i++)
		{
			if (xn[i] != -999.666) incon_out << std::right << std::setw(20) << std::setprecision(10) << std::scientific << xn[i];
		}
		if (xn[nPrimary - 1] != -999) incon_out << std::right << std::setw(20) << std::dec << (int)xn[nPrimary - 1];
		incon_out << std::endl;

		getline(MainInputFile, LineText);
		Line_Num += 1;
	}
	incon_out.close();

	LogFile << std::endl << "--->>> End reading INCON" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl << std::endl;
}

//Test of MODEL
void readMODEL()
{
	//using namespace Fluid_Proterties;
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Flow;

	//if (myid == iMaster) {
	//	MODEL_read = true;
	LogFile << "--->>> Start reading MODEL" << std::endl << std::endl;
	//}
	bool noAqu2, noGas, AquMixable, iRelaxedV, AccountingForDiffusion, nonIsothermal, realGasAccountingVapor;
	int ElemNameLength, NumGases, rockPermeabilitUnit;
	std::string eos_type = "PR";
	std::string ComponentName[10];
	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("HasAqu2");
	AryDescribles.push_back("HasGas");
	AryDescribles.push_back("AquMixable");
	AryDescribles.push_back("ElemNameLength");
	AryDescribles.push_back("RelaxedV");
	AryDescribles.push_back("AFDiffusion");
	AryDescribles.push_back("PermUnit");
	AryDescribles.push_back("nIsoTherm");
	AryDescribles.push_back("EosType");
	AryDescribles.push_back("RealGAcV");
	AryDescribles.push_back("GasCompo");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}


	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		switch (nIdx)
		{
		case 0:
		{
			if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { noAqu2 = false; }
			if (trim(strSub) == "FALSE" || trim(strSub) == "false" || trim(strSub) == "False") { noAqu2 = true; }
		}
		break;
		case 1:
		{
			if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { noGas = false; NumGases = 1; }
			if (trim(strSub) == "FALSE" || trim(strSub) == "false" || trim(strSub) == "False") { noGas = true; NumGases = 0; }
		}
		break;
		case 2:
		{
			if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { AquMixable = true; }
			if (trim(strSub) == "FALSE" || trim(strSub) == "false" || trim(strSub) == "False") { AquMixable = false; }
		}
		break;
		case 3:
		{
			if (strSub.length() > 0) ElemNameLength = atoi(strSub.c_str());
			if (ElemNameLength < 3 || ElemNameLength>13) ElemNameLength = 5;
		}
		break;
		case 4:
		{
			iRelaxedV = false;
			if (strSub.length() > 0)
			{
				if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { iRelaxedV = true; }
			}
		}
		break;
		case 5:
		{
			AccountingForDiffusion = false;
			if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { AccountingForDiffusion = true; };
		}
		break;
		case 6:
		{
			if (strSub.length() > 0) rockPermeabilitUnit = atoi(strSub.c_str());
			if (rockPermeabilitUnit < 0 || rockPermeabilitUnit>2) rockPermeabilitUnit = 0;
		}
		break;
		case 7:
		{
			if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { nonIsothermal = true; }
			if (trim(strSub) == "FALSE" || trim(strSub) == "false" || trim(strSub) == "False") { nonIsothermal = false; }
		}
		break;
		case 8:
		{
			std::string sTem = trim(strSub);
			boost::to_upper(sTem);
			if (sTem.length() > 2) eos_type = sTem;
		}
		break;
		case 9:
		{
			if (trim(strSub) == "FALSE" || trim(strSub) == "false" || trim(strSub) == "False") { realGasAccountingVapor = false; }
		}
		break;
		case 10:
		{
			strTmp = strSub;
			char* p = strtok(const_cast<char*>(strTmp.c_str()), ",");
			int nTmp = 0;
			std::string s1;
			while (p != NULL)
			{
				switch (nTmp)
				{
				case 0:
					s1 = p;
					ComponentName[0] = trim(s1);
					break;
				case 1:
					s1 = p;
					ComponentName[3] = trim(s1);
					break;
				case 2:
					s1 = p;
					ComponentName[4] = trim(s1);
					break;
				case 3:
					break;
					s1 = p;
					ComponentName[5] = trim(s1);
					break;
				}
				nTmp++;
				NumGases = nTmp;
				if (nTmp == 4) break;
				p = strtok(NULL, ",");
			}
			//
			if (NumGases < 0 || NumGases>4) {
				//if (myid == iMaster) {
				LogFile << "Warning: Gas component number (NumGases) must be between  0-4. the input NumGases=" << NumGases << std::endl;
				LogFile << "Simulation stop!" << std::endl;
				//}
				exit(0);
			}
			if (NumGases == 0) noGas = true;
			else noGas = false;
		}
		break;
		}
	}

	if (noAqu2) {
		ComponentName[2] = ComponentName[3];
		ComponentName[3] = ComponentName[4];
		ComponentName[4] = ComponentName[5];
	}

	//if (myid == iMaster) {
	LogFile << std::endl << "--->>> End reading MODEL" << std::endl << std::endl;
	LogFile << "===============================================================" << std::endl << std::endl;
	//}
}

// Test of RPCAP
void readRPCAP() {
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Geologic_Media_Properties;

	//if (myid == iMaster) {
	//	RPCAP_read = true;
	LogFile << std::endl << "--->>> Start reading RPCAP" << std::endl;
	//}
	struct PFMedium {
		std::string MediumName;
		int RelPermEquationNum;      // The number of the relative permeability function
		int PcapEquationNum;         // The number of the capillary function
//		int PhiPolyOrder;            // The order of the porosity vs.DP polynomial
//
		double DensG;                // The PF medium grain density(kg / m3)
		double Poros;                // The PF medium porosity
		double SpcHt;                // The PF medium specific heat(J / kg / C)
		double KThrW;                // The thermal conductivity of the saturated PF medium(W / m / C)
		double KThrD;                // The thermal conductivity of the dry PF medium(W / m / C)
		double Cbeta;                // high velocity non-Darcy parameter C_beta; =0 for Darcy flow.
									 //
		double Compr;                // The pore compressibility(1 / Pa)
		double Expan;                // The pore expansivity(1 / C)
		double Tortu;                //  The PF medium tortuosity factor
		double Klink;                // The Kinkenberg parameter(1 / Pa)

//		double Sw_irr;               // The irreducible water saturation
//		double Sg_irr;               // The irreducible gas saturation
//		double So_irr;               // The irreducible oil phase saturation
//		double Pc_max;               // The maximum capillary pressure

		double CritSat;              // The critical mobile saturation
		double permExpon;             //The exponent for permeability reduction in the presence of solid phases
		double beta;                 // Empirical parameter for geomechanics - induced permeability reduction
		double gama;                 // Empirical parameter for geomechanics - induced permeability reduction

//		double HiComp;               // Upper compressibility limit
//		double LoComp;               // Lower compressibility limit
//		double SatAtHiComp;          // Saturation corresponding to the high compressibility
//		double SatAtLoComp;          // Saturation corresponding to the low compressibility
//		double DeltaSat;             // Empirical parameter for adjusting the smoothness of the scanning curve
//		bool CompAsSolidSatFunction;	
//		double PhiCoeff[7];
		double perm[3];              // Intrinsic permeabilities along principal axes(m ^ 2)
		double RelPermParam[14];     // Parameters of the relative permeability functions
		double PcapParam[14];        // Parameters of the capillary pressure functions

		int nTableOW_K_PC;
		int nTableOG_K_PC;;
		std::vector<double> sw_WO;
		std::vector<double> kw_WO;
		std::vector<double> ko_WO;
		std::vector<double> pc_WO;

		std::vector<double> sg_OG;
		std::vector<double> kg_OG;
		std::vector<double> ko_OG;
		std::vector<double> pc_OG;
		std::vector<double> pc_WG;
		double swInitial_Vug;
		double swEnd_Vug;
	};
	class Rocks {
	public:
		std::vector<PFMedium> media;
		//		std::vector<PFMedium1> PoMed;
		double RPD[7], CPD[7];
		int DefaultRelPermType, DefaultCapPresType;
		int MaxWaterTab, MaxGasTab;
		//
		//Rocks();
		//void Set_Irreducible_Saturations();
	};
	Rocks ModelRocks;

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("DRPT");
	AryDescribles.push_back("DCPT");
	AryDescribles.push_back("DRPParameters");
	AryDescribles.push_back("DCPParameters");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		switch (nIdx)
		{
		case 0:
		{
			ModelRocks.DefaultRelPermType = atoi(strSub.c_str());
		}
		break;
		case 1:
		{
			ModelRocks.DefaultCapPresType = atoi(strSub.c_str());
		}
		break;
		case 2:
		{
			int i = 0;
			for (i = 0; i < 7; i++) ModelRocks.RPD[i] = 0.0;
			i = 0;
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				if (i >= 7) break;
				ModelRocks.RPD[i] = atof(p);
				i++;
				p = strtok(NULL, ",");
			}
		}
		break;
		case 3:
		{
			int i = 0;
			for (i = 0; i < 7; i++) ModelRocks.CPD[i] = 0.0;
			i = 0;
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				if (i >= 7) break;
				ModelRocks.CPD[i] = atof(p);
				i++;
				p = strtok(NULL, ",");
			}
		}
		break;
		}
	}

	//if (myid == iMaster) {
	LogFile << std::endl << "--->>> End reading RPCAP" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}
}

// Test of TIMBC
void readTIMBC() {                                    // "NeedToBeChecked"
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Geologic_Media_Properties;
	//using namespace Transient_Boundary_Conditions;

	//if (myid == iMaster) {
	//	TIMBC_read = true;
	LogFile << std::endl << "--->>> Start reading TIMBC" << std::endl;
	//}
	struct Transient_Boundary_Type {
		std::string ElemName;
		int ElemNum;
		int numTimePoint;
		short int primaryVariableNo;
		std::vector<double> timebc;
		std::vector<double> primaryV;
	};
	std::vector<Transient_Boundary_Type> tBC;
	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("Ele_Name");
	AryDescribles.push_back("N_PrimaryV");
	AryDescribles.push_back("TimeUnit");
	AryDescribles.push_back("timeBC");
	AryDescribles.push_back("primary");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}
	int nIdx = -1;
	double timeFactor = 1.0;;
	bool bfirst = true;
	struct Transient_Boundary_Type tBC_Tmp;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;
		if (nIdx == 0)
		{
			if (bfirst)
			{
				bfirst = false;
			}
			else
			{
				tBC.push_back(tBC_Tmp);
				tBC_Tmp.timebc.clear();
				tBC_Tmp.primaryV.clear();
			}
		}

		switch (nIdx)
		{
		case 0:
		{
			tBC_Tmp.ElemName = strSub;
		}
		break;
		case 1:
		{
			tBC_Tmp.primaryVariableNo = atoi(strSub.c_str());
		}
		break;
		case 2:
		{
			int j = atoi(strSub.c_str());
			timeFactor = 1.0;
			if (j == 1) timeFactor = 60.0;   // minutes
			if (j == 2) timeFactor = 3600.0;  // hours 
			if (j == 3) timeFactor = 3600.0 * 24;   // days
			if (j == 4) timeFactor = 3600.0 * 24 * 365.24;  // years 
		}
		break;
		case 3:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				tBC_Tmp.timebc.push_back(atof(p) * timeFactor);
				p = strtok(NULL, ",");
			}
		}
		break;
		case 4:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			while (p != NULL)
			{
				tBC_Tmp.primaryV.push_back(atof(p));
				p = strtok(NULL, ",");
			}
		}
		break;
		}
	}

	if (!bfirst)
	{
		tBC.push_back(tBC_Tmp);
		tBC_Tmp.timebc.clear();
		tBC_Tmp.primaryV.clear();
	}
	int numTransientBound = tBC.size();
	//if (myid == iMaster) {
	LogFile << " Time dependent boundary points read, numTransientBound=" << numTransientBound << std::endl;
	LogFile << std::endl << "--->>> End reading TIMBC" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}
}

// Test of INDOM
void readINDOM() {
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Geologic_Media_Properties;
	//using namespace DomainInitialConditions;
	struct Domain_Initial_Condition {
		int RockN;
		double domain_initial_cond[10];
	};
	Domain_Initial_Condition Domain_Initial_Conditions_tmp;
	std::vector<Domain_Initial_Condition> Domain_Initial_Conditions;
	struct PFMedium {
		std::string MediumName;
		int RelPermEquationNum;      // The number of the relative permeability function
		int PcapEquationNum;         // The number of the capillary function
//		int PhiPolyOrder;            // The order of the porosity vs.DP polynomial
//
		double DensG;                // The PF medium grain density(kg / m3)
		double Poros;                // The PF medium porosity
		double SpcHt;                // The PF medium specific heat(J / kg / C)
		double KThrW;                // The thermal conductivity of the saturated PF medium(W / m / C)
		double KThrD;                // The thermal conductivity of the dry PF medium(W / m / C)
		double Cbeta;                // high velocity non-Darcy parameter C_beta; =0 for Darcy flow.
									 //
		double Compr;                // The pore compressibility(1 / Pa)
		double Expan;                // The pore expansivity(1 / C)
		double Tortu;                //  The PF medium tortuosity factor
		double Klink;                // The Kinkenberg parameter(1 / Pa)

//		double Sw_irr;               // The irreducible water saturation
//		double Sg_irr;               // The irreducible gas saturation
//		double So_irr;               // The irreducible oil phase saturation
//		double Pc_max;               // The maximum capillary pressure

		double CritSat;              // The critical mobile saturation
		double permExpon;             //The exponent for permeability reduction in the presence of solid phases
		double beta;                 // Empirical parameter for geomechanics - induced permeability reduction
		double gama;                 // Empirical parameter for geomechanics - induced permeability reduction

//		double HiComp;               // Upper compressibility limit
//		double LoComp;               // Lower compressibility limit
//		double SatAtHiComp;          // Saturation corresponding to the high compressibility
//		double SatAtLoComp;          // Saturation corresponding to the low compressibility
//		double DeltaSat;             // Empirical parameter for adjusting the smoothness of the scanning curve
//		bool CompAsSolidSatFunction;	
//		double PhiCoeff[7];
		double perm[3];              // Intrinsic permeabilities along principal axes(m ^ 2)
		double RelPermParam[14];     // Parameters of the relative permeability functions
		double PcapParam[14];        // Parameters of the capillary pressure functions

		int nTableOW_K_PC;
		int nTableOG_K_PC;;
		std::vector<double> sw_WO;
		std::vector<double> kw_WO;
		std::vector<double> ko_WO;
		std::vector<double> pc_WO;

		std::vector<double> sg_OG;
		std::vector<double> kg_OG;
		std::vector<double> ko_OG;
		std::vector<double> pc_OG;
		std::vector<double> pc_WG;
		double swInitial_Vug;
		double swEnd_Vug;
	};
	class Rocks {
	public:
		std::vector<PFMedium> media;
		//		std::vector<PFMedium1> PoMed;
		double RPD[7], CPD[7];
		int DefaultRelPermType, DefaultCapPresType;
		int MaxWaterTab, MaxGasTab;
		//
		//Rocks();
		//void Set_Irreducible_Saturations();
	};
	Rocks ModelRocks;
	//if (myid == iMaster) {
	//	INDOM_read = true;
	LogFile << std::endl << "--->>> Start reading INDOM" << std::endl;
	//	if (ModelRocks.media.size() <= 0) {	//check if keyword ROCKS have bben read or not
	//		LogFile << std::endl << "!!!!! Error: Keyword INDOM should be read after ROCKS  " << std::endl;
	//		pc->StopAll();
	//	}
	//}

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("Rock_Name");
	AryDescribles.push_back("IniCond");

	std::string MediumName = "co211";
	bool contains_MediaName = false;
	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}
	int nIdx = -1;
	int rN = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		switch (nIdx)
		{
		case 0:
		{
			if (strSub.length() < 5)
			{
				contains_MediaName = false;
				cout << strSub.length();
				break;
			}
			strSub = strSub.substr(0, 5);
			for (size_t i = 0; i < 2; i++) {
				if (strSub.compare(MediumName) != 0) continue;
				contains_MediaName = true;
				rN = i;
				break;
			}
		}
		break;
		case 1:
		{
			if (contains_MediaName)
			{
				Domain_Initial_Conditions_tmp.RockN = rN;
				char* p = strtok(const_cast<char*>(strSub.c_str()), ",");

				int i = 0;
				while (p != NULL)
				{
					if (i > 9) break;
					Domain_Initial_Conditions_tmp.domain_initial_cond[i] = atof(p);
					p = strtok(NULL, ",");
					i++;
				}
				Domain_Initial_Conditions.push_back(Domain_Initial_Conditions_tmp);
			}
			else
				/*{
					if (myid == iMaster) {
						LogFile << std::endl << "Warning: read unknown media name: " << strSub << std::endl;
					}
				}*/
				contains_MediaName = false;
		}
		break;
		}
	}
	int numDomainIniCondition = Domain_Initial_Conditions.size();
	//if (myid == iMaster) {
	LogFile << std::endl << "Have read " << numDomainIniCondition << " domain initial condition data" << std::endl;
	LogFile << std::endl << "--->>> End reading INDOM" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}

}

//Test of DIFFU
void readDIFFU() {
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	//using namespace Diffusion_Parameters;

	//if (myid == iMaster) {
	//	DIFFU_read = true;
	LogFile << std::endl << "--->>> Start reading DIFFU" << std::endl;
	//}
	int NumPhases = 4, PhaseNum = 4, ComponentNum = 4;
	std::vector<std::string> AryDescribles;
	std::vector<std::vector<double>> diffusivity;
	AryDescribles.push_back("NoDif");
	AryDescribles.push_back("DiffC");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}
	int nIdx = -1;
	int nReadNum = 0;
	std::vector<double> diffusivity_tmp1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		switch (nIdx)
		{
		case 0:
		{
			for (int i = 0; i < NumPhases; i++) diffusivity_tmp1.push_back(0.0);
			diffusivity.push_back(diffusivity_tmp1);
			diffusivity_tmp1.clear();
		}
		break;
		case 1:
		{
			char* p = strtok(const_cast<char*>(strSub.c_str()), ",");
			int j = 0;
			while (p != NULL)
			{
				if (j >= NumPhases) break;
				diffusivity_tmp1.push_back(atof(p));
				p = strtok(NULL, ",");
				j++;
			}
			diffusivity.push_back(diffusivity_tmp1);
			diffusivity_tmp1.clear();
		}
		break;
		}
	}
	int dif_com = diffusivity.size();
	//if (myid == iMaster) {
	LogFile << std::endl << "Have read diffusion coefficients for " << dif_com << " components" << std::endl;
	LogFile << std::endl << "--->>> End reading DIFFU" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl;
	//}

	for (int j = 0; j < PhaseNum; j++) diffusivity_tmp1.push_back(0.0);
	if (dif_com < ComponentNum) {
		for (int j = 0; j < ComponentNum - dif_com; j++) diffusivity.push_back(diffusivity_tmp1);
	}
}

namespace Hydrate_Related_Parameters
{
	// std::vector<std::vector<double>> xhydrsolid;
	std::vector<HydrateGas> hydrGasP;
	int NumHydrCom;
	int N_ThC;
	int N_SpH;
	int N_Rho;
	int EquationOption;
	int eqTableLength;
	int eqSectionNum;
	//int SandFlow_flag;//????????????--> MODEL block
	double P_quad_hyd;//
	double MW_Hydrate;// All of the gas hydrate, not only methane hydrate.
	double InitialHydrateMass;//???????
	double RefDepartureEnth;//???????
	std::vector<double> p_ThC, p_SpH, p_Rho;
	std::vector<double> tempHi, tempLow;
	std::vector<double> eqPressure, eqTemper;
	std::vector<std::vector<double>> acof;
	InhibitorParameters inhibp;

	// Bingham Fluid?????????????
	//bool pressure_gradiant_eff;
	//double EffectDP_fac1;
	//double EffectDP_fac2;

	const double MW_CH4 = 1.604300e1;//
	const double Tc_qup_hyd = 1.0e-2;//
	const double Tk_qup_hyd = 273.16e0;//
	std::string gas_viscosity_option;		// Selection to compute the viscosity of gas.
	std::string gas_cubic_EOS;				// The cubic equation of state.
	bool Lee_Kessler_method;		// Flag indicating whether the Lee-Kessker method is used to compute the departure enthalpy.
	bool noInhib;
}
// ReadHYDRA
//void readHYDRA() {
//	using namespace Hydrate_Related_Parameters;
//
//	//if (myid == iMaster) {
//	//	HYDRA_read = true;
//	LogFile << std::endl << "--->>> Start reading HYDRA" << std::endl;
//	//}
//
//	// Parameters Initialization
//	NumHydrCom = 1;
//	N_ThC = 1;
//	N_SpH = 1;
//	N_Rho = 1;
//	p_ThC = {}, p_SpH = {}, p_Rho = {};
//	EquationOption = 1;
//
//	std::vector<std::string> AryDescribles;
//	AryDescribles.push_back("NumHydrateCom");
//	AryDescribles.push_back("N_ThC");
//	//AryDescribles.push_back("An_ThC");
//	AryDescribles.push_back("N_SpH");
//	//AryDescribles.push_back("Bn_SpH");
//	AryDescribles.push_back("N_Rho");
//	//AryDescribles.push_back("Cn_Rho");
//
//	AryDescribles.push_back("inhibitor_flag");
//	AryDescribles.push_back("Inhibitor_Name");
//	AryDescribles.push_back("Max_TShift");
//	AryDescribles.push_back("Y_atMaxTShift");
//	AryDescribles.push_back("InhibitorMW");
//	AryDescribles.push_back("InhibitorDens");
//	AryDescribles.push_back("InhibitorEnthSol");
//	AryDescribles.push_back("InhibitorCpCoeff");
//	AryDescribles.push_back("EquationOption");
//	AryDescribles.push_back("SandFlow_flag");
//	//AryDescribles.push_back("Reaction_Type");
//
//	std::string strTmp;
//	for (int i = 0; i < AryDescribles.size(); i++)
//	{
//		strTmp = AryDescribles[i];
//		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
//		AryDescribles[i] = strTmp;
//	}
//
//	int nIdx = -1;
//	std::string strSub = "";
//	while (true)
//	{
//		nIdx = -1;
//		strSub = "";
//		int nRes = readSubStr(strSub, nIdx, AryDescribles);
//		if (nRes == -1) break;
//		if (nRes > 0) continue;
//
//		int iTem;
//		double rTem;
//		switch (nIdx)
//		{
//		case 0:
//		{
//			if (strSub.length() > 0) {
//				NumHydrCom = atoi(strSub.c_str());
//				if (NumHydrCom < 1) NumHydrCom = 1;
//				//cout << NumHydrCom << endl;
//				for (int i = 0; i < NumHydrCom; i++) {
//					//getline(MainInputFile, LineText);
//					LineText = readLineText();	// 从输入文件中读取数据。
//					if (!LineText.empty())	// 如果LineText非空，则第一句用来将非空字符前的空格去掉，第二局是将其后面的空格去掉。
//					{
//						LineText.erase(0, LineText.find_first_not_of(" "));
//						LineText.erase(LineText.find_last_not_of(" ") + 1);
//					}
//					//strTmp = LineText;
//					//char* p = strtok(const_cast<char*>(strTmp.c_str()), ",");
//					transform(LineText.begin(), LineText.end(), LineText.begin(), ::toupper);
//					// Method to split input parameters by ","
//					std::vector<std::string> hydratep;
//					hydratep.clear();
//					boost::split(hydratep, LineText, boost::is_any_of(","));
//					struct HydrateGas hydrGasP_Tem;
//					cout << hydratep[0].c_str() << endl;
//
//					hydrGasP_Tem.NameG = hydratep[0].c_str();
//					hydrGasP_Tem.HydrNum = atoi(hydratep[1].c_str());
//					hydrGasP_Tem.MoleFracG = atof(hydratep[2].c_str());
//					hydrGasP.push_back(hydrGasP_Tem);
//				}
//				//cout << hydrGasP[0].HydrNum << endl;
//				//cout << hydrGasP[1].HydrNum << endl;
//				////range check!
//				//Max_NumNRIterations = rangeCheck(iTem, "Maximum number of 0 iterations", 1, 100, 8);
//			}
//		}
//		break;
//		case 1:
//		{
//			if (strSub.length() > 0) {
//				N_ThC = atoi(strSub.c_str());
//				if (N_ThC < 1) N_ThC = 1;
//				LineText = readLineText();
//				if (!LineText.empty())
//				{
//					LineText.erase(0, LineText.find_first_not_of(" "));
//					LineText.erase(LineText.find_last_not_of(" ") + 1);
//				}
//				std::vector<std::string> p;
//				p.clear();
//				boost::split(p, LineText, boost::is_any_of(","));
//				for (int i = 0; i < N_ThC; i++) {
//					p_ThC.push_back(atof(p[i].c_str()));
//				}
//				cout << p_ThC[0] << endl;
//				cout << p_ThC[1] << endl;
//				cout << p_ThC[2] << endl;
//			}
//		}
//		break;
//		case 2:
//		{
//			if (strSub.length() > 0) {
//				N_SpH = atoi(strSub.c_str());
//				if (N_SpH < 1) N_SpH = 1;
//				LineText = readLineText();
//				if (!LineText.empty())
//				{
//					LineText.erase(0, LineText.find_first_not_of(" "));
//					LineText.erase(LineText.find_last_not_of(" ") + 1);
//				}
//				std::vector<std::string> p;
//				p.clear();
//				boost::split(p, LineText, boost::is_any_of(","));
//				for (int i = 0; i < N_SpH; i++) {
//					p_SpH.push_back(atof(p[i].c_str()));
//				}
//				cout << p_SpH[0] << endl;
//				cout << p_SpH[1] << endl;
//				cout << p_SpH[2] << endl;
//			}
//		}
//		break;
//		case 3:
//		{
//			if (strSub.length() > 0) {
//				N_Rho = atoi(strSub.c_str());
//				if (N_Rho < 1) N_Rho = 1;
//				LineText = readLineText();
//				if (!LineText.empty())
//				{
//					LineText.erase(0, LineText.find_first_not_of(" "));
//					LineText.erase(LineText.find_last_not_of(" ") + 1);
//				}
//				std::vector<std::string> p;
//				p.clear();
//				boost::split(p, LineText, boost::is_any_of(","));
//				for (int i = 0; i < N_Rho; i++) {
//					p_Rho.push_back(atof(p[i].c_str()));
//				}
//				cout << p_Rho[0] << endl;
//				cout << p_Rho[1] << endl;
//				cout << p_Rho[2] << endl;
//			}
//		}
//		break;
//		case 4:
//		{
//			if (strSub.length() > 0) {// 0 (false) or 1 (true)
//				int p = atoi(strSub.c_str());
//				if (p == 0 || p == 1) inhibp.inhibitor_flag = p;
//				else inhibp.inhibitor_flag = 0;		// default is 0, means inhibitor does not exist.
//
//				//Flow::fs.noinhib = false????????????
//			}
//		}
//		break;
//		case 5:
//		{
//			inhibp.Inhibitor_Name = "";
//			if (strSub.length() > 0) inhibp.Inhibitor_Name = strSub.c_str();
//			cout << inhibp.Inhibitor_Name << endl;
//		}
//		break;
//		case 6:
//		{
//			if (strSub.length() > 0) inhibp.Max_TShift = atof(strSub.c_str());
//			if (inhibp.Max_TShift < 0) inhibp.Max_TShift = 0.0;
//		}
//		break;
//		case 7:
//		{
//			if (strSub.length() > 0) inhibp.Y_atMax_TShift = atof(strSub.c_str());
//			if (inhibp.Y_atMax_TShift < 0) inhibp.Y_atMax_TShift = 0.0;
//		}
//		break;
//		case 8:
//		{
//			if (strSub.length() > 0) inhibp.InhibMW = atof(strSub.c_str());
//			if (inhibp.InhibMW < 0) inhibp.InhibMW = 0.0;
//		}
//		break;
//		case 9:
//		{
//			if (strSub.length() > 0) inhibp.InhibDens = atof(strSub.c_str());
//			if (inhibp.InhibDens < 0) inhibp.InhibDens = 0.0;
//		}
//		break;
//		case 10:
//		{
//			if (strSub.length() > 0) inhibp.InhibEnthSol = atof(strSub.c_str());
//			if (inhibp.InhibEnthSol < 0) inhibp.InhibEnthSol = 0.0;
//		}
//		break;
//		case 11:
//		{
//			if (strSub.length() > 0) {
//				std::vector<std::string> p;
//				p.clear();
//				boost::split(p, strSub, boost::is_any_of(","));
//				if (p.size() == 3) {
//					for (int i = 0; i < 3; i++) {
//						inhibp.InhibCpCoeff[i] = atof(p[i].c_str());
//						if (inhibp.InhibCpCoeff[i] < 0) inhibp.InhibCpCoeff[i] = 0.0;
//					}
//				}
//
//				cout << inhibp.InhibCpCoeff[0] << endl;
//				cout << inhibp.InhibCpCoeff[1] << endl;
//				cout << inhibp.InhibCpCoeff[2] << endl;
//
//
//
//				////inhibp.InhibCpCoeff = atof(strSub.c_str());
//				////if (inhibp.InhibCpCoeff < 0) inhibp.InhibCpCoeff = 0.0;
//			}
//		}
//		break;
//		case 12:
//		{
//			if (strSub.length() > 0) EquationOption = atoi(strSub.c_str());
//			if (EquationOption < 0) EquationOption = 2;		// 0 for Moridis, 1 for Kamath, 2 for combining two methods, ??????? or other methods (查表等).
//		}
//		break;
//		case 13:
//		{
//			if (strSub.length() > 0) {// 0 (false) or 1 (true)
//				int p = atoi(strSub.c_str());
//				if (p == 0 || p == 1) SandFlow_flag = p;
//				else SandFlow_flag = 0;		// default is 0, means inhibitor does not exist.
//
//				//Flow::fs.noSand = false
//			}
//		}
//		break;
//		default:
//			break;
//		}
//	}
//
//	// Salt_Present?????????????
//
//}

void readHYDRA()
{
	//using namespace General_Control_Parameters;
	//using namespace Basic_Parameters;
	//using namespace General_External_File_Units;
	using namespace Hydrate_Related_Parameters;
	//using namespace Diffusion_Parameters;

	//if (myid == iMaster) {
	//	HYDRA_read = true;
	LogFile << "--->>> Start reading HYDRA_" << std::endl << std::endl;
	//}

	// Parameters Initialization
	NumHydrCom = 1;
	N_ThC = 1;
	N_SpH = 1;
	N_Rho = 1;
	noInhib = true;
	inhibp.Inhibitor_Name = " ";
	inhibp.Max_TShift = 0.0;
	inhibp.Y_atMax_TShift = 0.0;
	inhibp.InhibMW = 0.0;
	inhibp.InhibDens = 0.0;
	inhibp.InhibEnthSol = 0.0;
	for (int i = 0; i < 3; i++) inhibp.InhibCpCoeff[i] = 0.0;
	EquationOption = 2;
	int nSub;

	//p_ThC = {}, p_SpH = {}, p_Rho = {};

	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("NumHydrateCom");
	AryDescribles.push_back("N_ThC");
	AryDescribles.push_back("N_SpH");
	AryDescribles.push_back("N_Rho");
	AryDescribles.push_back("Has_Inhibitor");
	AryDescribles.push_back("Inhibitor_Name");
	AryDescribles.push_back("Max_TShift");
	AryDescribles.push_back("Y_atMaxTShift");
	AryDescribles.push_back("InhibitorMW");
	AryDescribles.push_back("InhibitorDens");
	AryDescribles.push_back("InhibitorEnthSol");
	AryDescribles.push_back("InhibitorCpCoeff");
	AryDescribles.push_back("EquationOption");
	//AryDescribles.push_back("SandFlow_flag");
	//AryDescribles.push_back("Reaction_Type");

	std::string strTmp;
	for (int i = 0; i < AryDescribles.size(); i++)
	{
		strTmp = AryDescribles[i];
		transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);
		AryDescribles[i] = strTmp;
	}

	int nIdx = -1;
	std::string strSub = "";
	while (true)
	{
		nIdx = -1;
		strSub = "";
		int nRes = readSubStr(strSub, nIdx, AryDescribles);
		if (nRes == -1) break;
		if (nRes > 0) continue;

		int iTem;
		double rTem;
		switch (nIdx)
		{
		case 0:
		{
			if (strSub.length() > 0) {
				NumHydrCom = atoi(strSub.c_str());
				if (NumHydrCom < 1) NumHydrCom = 1;

				/* if NumHydrCom > 1, the thing must be clear that the number of hydrate components must be corresponding to the number of gas in GasLaw in rGasPar...
					which means: 0 for CH4, 2/3(!noAqu2) for CO2, 3/4(!noAqu2) for N2...
					Tips: Line 289-292, EOS_of_NCG cpp file.
				*/

				for (int i = 0; i < NumHydrCom; i++) {
					LineText = readLineText();
					transform(LineText.begin(), LineText.end(), LineText.begin(), ::toupper);
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					std::vector<std::string> hydratep;
					hydratep.clear();
					boost::split(hydratep, LineText, boost::is_any_of(","));
					if (hydratep.size() < 3) {
						LogFile << "Hydrate Name, HydrNum or MoleFracG is/are missing, Please input all of them." << std::endl;
						exit(0);
					}
					struct HydrateGas hydrGasP_Tem;
					hydrGasP_Tem.NameG = hydratep[0].c_str();
					hydrGasP_Tem.HydrNum = atof(hydratep[1].c_str());
					hydrGasP_Tem.MoleFracG = atof(hydratep[2].c_str());
					hydrGasP.push_back(hydrGasP_Tem);
				}
			}
		}
		break;
		case 1:
		{
			if (strSub.length() > 0) {
				N_ThC = atoi(strSub.c_str());
				if (N_ThC < 1) N_ThC = 1;
				LineText = readLineText();
				if (!LineText.empty())
				{
					LineText.erase(0, LineText.find_first_not_of(" "));
					nSub = LineText.find("/");
					if (nSub != -1) LineText.erase(nSub);
					LineText.erase(LineText.find_last_not_of(" ") + 1);
				}
				std::vector<std::string> p;
				p.clear();
				boost::split(p, LineText, boost::is_any_of(","));
				if (p.size() != N_ThC) {
					LogFile << "The number of p_ThC is " << N_ThC << ", please input all of them." << std::endl;
					exit(0);
				}
				for (int i = 0; i < N_ThC; i++) {
					p_ThC.push_back(atof(p[i].c_str()));
				}
			}
		}
		break;
		case 2:
		{
			if (strSub.length() > 0) {
				N_SpH = atoi(strSub.c_str());
				if (N_SpH < 1) N_SpH = 1;
				LineText = readLineText();
				if (!LineText.empty())
				{
					LineText.erase(0, LineText.find_first_not_of(" "));
					nSub = LineText.find("/");
					if (nSub != -1) LineText.erase(nSub);
					LineText.erase(LineText.find_last_not_of(" ") + 1);
				}
				std::vector<std::string> p;
				p.clear();
				boost::split(p, LineText, boost::is_any_of(","));
				if (p.size() != N_SpH) {
					LogFile << "The number of p_SpH is " << N_SpH << ", please input all of them." << std::endl;
					exit(0);
				}
				for (int i = 0; i < N_SpH; i++) {
					p_SpH.push_back(atof(p[i].c_str()));
				}
			}
		}
		break;
		case 3:
		{
			if (strSub.length() > 0) {
				N_Rho = atoi(strSub.c_str());
				if (N_Rho <= 0) N_Rho = 0;
				if (N_Rho > 0) {
					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					std::vector<std::string> p;
					p.clear();
					boost::split(p, LineText, boost::is_any_of(","));
					if (p.size() != N_Rho) {
						LogFile << "The number of p_Rho is " << N_Rho << ", please input all of them." << std::endl;
						exit(0);
					}
					for (int i = 0; i < N_Rho; i++) {
						p_Rho.push_back(atof(p[i].c_str()));
					}
				}
			}
		}
		break;
		case 4:
		{
			if (strSub.length() > 0) {// 0 (false) or 1 (true)
				if (trim(strSub) == "TRUE" || trim(strSub) == "true" || trim(strSub) == "True") { noInhib = false; }
			}
		}
		break;
		case 5:
		{
			if (strSub.length() > 0) inhibp.Inhibitor_Name = strSub.c_str();
		}
		break;
		case 6:
		{
			if (strSub.length() > 0) inhibp.Max_TShift = atof(strSub.c_str());
			if (inhibp.Max_TShift < 0) inhibp.Max_TShift = 0.0;
		}
		break;
		case 7:
		{
			if (strSub.length() > 0) inhibp.Y_atMax_TShift = atof(strSub.c_str());
			if (inhibp.Y_atMax_TShift < 0) inhibp.Y_atMax_TShift = 0.0;
		}
		break;
		case 8:
		{
			if (strSub.length() > 0) inhibp.InhibMW = atof(strSub.c_str());
			if (inhibp.InhibMW < 0) inhibp.InhibMW = 0.0;
		}
		break;
		case 9:
		{
			if (strSub.length() > 0) inhibp.InhibDens = atof(strSub.c_str());
			if (inhibp.InhibDens < 0) inhibp.InhibDens = 0.0;
		}
		break;
		case 10:
		{
			if (strSub.length() > 0) inhibp.InhibEnthSol = atof(strSub.c_str());
			if (inhibp.InhibEnthSol < 0) inhibp.InhibEnthSol = 0.0;
		}
		break;
		case 11:
		{
			if (strSub.length() > 0) {
				std::vector<std::string> p;
				p.clear();
				boost::split(p, strSub, boost::is_any_of(","));
				if (p.size() != 3) {
					LogFile << "The number of inhibitor coefficients should be 3, please input all of them." << std::endl;
					exit(0);
				}
				else {
					for (int i = 0; i < 3; i++) {
						inhibp.InhibCpCoeff[i] = atof(p[i].c_str());
						if (inhibp.InhibCpCoeff[i] < 0) inhibp.InhibCpCoeff[i] = 0.0;
					}
				}
			}
		}
		break;
		case 12:
		{
			// 1 for Kamath, 0/2 for Moridis, 3 for table interpolation, 4 for user provided functions.
			if (strSub.length() > 0) EquationOption = atoi(strSub.c_str());

			if (EquationOption < 0 || EquationOption > 4) EquationOption = 2;

			else if (EquationOption == 3) {
				eqTableLength = 0;
				LineText = readLineText();
				if (!LineText.empty())
				{
					LineText.erase(0, LineText.find_first_not_of(" "));
					nSub = LineText.find("/");
					if (nSub != -1) LineText.erase(nSub);
					LineText.erase(LineText.find_last_not_of(" ") + 1);
				}
				eqTableLength = atoi(LineText.c_str());
				if (eqTableLength < 5) {
					EquationOption = 0;
					LogFile << "The Moridis [2003] will be used." << endl;
				}
				else {
					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					std::vector<std::string> p;
					p.clear();
					boost::split(p, LineText, boost::is_any_of(","));
					if (p.size() != eqTableLength) {
						LogFile << "The eqTableLength must be " << eqTableLength << ", please input all of them." << std::endl;
						exit(0);
					}
					for (int i = 0; i < eqTableLength; i++) {
						eqPressure.push_back(atof(p[i].c_str()));
					}

					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					p.clear();
					boost::split(p, LineText, boost::is_any_of(","));

					if (p.size() != eqTableLength) {
						LogFile << "The eqTableLength must be " << eqTableLength << ", please input all of them." << std::endl;
						exit(0);
					}
					for (int i = 0; i < eqTableLength; i++) {
						eqTemper.push_back(atof(p[i].c_str()));
					}
				}			
			}
			else if (EquationOption == 4) {
				eqSectionNum = 0;
				LineText = readLineText();
				if (!LineText.empty())
				{
					LineText.erase(0, LineText.find_first_not_of(" "));
					nSub = LineText.find("/");
					if (nSub != -1) LineText.erase(nSub);
					LineText.erase(LineText.find_last_not_of(" ") + 1);
				}
				eqSectionNum = atoi(LineText.c_str());
				if (eqSectionNum < 1) {
					EquationOption = 0;
					LogFile << "The Moridis [2003] will be used." << endl;
				}
				else if (eqSectionNum > 10) {
					LogFile << "The maximum allowed function section number is 10, stop simulation!" << endl;
					exit(0);
				}
				else {
					// Input tempLow data.
					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					std::vector<std::string> p;
					p.clear();
					boost::split(p, LineText, boost::is_any_of(","));
					if (p.size() != eqSectionNum) {
						LogFile << "The eqSectionNum must be " << eqSectionNum << ", please input all of them." << std::endl;
						exit(0);
					}
					for (int i = 0; i < eqSectionNum; i++) {
						tempLow.push_back(atof(p[i].c_str()));
					}

					// Input tempHi data.
					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					p.clear();
					boost::split(p, LineText, boost::is_any_of(","));
					if (p.size() != eqSectionNum) {
						LogFile << "The eqSectionNum must be " << eqSectionNum << ", please input all of them." << std::endl;
						exit(0);
					}
					for (int i = 0; i < eqSectionNum; i++) {
						tempHi.push_back(atof(p[i].c_str()));
					}

					// Input acof data.
					LineText = readLineText();
					if (!LineText.empty())
					{
						LineText.erase(0, LineText.find_first_not_of(" "));
						nSub = LineText.find("/");
						if (nSub != -1) LineText.erase(nSub);
						LineText.erase(LineText.find_last_not_of(" ") + 1);
					}
					p.clear();
					boost::split(p, LineText, boost::is_any_of(";"));
					if (p.size() != eqSectionNum) {
						LogFile << "The eqSectionNum must be " << eqSectionNum << ", please input all of them." << std::endl;
						exit(0);
					}
					std::vector<double> tem_acof;
					for (int i = 0; i < eqSectionNum; i++) {
						std::vector<std::string> q;
						q.clear();
						boost::split(q, p[i], boost::is_any_of(","));
						if (q.size() != 7) {
							LogFile << "The row of acof data must be 7, please input all of them." << std::endl;
							exit(0);
						}
						for (int j = 0; j < 7; j++) {
							tem_acof.push_back(atof(q[j].c_str()));
						}
						acof.push_back(tem_acof);
						tem_acof.clear();
					}
				}
			}
		}
		break;
		//case 13:
		//{
		//	SandFlow_flag = 0;		// default is 0, means inhibitor does not exist.
		//	if (strSub.length() > 0) {// 0 (false) or 1 (true)
		//		int p = atoi(strSub.c_str());
		//		if (p == 0 || p == 1) SandFlow_flag = p;

		//		//Flow::fs.noSand = false
		//	}
		//}
		//break;
		default:
			break;
		}

		// Salt_Present?????????????
	}
	//if (myid == iMaster) {
	LogFile << "--->>> End reading HYDRA_" << std::endl << std::endl;
	LogFile << "===============================================================" << std::endl << std::endl;
	//}}
}

void ReadFormattedTextLine(std::string TextLine, std::string fmt, std::string strs[], double dbls[], int ints[]) {
	//using namespace General_External_File_Units;
	unsigned int CurrentLoc = 0;
	unsigned int CurrentD = 0, CurrentS = 0, CurrentI = 0;
	std::size_t stringLength = TextLine.length();

	for (std::string::iterator it = fmt.begin(); it != fmt.end(); ++it) {

		if (*it == 'D') {
			++it;
			int i10 = (*it - '0') * 10;
			++it;
			int i0 = (*it - '0') + i10;

			if ((CurrentLoc + i0) > stringLength) TextLine.append(i0, ' ');
			std::string TemString = TextLine.substr(CurrentLoc, i0);
			std::string SpaceString;
			SpaceString.append(i0, ' ');

			int item = 0;
			std::size_t found;
			do {
				found = TemString.find_first_not_of("0123456789+-Ee. ");
				if (found != std::string::npos) {
					int item = int(TemString.at(found));
					if (item == 13) TemString.replace(found, 1, 1, ' ');
					else if (item == 9) TemString.replace(found, 1, 1, ' ');
					else {
						//LogFile << "error in reading file '" << readingFileName << "' at the text line: " << readingLineNumber;
						LogFile << " for reading expected double variable at: '" << TemString << "'" << std::endl;
						std::cout << "Abnormal stop, please see the log file for details." << std::endl;
						exit(1);
					}
				}
			} while (found != std::string::npos);
			if (TemString == SpaceString) {
				dbls[CurrentD] = 0.0;
			}
			else {
				dbls[CurrentD] = std::stod(TemString);
			}

			CurrentD = CurrentD + 1;
			CurrentLoc = CurrentLoc + i0;
		}
		else if (*it == 'S') {
			++it;
			int i10 = (*it - '0') * 10;
			++it;
			int i0 = (*it - '0') + i10;
			if ((CurrentLoc + i0) > stringLength) {
				TextLine.append(CurrentLoc + i0 - stringLength, ' ');
			}
			strs[CurrentS] = TextLine.substr(CurrentLoc, i0);

			CurrentS = CurrentS + 1;
			CurrentLoc = CurrentLoc + i0;
		}
		else if (*it == 'I') {
			++it;
			int i10 = (*it - '0') * 10;
			++it;
			int i0 = (*it - '0') + i10;
			if ((CurrentLoc + i0) > stringLength) TextLine.append(i0, ' ');
			std::string TemString = TextLine.substr(CurrentLoc, i0);
			std::string SpaceString;
			SpaceString.append(i0, ' ');
			size_t found;
			do {
				found = TemString.find_first_not_of("0123456789+- ");
				if (found != std::string::npos) {
					int item = int(TemString.at(found));
					if (item == 13) TemString.replace(found, 1, 1, ' ');
					else if (item == 9) TemString.replace(found, 1, 1, ' ');
					else {
						//LogFile << "error in reading file '" << readingFileName << "' at the text line: " << readingLineNumber;
						LogFile << ", for reading expected integer variable at: '" << TemString << "'" << std::endl;
						std::cout << "Abnormal stop, please see the log file for details." << std::endl;
						exit(1);
					}
				}
			} while (found != std::string::npos);

			if (TemString == SpaceString) {
				ints[CurrentI] = 0;
			}
			else {
				ints[CurrentI] = std::stoi(TemString);
			}
			CurrentI = CurrentI + 1;
			CurrentLoc = CurrentLoc + i0;
		}
	}

}

char* ReadGetLine(char* szLine, int nMaxCount, FILE* pFile)
{
	char* pResult = fgets(szLine, nMaxCount, pFile);
	if (pResult != NULL)
	{
		char* find = strchr(szLine, '\n');
		if (find != NULL) *find = '\0';
		return pResult;
	}
	memset(szLine, 0, sizeof(char) * nMaxCount);
	return NULL;
}

// readHYDRA for read_main_inputs.cpp/hpp
void readHYDRATE() {

}

// Function for calculating component hydrate mass fraction of gas and H2O.
void HydrateMW(const int& GasMW, const int& HydrNum, double& HydrMW, double& MassFracGas, double& MassFracH2O) {
	//T+H 922-966
	double H2OMW = 18.0;
	HydrMW = GasMW + H2OMW * HydrNum;
	MassFracGas = GasMW / HydrMW;
	MassFracH2O = 1.0e0 - GasMW;
	//...... Convert to g / mol
	HydrMW = 1.0e3 * HydrMW;
}

//void readIinitialConditionFile() {
//
//	std::string fileName;
//	std::string fileExtName;
//	std::string keyWord;
//
//	std::vector<int> workingList;
//	double dbls[8];
//	//double dbls[10];		// hydrate parameters input by Bingbo.
//	std::string strs[8];
//	int ints[8];
//	char szLine[128];
//	int Incon_Portions_Number = 1;
//	int readingLineNumber = 0;
//
//	fileName = "INCON";
//	std::string readingFileName = fileName;
//	//workingList.push_back(mpi_size - 1);
//	FILE* ReadInconFile;
//	ReadInconFile = fopen(fileName.c_str(), "r");
//	if (ReadInconFile == NULL)
//	{
//		LogFile << "Open File: " << fileName << " error. You may need to add keyword 'INCON' in the main input." << std::endl;
//		exit(1);
//	}
//
//	ReadGetLine(szLine, 128, ReadInconFile);
//	readingLineNumber++;
//	ReadFormattedTextLine(szLine, "S05S05S05I05", strs, dbls, ints);
//
//
//	keyWord = strs[0];
//	if (strs[2] == "SAVE=" || strs[2] == "save=" || strs[2] == "Save=") {
//		if (ints[0] > 1) {
//			Incon_Portions_Number = ints[0];
//		}
//	}
//
//	int iTem, myNum;
//	pc.GetReduce(Incon_Portions_Number, iTem, "MAX");
//	Incon_Portions_Number = iTem;
//	workingList.resize(Incon_Portions_Number);
//	for (int i = 0; i < Incon_Portions_Number; i++)
//		workingList[i] = mpi_size - 1 - i;
//	bool inList = false;
//	for (int i = 0; i < Incon_Portions_Number; i++) {
//		if (myid == workingList[i]) {
//			inList = true;
//			myNum = i;
//		}
//	}
//
//	std::vector<std::string> elemName;
//	std::vector<double>  por, perm1, perm2, perm3;
//	std::vector<std::string> stateIdx;
//	std::vector < std::vector <double>> primaryVar;
//	//std::vector < std::vector <double>> hydrateVar;		// hydrate parameters input by Bingbo.
//	std::vector <double> temVar;
//	temVar.resize(NumCom + 1);
//	//temVar.resize(NumCom + 4);		// hydrate parameters input by Bingbo.
//	int i1 = Tot_NumTimeSteps;
//	int i2 = Tot_NumNRIterations;
//	double d1 = InitialTimeStep;
//	double d2 = TimeShift;
//
//	if (inList) {
//		long long int item1 = myNum;
//		fileExtName = "." + std::to_string(item1);
//
//		if (Incon_Portions_Number <= mpi_size) {
//			if (myid != workingList[0]) {
//				fileName = fileName + fileExtName;
//				ReadInconFile = fopen(fileName.c_str(), "r");
//				if (ReadInconFile == NULL)
//				{
//					LogFile << "Open File: " << fileName << " error. " << std::endl;
//					pc.StopAll();
//				}
//			}
//		}
//		else {
//			if (myid == workingList[0]) {
//				LogFile << "The INCON file seperated into too many files, please use more cores/cpus to run the simulation or ",
//					LogFile << "seperate the INCON with less files.\n";
//			}
//			pc.StopAll();
//		}
//
//		if (myid != workingList[0]) { ReadGetLine(szLine, mnReadMaxLen, ReadInconFile);; readingLineNumber++; }  //read the first line
//		std::string fmt = "S" + ElemNaLth[ElemNameLength] + "S" + ElemNaLth[15 - ElemNameLength] + "D15D15D15D15S03";
//		std::string TemString;
//		std::vector<std::string>Line_Split_tmp;
//		std::string LineText;
//		int NumComP1 = NumCom + 1;
//		while (ReadGetLine(szLine, mnReadMaxLen, ReadInconFile) != NULL) {
//			readingLineNumber++;
//			ReadFormattedTextLine(szLine, fmt, strs, dbls, ints);
//			TemString = strs[0].substr(0, 3);
//			if (TemString == "   ") break;
//			if (TemString == "+++") break;
//
//			elemName.push_back(strs[0]);
//			stateIdx.push_back(strs[2]);
//			por.push_back(dbls[0]);
//			perm1.push_back(dbls[1]);
//			perm2.push_back(dbls[2]);
//			perm3.push_back(dbls[3]);
//
//			ReadGetLine(szLine, 256, ReadInconFile);
//			readingLineNumber++;
//			Line_Split_tmp.clear();
//			LineText = szLine;
//			cf.SplitString(LineText, Line_Split_tmp, " ");
//			int item = Line_Split_tmp.size();
//			if (item > NumComP1) item = NumComP1;
//			if (item < NumComP1) {
//				LogFile << "The primary variables provided in file INCON is not enough at line: " << readingLineNumber << std::endl;
//				LogFile << "The number of primary varible must equals to Component_Number + 1 (including temperature), stop simulation!" << std::endl;
//				pc.StopAll();
//			}
//
//			for (int i = 0; i < item; i++) temVar[i] = cf.stringToNum<double>(Line_Split_tmp[i], readingFileName, readingLineNumber);
//			primaryVar.push_back(temVar);
//		}
//		if (TemString == "+++") {
//			ReadGetLine(szLine, mnReadMaxLen, ReadInconFile);
//			readingLineNumber++;
//			ReadFormattedTextLine(szLine, "I10I10S05D20D20", strs, dbls, ints);
//			i1 = ints[0];
//			i2 = ints[1];
//			d1 = dbls[0];
//			d2 = dbls[1];
//			if (d1 <= 0.0) d1 = 1000.0;
//		}
//		fclose(ReadInconFile);
//	}                                    // end of  if (inList)
//	pc.GetReduce(i1, Tot_NumTimeSteps, "MAX");
//	pc.GetReduce(i2, Tot_NumNRIterations, "MAX");
//	pc.GetReduce(d1, InitialTimeStep, "MAX");
//	pc.GetReduce(d2, TimeShift, "MAX");
//
//	std::vector<std::string>temElemName;
//	std::vector<bool> notDone;
//	std::vector<int> sendToLocalNum;
//	std::vector<int> sendToLocalIdx;
//
//	int localGridNum = ModelMesh.elem.size();
//	notDone.resize(localGridNum, true);
//	sendToLocalNum.resize(mpi_size);
//
//	ElemState.resize(localGridNum);
//	ElemMedia1.resize(localGridNum);
//	for (int j = 0; j < localGridNum; j++) {
//		ElemMedia1[j].porosity = 0.0; ElemMedia1[j].perm[0] = 0.0;
//		ElemMedia1[j].perm[1] = 0.0; ElemMedia1[j].perm[2] = 0.0;
//	}
//	xk.resize(localGridNum, std::vector<double>(NumCom + 1, 0.0));
//	for (int j = 0; j < localGridNum; j++) ElemState[j].index1 = 100;
//
//	for (int i = 0; i < Incon_Portions_Number; i++) {
//		for (int j = 0; j < mpi_size; j++) sendToLocalNum[j] = 0;
//
//		int portionInconSize = elemName.size();
//		pc.BroadcastDatum(portionInconSize, workingList[i]);
//
//		temElemName.resize(portionInconSize);
//		if (myid == workingList[i]) {
//			sendToLocalIdx.resize(portionInconSize, -1);
//			for (int j = 0; j < portionInconSize; j++) temElemName[j] = elemName[j];
//		}
//		pc.BroadcastStrVector(temElemName, portionInconSize, workingList[i]);
//
//		std::vector<int> index;
//		std::vector<int> localIdx;
//		std::vector<int> localLocation;
//
//		for (int j = 0; j <= portionInconSize; j++)
//			index.push_back(j);
//		temElemName.push_back("  ");
//		for (int j = portionInconSize; j > 0; j--) temElemName[j] = temElemName[j - 1];
//		temElemName[0] = "   ";
//		cf.quickSort(portionInconSize, temElemName, index);
//		if (portionInconSize == 0) temElemName.push_back("  ");
//
//		int iJ0 = 1;
//		for (int j = 0; j < localGridNum; j++) {
//			int iJ = -1;
//			if (ModelMesh.elem[j].name != temElemName[iJ0]) {
//				if (notDone[j]) cf.binarySearch(portionInconSize, temElemName, index, ModelMesh.elem[j].name, iJ);
//			}
//			else iJ = iJ0;
//			if (iJ > 0) {
//				localIdx.push_back(iJ - 1);
//				localLocation.push_back(j);
//				notDone[j] = false;
//				if (iJ < portionInconSize) iJ0 = iJ + 1;
//			}
//		}
//		//
//
//		int localReceiveN = localIdx.size();
//		index.clear();
//		std::vector<int> TemIntV;
//		for (int j = 0; j <= localReceiveN; j++) index.push_back(j);
//		TemIntV.push_back(-1);
//		for (int j = 0; j < localReceiveN; j++) TemIntV.push_back(localIdx[j]);
//		cf.quickSort_I(localReceiveN, TemIntV, index);
//		for (int j = 0; j < localReceiveN; j++) localIdx[j] = TemIntV[index[j + 1]];
//		for (int j = 0; j < localReceiveN; j++) TemIntV[j] = localLocation[j];
//		for (int j = 0; j < localReceiveN; j++) localLocation[j] = TemIntV[index[j + 1] - 1];
//		TemIntV.clear();
//		//
//		if (myid != workingList[i]) pc.NBSendI(localReceiveN, workingList[i], workingList[i] * 10000 + myid);
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++)
//				if (myid != j) pc.NBReceiveI(sendToLocalNum[j], j, myid * 10000 + j);
//			sendToLocalNum[workingList[i]] = localReceiveN;
//		}
//		if (myid != workingList[i]) {
//			if (localReceiveN > 0) pc.NBSendIntVector(localIdx, localReceiveN, workingList[i], workingList[i] * 10000 + myid + TenM);
//		}
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (myid != j) {
//					std::vector<int> temIdx;
//					if (sendToLocalNum[j] > 0) {
//						temIdx.resize(sendToLocalNum[j], -1);
//						pc.NBReceiveIntVector(temIdx, sendToLocalNum[j], j, myid * 10000 + j + TenM);
//						for (int k = 0; k < sendToLocalNum[j]; k++)  sendToLocalIdx[temIdx[k]] = j;
//					}
//				}
//				else {
//					for (int k = 0; k < localReceiveN; k++) sendToLocalIdx[localIdx[k]] = j;
//				}
//			}
//		}
//		index.clear();
//		localIdx.clear();
//		temElemName.clear();
//
//		ModelParallelWork mpw;
//		int sendN;
//		int tag;
//		std::vector<std::string> workArray;
//
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (workingList[i] == j) continue;
//				sendN = sendToLocalNum[j];
//				tag = workingList[i] * 10000 + j;
//
//				if (sendN > 0) mpw.SendLocalStrVariables(stateIdx, sendToLocalIdx, tag, sendN, j);
//			}
//			for (int k = 0; k < portionInconSize; k++)
//				if (sendToLocalIdx[k] == myid)
//					workArray.push_back(stateIdx[k]);
//		}
//		else {
//			tag = workingList[i] * 10000 + myid;
//			if (localReceiveN > 0) mpw.ReceiveLocalStrVariables(workArray, 3, tag, localReceiveN, workingList[i]);
//		}
//		for (int j = 0; j < localReceiveN; j++) {
//			ElemState[localLocation[j]].index1 = 6;
//			for (int k = 0; k < fs.NumberOfStates; k++) {
//				if (workArray[j] == fs.StatName[k]) ElemState[localLocation[j]].index1 = k + 1;
//			}
//		}
//		workArray.clear();
//
//		std::vector<double> workArray1;
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (workingList[i] == j) continue;
//				sendN = sendToLocalNum[j];
//				tag = workingList[i] * 10000 + j + TenM;
//				if (sendN > 0) mpw.SendLocalDblVariables(perm1, sendToLocalIdx, tag, sendN, j);
//			}
//			for (int k = 0; k < portionInconSize; k++)
//				if (sendToLocalIdx[k] == myid) workArray1.push_back(perm1[k]);
//		}
//		else {
//			tag = workingList[i] * 10000 + myid + TenM;
//			if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//		}
//
//		for (int j = 0; j < localReceiveN; j++) {
//			ElemMedia1[localLocation[j]].perm[0] = workArray1[j];
//		}
//		workArray1.clear();
//		//
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (workingList[i] == j) continue;
//				sendN = sendToLocalNum[j];
//				tag = workingList[i] * 10000 + j + TenM * 2;
//				if (sendN > 0) mpw.SendLocalDblVariables(perm2, sendToLocalIdx, tag, sendN, j);
//			}
//			for (int k = 0; k < portionInconSize; k++)
//				if (sendToLocalIdx[k] == myid) workArray1.push_back(perm2[k]);
//		}
//		else {
//			tag = workingList[i] * 10000 + myid + TenM * 2;
//			if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//		}
//
//		for (int j = 0; j < localReceiveN; j++) {
//			ElemMedia1[localLocation[j]].perm[1] = workArray1[j];
//		}
//		workArray1.clear();
//		//
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (workingList[i] == j) continue;
//				sendN = sendToLocalNum[j];
//				tag = workingList[i] * 10000 + j + TenM * 3;
//				if (sendN > 0) mpw.SendLocalDblVariables(perm3, sendToLocalIdx, tag, sendN, j);
//			}
//			for (int k = 0; k < portionInconSize; k++)
//				if (sendToLocalIdx[k] == myid) workArray1.push_back(perm3[k]);
//		}
//		else {
//			tag = workingList[i] * 10000 + myid + TenM * 3;
//			if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//		}
//
//		for (int j = 0; j < localReceiveN; j++) {
//			ElemMedia1[localLocation[j]].perm[2] = workArray1[j];
//		}
//		workArray1.clear();
//		//
//		if (myid == workingList[i]) {
//			for (int j = 0; j < mpi_size; j++) {
//				if (workingList[i] == j) continue;
//				sendN = sendToLocalNum[j];
//				tag = workingList[i] * 10000 + j + TenM * 4;
//				if (sendN > 0) mpw.SendLocalDblVariables(por, sendToLocalIdx, tag, sendN, j);
//			}
//			for (int k = 0; k < portionInconSize; k++)
//				if (sendToLocalIdx[k] == myid) workArray1.push_back(por[k]);
//		}
//		else {
//			tag = workingList[i] * 10000 + myid + TenM * 4;
//			if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//		}
//
//		for (int j = 0; j < localReceiveN; j++) {
//			ElemMedia1[localLocation[j]].porosity = workArray1[j];
//		}
//		workArray1.clear();
//		//
//		std::vector<double> rTem;
//		for (int l = 0; l < NumCom + 1; l++) {
//			if (myid == workingList[i]) {
//				// for (int j = 0; j < mpi_size; j++) {
//				//	 if (workingList[i] == j) continue;
//				//	 rTem.resize(portionInconSize);
//				//	 for (int k = 0; k < portionInconSize; k++) rTem[k] = primaryVar[k][l];
//				//	 sendN = sendToLocalNum[j];
//				//	 tag = workingList[i] * 10000 + j + TenM * (5 + l);
//				//	 if (sendN > 0) {
//				//		 mpw.SendLocalDblVariables(rTem, sendToLocalIdx, tag, sendN, j);
//				//	 }
//				// }
//				// for (int k = 0; k < portionInconSize; k++) {
//				//	 if (sendToLocalIdx[k] == myid) workArray1.push_back(primaryVar[k][l]);
//				// }
//			 //}
//			 //else {
//				// tag = workingList[i] * 10000 + myid + TenM * (5 + l);
//				// if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//			 //}
//
//			 //for (int j = 0; j < localReceiveN; j++) {
//				// xk[localLocation[j]][l] = workArray1[j];
//			 //}
//			 //workArray1.clear();
//		  // }
//
//		  // // hydrate parameters input by Bingbo.
//		  // //std::vector<double> rTem;
//		  // rTem.clear();
//		  // for (int l = 0; l < 3; l++) {
//			 //  if (myid == workingList[i]) {
//				for (int j = 0; j < mpi_size; j++) {
//					if (workingList[i] == j) continue;
//					rTem.resize(portionInconSize);
//					for (int k = 0; k < portionInconSize; k++) rTem[k] = primaryVar[k][l];
//					sendN = sendToLocalNum[j];
//					tag = workingList[i] * 10000 + j + TenM * (5 + l);
//					if (sendN > 0) {
//						mpw.SendLocalDblVariables(rTem, sendToLocalIdx, tag, sendN, j);
//					}
//				}
//				for (int k = 0; k < portionInconSize; k++) {
//					if (sendToLocalIdx[k] == myid) workArray1.push_back(primaryVar[k][l]);
//				}
//			}
//			else {
//				tag = workingList[i] * 10000 + myid + TenM * (5 + l);
//				if (localReceiveN > 0) mpw.ReceiveLocalDblVariables(workArray1, tag, localReceiveN, workingList[i]);
//			}
//
//			for (int j = 0; j < localReceiveN; j++) {
//				xk[localLocation[j]][l] = workArray1[j];
//			}
//			workArray1.clear();
//		}
//		localLocation.clear();
//
//	}
//	//
//	double tem1 = ElemMedia1[0].porosity + ElemMedia1[0].perm[0] + ElemMedia1[0].perm[1] + ElemMedia1[0].perm[2];
//	double tem;
//	pc.GetReduce(tem1, tem, "SUM");
//	if (tem == 0.0) SVParamReadFromINCON = false;
//	//			for (int i = 0; i < localGridNum; i++) std::cout <<i<<"  "<< ModelMesh.elem[i].name <<"  "<< xk[i][0] << "  " << xk[i][1] << "  " << xk[i][2] << std::endl;
//	return;
//}
double Gas_Solution_Heat(int& gas_id, double& tc, double& salt_molality, double& HK_GinH2O, double& HK_O2inH2O) {
	//   Heat of solution of non - condensible gases into water as a function of temperature and salt concentration
	double MW_Gas, t2, t3, t4, t5, tk;
	double dkhdt, dkhwdt, dkho2;
	if (gas_id < 0) MW_Gas = 2.896e1;
	else  MW_Gas = 16.0;
	t2 = tc * tc;
	t3 = t2 * tc;
	t4 = t2 * t2;
	t5 = t3 * t2;
	tk = tc + 273.15;
	if (HK_GinH2O <= 0.0) return 0.0;
	int gas_index;
	if (gas_id < 0) gas_index = -1;
	else  gas_index = 0;
	double FSalt, FSaltO, LOG_10 = log(10.0);
	switch (gas_index)
	{
	case (-1):
		//  AIR: D'Amore and Truesdell [1988], Cramer [1982], Cygan [1991]
		dkhwdt = 1.01325e5 * (1586.03 - 1.187560e01 * tc - 2.094846e-1 * t2 + 2.041320e-3 * t3
			- 6.069400e-6 * t4 + 6.002460e-9 * t5);
		dkho2 = 1.0e5 * (6.10628e2 + 14.01464 * tc - 4.178970e-1 * t2 + 28.554e-4 * t3
			- 7.7108e-6 * t4 + 7.3914e-9 * t5);

		if (salt_molality > 1.0e-12) {
			double dsogas = -2.36905e-3 + 4.84876e-05 * tc - 21.90402e-08 * t2 + 34.34892e-11 * t3;
			double dsoo2 = -1.16909e-3 + 11.1037e-06 * tc - 26.26329e-09 * t2 + 39.66268e-12 * t3;
			FSalt = LOG_10 * salt_molality * dsogas;
			FSaltO = LOG_10 * salt_molality * dsoo2;
		}
		else {
			FSalt = 0.0;
			FSaltO = 0.0;
		}
		dkhwdt = dkhwdt / HK_GinH2O + FSalt;
		dkho2 = dkho2 / HK_O2inH2O + FSaltO;
		dkhdt = 0.79 * dkhwdt + 0.21 * dkho2;
		break;
	case (0):
		//  CH4: CRAMER, 1982.
		dkhwdt = 1.0e5 * (671.091 + 13.74134 * tc - 5.192370e-01 * t2 + 4.386080e-03 * t3 - 15.97995e-06 * t4
			+ 26.77032e-09 * t5 - 16.82058e-12 * t3 * t3);
		if (salt_molality > 1.0e-12) {
			double dsogas = -1.40166e-3 + 2.6472e-05 * tc - 14.57199e-08 * t2 + 31.51868e-11 * t3
				- 27.62930e-14 * t4;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (1):
		//  C2H6: de Hemptinne, Dhima, and Shakir [2000]
		dkhwdt = 1.0e6 * (1.05e4 / tk - 2.78321e1) / tk;
		if (salt_molality > 1.0e-12) {
			double dsogas = -1.40166e-3 + 2.6472e-05 * tc - 14.57199e-08 * t2 + 31.51868e-11 * t3 - 27.62930e-14 * t4;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (2):
		// C3H8: de Hemptinne, Dhima, and Shakir [2000]
		dkhwdt = 1.0e6 * (1.265e4 / tk - 3.39873e1) / tk;
		if (salt_molality > 1.0e-12) {
			double dsogas = -1.40166e-3 + 2.6472e-05 * tc - 14.57199e-08 * t2 + 31.51868e-11 * t3 - 27.62930e-14 * t4;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (3):
		//  H2S: de Hemptinne, Dhima, and Shakir [2000]
		dkhwdt = 1.0e6 * (1.19995e04 / tk - 3.06136e1) / tk;
		dkhdt = dkhwdt / HK_GinH2O;
		break;
	case (4):
		//  CO2: CRAMER, 1982.
		dkhwdt = 1.0e5 * (19.6025 + 1.6411480 * tc - 22.22022e-3 * t2 + 8.735200e-5 * t3
			- 11.04995e-8 * t4);
		if (salt_molality > 1.0e-12) {
			double dsogas = -7.17823e-4 + 9.87708e-06 * tc - 3.11478e-08 * t2 + 4.32932e-11 * t3;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (5):
		//  N2: D'Amore and Truesdell [1988], Cygan [1991]
		dkhwdt = 1.01325e5 * (1586.03 - 11.8756 * tc - 20.94846e-2 * t2 + 20.4132e-4 * t3 - 6.0694e-6 * t4
			+ 6.00246e-9 * t5);
		if (salt_molality > 1.0e-12) {
			double dsogas = -2.36905e-3 + 4.84876e-05 * tc - 21.90402e-08 * t2 + 34.34892e-11 * t3;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (6):
		// O2: Cramer [1982], Cygan [1991]
		dkhwdt = 1.0e5 * (610.628 + 1.401464e1 * tc - 4.17897e-1 * t2 + 2.85540e-3 * t3
			- 7.71080e-6 * t4 + 7.39140e-9 * t5);

		if (salt_molality > 1.0e-12) {
			double dsogas = -1.16909e-03 + 1.11037e-05 * tc - 2.626329e-08 * t2 + 3.9662684e-11 * t3;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (7):
		//  H2O: 
		break;
	case (8):
		//  CH3CH2OH:  not ready, neglect the heat effect for this type gas disolution
		dkhdt = 0.0;
		break;
	case (9):
		//  H2: D'Amore and Truesdell [1988]
		dkhwdt = 1.01325e5 * (761.981 - 17.10334 * tc + 3.62286e-03 * t2 + 11.4208e-04 * t3
			- 7.1033e-06 * t4 + 16.875e-09 * t5 - 14.33509e-12 * t3 * t3);
		dkhdt = dkhwdt / HK_GinH2O;
		break;
	case (10):
		// n-C4H10: Carroll et al., Fluid Phase Equilibria, 140, 157-169 [1997]
		dkhwdt = 1.0e9 * (0.017687 + (13546.17 / tk - 42.85973) / tk);
		if (salt_molality > 1.0e-12) {
			double dsogas = -1.40166e-3 + 2.6472e-05 * tc - 14.57199e-08 * t2 + 31.51868e-11 * t3
				- 27.62930e-14 * t4;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	case (11):
		//  i-C4H10: de Hemptinne, Dhima, and Shakir [2000]
		dkhwdt = 1.0e6 * (2.0940e4 / tk - 5.3627e01) / tk;
		if (salt_molality > 1.0e-12) {
			double dsogas = -1.40166e-3 + 2.6472e-05 * tc - 14.57199e-08 * t2 + 31.51868e-11 * t3
				- 27.62930e-14 * t4;
			FSalt = LOG_10 * salt_molality * dsogas;
		}
		else FSalt = 0.0;
		dkhdt = dkhwdt / HK_GinH2O + FSalt;
		break;
	default:
		dkhdt = 0.0;
	}
	return -8314.56 * tk * tk * dkhdt / MW_Gas;
}
double HenryC_inv(int g_id, double& tx, double& smol) {
	static short iCH3CH2OH = 0;
	double t2, t3, t4, t5;
	double rkh, rkho2, rkhw, sogas, soo2;
	t2 = tx * tx;
	t3 = t2 * tx;
	t4 = t2 * t2;
	t5 = t3 * t2;
	double T_k = tx + 273.15;
	bool NaclImpact = true;
	if (smol <= 1.0e-10) NaclImpact = false;
	int gas_index;
	if (g_id < 0) gas_index = -1;
	else  gas_index = 0;
	// "CH4","C2H6","C3H8","H2S","CO2","N2","O2","H2O","CH3CH2OH","H2","N-C4H10", "I-C4H10", "AIR", "NH3", "C2H2","C2H4"
	switch (gas_index)
	{
	case (0):
		//  CH4: CRAMER, 1982.
		rkh = (24582.4E0 + 671.091E0 * tx + 6.87067E0 * t2 - .173079E0 * t3 + 1.09652E-03 * t4
			- 3.19599E-06 * t2 * t3 + 4.46172E-9 * t3 * t3 - 2.40294E-12 * t4 * t3) * 1.E5;
		if (NaclImpact) {
			double sogas = 1.64818E-1 - 1.40166E-3 * tx + 1.32360E-5 * t2 - 4.85733E-8 * t3
				+ 7.87967E-11 * t4 - 5.52586E-14 * t5;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;

	case(1):
		//  C2H6, de Hemptinne, Dhima, and Shakir [2000]
		rkh = 1.0e6 * exp(2.01791e2 - 1.05e04 / T_k - 2.78321e1 * log(T_k));
		if (NaclImpact) {
			double sogas = 1.64818e-01 - 1.40166e-03 * tx + 1.32360e-05 * t2
				- 4.85733e-08 * t3 + 7.87967e-11 * t4 - 5.52586e-14 * t5;
			rkh = rkh * pow(10.0, smol * (sogas + 3.2485e-2));
		}
		break;
	case(2):
		// C3H8: de Hemptinne, Dhima, and Shakir [2000]
		rkh = 1.0e6 * exp(2.44299e2 - 1.265e4 / T_k - 3.39873e1 * log(T_k));
		if (NaclImpact) {
			double sogas = 1.64818e-01 - 1.40166e-03 * tx + 1.32360e-05 * t2
				- 4.85733e-08 * t3 + 7.87967e-11 * t4 - 5.52586e-14 * t5;
			rkh = rkh * pow(10.0, smol * (sogas + 5.31775e-2));
		}
		break;
	case(3):
		// H2S: de Hemptinne, Dhima, and Shakir[2000]
		rkh = 1.0e6 * exp(2.18415e2 - 1.19995e4 / T_k - 3.06136e1 * log(T_k));
		break;
	case(4):
		//  CO2: CRAMER, 1982.
		rkh = (783.666E0 + 19.6025E0 * tx + 0.820574E0 * t2 - 7.40674E-03 * t3 + 2.18380E-05 * t4
			- 2.20999E-08 * t5) * 1.E5;
		if (NaclImpact) {
			double sogas = 1.19784E-1 - 7.17823E-4 * tx + 4.93854E-6 * t2 - 1.03826E-8 * t3 + 1.08233E-11 * t4;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;
	case(5):
		//  N2: D'AMORE AND TRUESDELL (1988),  CYGAN (1991)
		rkh = (51372.6E0 + 1586.03E0 * tx - 5.93780E0 * t2 - 6.98282E-2 * t3 + 5.10330E-4 * t4
			- 1.21388E-6 * t2 * t3 + 1.00041E-9 * t3 * t3) * 1.01325E5;
		if (NaclImpact) {
			double sogas = .183369E0 - 2.36905E-3 * tx + 2.42438E-5 * t2 - 7.30134E-8 * t3 + 8.58723E-11 * t4;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;
	case(6):
		//	O2: Cramer [1982], Cygan [1991]
		rkh = 1.0e5 * (26234.0 + 610.628 * tx + 7.00732 * t2 - 1.39299e-1 * t3 + 7.13850e-4 * t4
			- 1.54216e-6 * t5 + 1.23190e-9 * t3 * t3);
		if (NaclImpact) {
			double sogas = 0.16218 - 1.16909e-03 * tx + 5.55185e-6 * t2 - 8.75443e-09 * t3 + 9.91567e-12 * t4;
			rkh = rkh * pow(10.0, smol * sogas);
		}

		break;
	case(7):
		// H2O
		rkh = 1.0;
		break;
	case(8):
		// CH3CH2OH
		rkh = (24582.4E0 + 671.091E0 * tx + 6.87067E0 * t2 - .173079E0 * t3 + 1.09652E-03 * t4
			- 3.19599E-06 * t2 * t3 + 4.46172E-9 * t3 * t3 - 2.40294E-12 * t4 * t3) * 1.E5;
		if (NaclImpact) {
			double sogas = 1.64818E-1 - 1.40166E-3 * tx + 1.32360E-5 * t2 - 4.85733E-8 * t3
				+ 7.87967E-11 * t4 - 5.52586E-14 * t5;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;
	case(9):
		//  H2: D'AMORE AND TRUESDELL, 1988
		rkh = (57106.3E0 + 761.981E0 * tx - 8.55167E0 * t2 + 1.20762E-3 * t3 + 2.85520E-04 * t4
			- 1.42066E-06 * t2 * t3 + 2.81250E-9 * t3 * t3 - 2.04787E-12 * t4 * t3) * 1.01325E5;
		break;
	case(10):
		// n - C4H10 : Carroll et al., Fluid Phase Equilibria, 140, 157 - 169[1997]
		rkh = 1.0e9 * exp(285.8038 + 0.017687 * T_k - 13546.17 / T_k - 42.85973 * log(T_k));
		if (NaclImpact) {
			double sogas = 1.64818e-01 - 1.40166e-03 * tx + 1.32360e-05 * t2 - 4.85733e-08 * t3 + 7.87967e-11 * t4 - 5.52586e-14 * t5;
			rkh = rkh * pow(10.0, smol * (sogas + 7.16717e-2));       //Soreide and Wilson
		}
		break;
	case(11):
		// i-C4H10: Tsonopoulos and Wilson, Fluid Phase Equilib. 29, 391Ð414 [1983]
		rkh = 1.0e6 * exp(3.8464e2 - 2.0940e4 / T_k - 5.3627e01 * log(T_k));
		if (NaclImpact) {
			double sogas = 1.64818e-01 - 1.40166e-03 * tx + 1.32360e-05 * t2 - 4.85733e-08 * t3 + 7.87967e-11 * t4 - 5.52586e-14 * t5;
			rkh = rkh * pow(10.0, smol * (sogas + 6.666343e-2));       //Soreide and Wilson
		}
		break;
	case(12):
		//  AIR: D'AMORE AND TRUESDELL (1988), CRAMER (1982), CYGAN (1991)
		rkhw = (51372.6E0 + 1586.03E0 * tx - 5.93780E0 * t2 - 6.98282E-2 * t3
			+ 5.10330E-4 * t4 - 1.21388E-6 * t2 * t3 + 1.00041E-9 * t3 * t3) * 1.01325e5;
		rkho2 = (26234.0E0 + 610.628E0 * tx + 7.00732E0 * t2 - .139299E0 * t3
			+ 7.13850E-4 * t4 - 1.54216E-6 * t2 * t3 + 1.23190E-9 * t3 * t3) * 1.01325e5;
		if (NaclImpact) {
			sogas = .183369E0 - 2.36905E-3 * tx + 2.42438E-5 * t2 - 7.30134E-8 * t3
				+ 8.58723E-11 * t4;
			soo2 = 0.16218E0 - 1.16909E-3 * tx + 5.55185E-6 * t2 - 8.75443E-9 * t3 + 9.91567E-12 * t4;
			rkh = 0.79 * rkhw * pow(10.0, smol * sogas) + 0.21 * rkho2 * pow(10.0, smol * soo2);
		}
		else {
			rkh = 0.79 * rkhw + 0.21 * rkho2;
		}
		break;
	case(13):
		// NH3 : D'AMORE AND TRUESDELL (1988), only up to 320°C 
		rkh = 2.03924 + 9.51358e-02 * tx - 2.85096e-03 * t2 + 3.69927e-05 * t3 - 5.93259e-08 * t4;
		break;
	case(14):
		//  C2H2, not ready using CH4 data
		rkh = (24582.4E0 + 671.091E0 * tx + 6.87067E0 * t2 - .173079E0 * t3 + 1.09652E-03 * t4
			- 3.19599E-06 * t2 * t3 + 4.46172E-9 * t3 * t3 - 2.40294E-12 * t4 * t3) * 1.E5;
		if (NaclImpact) {
			sogas = 1.64818E-1 - 1.40166E-3 * tx + 1.32360E-5 * t2 - 4.85733E-8 * t3
				+ 7.87967E-11 * t4 - 5.52586E-14 * t5;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;
	case(15):
		//  C2H4, not ready using CH4 data
		rkh = (24582.4E0 + 671.091E0 * tx + 6.87067E0 * t2 - .173079E0 * t3 + 1.09652E-03 * t4
			- 3.19599E-06 * t2 * t3 + 4.46172E-9 * t3 * t3 - 2.40294E-12 * t4 * t3) * 1.E5;
		if (NaclImpact) {
			sogas = 1.64818E-1 - 1.40166E-3 * tx + 1.32360E-5 * t2 - 4.85733E-8 * t3
				+ 7.87967E-11 * t4 - 5.52586E-14 * t5;
			rkh = rkh * pow(10.0, smol * sogas);
		}
		break;
	}
	return rkh;
}
//

// main program
int main() {
	double tx = 1.2, a;
	double Hc_O2inH2O = 1.0;
	int gid = 0;
	double xs = 0.0, Hc_GinH2O = 3e-10;
	double HK_GinH2O = 1.0 / Hc_GinH2O;
	//double hc = HenryC_inv(gid, tx, xs);
	//double EnthalpyCH4Solution = realG.Gas_Solution_Heat(gid, tx, Xmol_iA, Hc_GinH2O, Hc_O2inH2O);
	double EnthalpyCH4Solution = Gas_Solution_Heat(gid, tx, xs, HK_GinH2O, Hc_O2inH2O);
	a = fabs(-tx);
	cout << "a = " << a << endl;

	//int NumElemTot = 3, NumCom = 3, NumComPlus1 = NumCom + 1, Nacl=4, Heat = 3;
	////std::vector<std::vector<double>> xk;
	//double xk[3][4];
	//xk[0][0] = 2.0e6;
	//xk[0][1] = 2.002e2;
	//xk[0][2] = 2.0e-1;
	//xk[0][3] = 2.0e1;
	//xk[0][4] = 2.0e-1;
	//xk[1][0] = 5.0e6;
	//xk[1][1] = 5.0e2;
	//xk[1][2] = 5.0e-1;
	//xk[1][3] = 5.0e1;
	//xk[1][3] = 5.0e-1;
	//xk[2][0] = 8.0e6;
	//xk[2][1] = 8.0e2;
	//xk[2][2] = 8.0e-1;
	//xk[2][3] = 8.0e1;
	//xk[2][3] = 8.0e-1;
	//struct State_Conditions {
	//	unsigned char index0, index1;		// Element previous and current phase state index
	//	unsigned char changes;				// Number of state changes per Newtonian Iteration
	//	double pres0, pres;					// Pressure
	//	double temp0, temp;					// Temperature
	//	double param1, param2;				// Auxiliary parameters stored at state point
	//	double iniPorosity;					//system initial porosity
	//	double alphaS;						// sum of saturation for all phases, for use in Relaxed volume method
	//	double bulbP0, bulbP;
	//};
	//std::vector<State_Conditions> ElemState;
	//int pm1 = 0, wi = 0;
	//short ipos = 2;
	//bool noAqu2 = true, noNacl = true, CompositionModel = false, nonIsothermal = true;
	//double tem = 0.0, intpart;
	//if (noAqu2) ipos = 1;
	//for (int i = 0; i < 1; i++) {
	//	if (xk[i][0] <= 1000.0 || xk[i][0] > 2.099e+08) pm1++;		// for H2O system.
	//	if (xk[i][1] >= 100.0 && xk[i][1] <= 701.0) {
	//		double tem1 = modf(xk[i][1], &intpart);
	//		short idx = int(intpart / 100.0 + 0.01);
	//		//ElemState[i].index1 = idx;
	//		xk[i][1] = tem1;
	//	}
	//	if (!noAqu2 && CompositionModel) { if (xk[i][1] < 0.0 || xk[i][1]>1.0) pm1++; }
	//	bool tb1 = (xk[i][ipos] >= 0.0 && xk[i][ipos] <= 1.0);
	//	bool tb2 = (xk[i][ipos] >= 10.0 && xk[i][ipos] <= 11.0);
	//	if (!tb1 && !tb2) pm1++;
	//	for (int j = ipos + 1; j < NumCom; j++) {
	//		if (noNacl || (!noNacl && j != Nacl)) tb1 = (xk[i][j] >= 0.0 && xk[i][j] <= 1.0);
	//		if (!tb1) pm1++;
	//	}
	//	if (nonIsothermal) {
	//		if (xk[i][Heat] <= -50.0 || xk[i][Heat] > 800) pm1++;
	//	}
	//	if (!noNacl) {
	//		tb1 = (xk[i][Nacl] >= 0.0 && xk[i][Nacl] <= 1.0);
	//		tb2 = (xk[i][Nacl] >= 10.0 && xk[i][Nacl] <= 11.0);
	//		if (!tb1 && !tb2) pm1++;
	//	}

	//	tem = tem + xk[i][Heat];

	//	if (pm1 > 0) {
	//		LogFile << "Initialization wrong at gridblock '" << i;
	//		LogFile << "' with";
	//		for (int j = 0; j < NumComPlus1; j++) LogFile << "  xk_" << j + 1 << "=" << xk[i][j];
	//		LogFile << std::endl;
	//		wi++;
	//		pm1 = 0;
	//	}
	//}

	//if (tem == 0.0) {
	//	LogFile << "Temperature of the flow system has not been initialized (all are 0.0), stop simulation!" << std::endl;
	//	wi = 1;
	//}
	//if (wi > 0) exit(1);

	//struct Fluid {
	//	std::string FluidName;
	//	double Density;              // Fluid density(kg / m3)
	//	double Beta;                 // Fluid compressibility factor,  beta. rho_air= beta*p_air (kg/m**3/Pa)
	//	double Visco;                // Fluid viscosity (Pa.s)
	//	double ThermalExansionC;     // Thermal expansion coefficient(1/ C^o)
	//	double BinghamG;             // Bingham parameter G   (Pa/m)
	//	double PowerLawN;            // Power-law n;
	//	double PowerLawH;            // Power-law H;
	//	double MinViscoPLFluid;      // Min.viscosity for power - law fluid, default = 1.e-5 (Pa.s)
	//	double ThermalCondc;         // Fluid thermal conductivity (W / m / C)
	//	short FluidType;
	//	double ref_P, ref_T;
	//	double parameters[10];
	//};
	//std::vector<Fluid> Fluids;
	//// Write equation of state informations
	//bool noGasSim = false;
	//int ien = 0;
	//if (noGasSim) ien = 1;
	//int ig = ien;
	//int NumPhases = 2;
	//for (int i = ig; i < NumPhases - 1; i++) LogFile << Fluids[i].FluidName << ", ";
	//if (!noNacl) LogFile << " Nacl, ";
	//LogFile << "and " << Fluids[NumPhases - 1].FluidName << " Fluid Mixtures";
	//if (nonIsothermal) LogFile << " Under Non-isothermal Condition." << std::endl;
	//else LogFile << " Under Isothermal Condition." << std::endl;

	//int NumGases = 1;
	//bool noGas = false;
	////if (!noGasSim) {
	////	LogFile << "Number of Gas Component: " << NumGases << std::endl;
	////	LogFile << "Gases includes: ";
	////	for (int i = 0; i < NumGases - 1; i++) LogFile << ComponentName[GetGasIdx(i)] << ", ";
	////	LogFile << ComponentName[GetGasIdx(NumGases - 1)];
	////	if (General_Control_Parameters::realGasAccountingVapor) LogFile << ", and Vapor.";
	////	LogFile << std::endl;
	////}
	////if (AccountingForDiffusion)  LogFile << "Diffusion is accounted in the simulation." << std::endl;
	////else LogFile << "No diffusion is accounted in the simulation." << std::endl;

	////LogFile << "Number of Fluid Components: " << NumCom << " (";
	////for (int i = 0; i < NumComPlus1 - 1; i++) LogFile << ComponentName[i] << ", ";
	////LogFile << "and " << ComponentName[NumComPlus1 - 1];
	////LogFile << ")" << std::endl;
	////LogFile << "Number of Equations per Gridblocks: " << NumEqu << std::endl;
	////LogFile << "Number of Phase Can be Present: " << NumAllPhases << " (";
	////LogFile << "Gas, ";
	////LogFile << "Aqueous";
	////if (Fluids[Gas].FluidName == "CO2") LogFile << ", Liquid CO2";
	////if (HasSolidPhase) LogFile << ", Solid";
	////if (!noAqu2) LogFile << " and Aqueous2";
	////LogFile << ")" << std::endl;


	//std::string GasFluidName = "HYDRATE";


	//if (Fluids[0].FluidName != GasFluidName && Fluids[0].FluidName != "MIXGASES") {
	//	for (int i = 0; i < NumElemTot; i++) {
	//		if (noGas) ElemState[i].index1 = 1;
	//		else {
	//			if (xk[i][ipos] >= 10.0) ElemState[i].index1 = 2;
	//			else ElemState[i].index1 = 1;
	//		}
	//	}
	//}
	//// f9 conditional statements
	//int a, b, c;
	//cout << "Input value of a: \n";
	//cin >> a;
	//cout << "Input value of b: \n";
	//cin >> b;
	////
	//c = (a - b) < 0 ? (a + b) : (b - a);
	//std::vector<std::string> stateIdx;
	//stateIdx.push_back("0");
	//stateIdx.push_back("1");
	//if (stateIdx[1] == 1) {
	//	cout << "ok" << endl;
	//}




	//// Test int/bool
	//int a = 2;
	//while (a) {
	//	cout << "int == bool --> a != 0 means true; only if a == 0, it means false" << endl;
	//	a--;
	//}

	//



	//int N_an;
	//cin >> N_an;
	//std::vector<double> an;
	//for (int i = 0; i < N_an; i++) {
	//	an.push_back(i + 4);
	//}
	//double b = an[1];
	//cout << "b = " << b << endl;
	//std::vector <double> temVar;
	//temVar.push_back(1);
	//temVar.push_back(2);
	//temVar.push_back(3);
	//temVar.resize(7);


	//
	LogFile.open("test.log");
	MainInputFile.open("test.dat");
	/*std::string fileName;
	FILE* ReadInputFile = NULL;
	fileName = "test.dat";
	ReadInputFile = fopen(fileName.c_str(), "r");

	char szLine[128];
	double dbls[10];
	std::string strs[8];
	std::vector < std::vector <double>> primaryVar;
	std::vector < std::vector <double>> hydrateVar;

	int ints[8];
	int NumCom = 6;
	ReadGetLine(szLine, 128, ReadInputFile);
	ReadFormattedTextLine(szLine, "D10D10D10D10D10D10D10D10D10D10", strs, dbls, ints);
	for (int i = 0; i < NumCom + 1; i++) temVar[i] = dbls[i];
	primaryVar.push_back(temVar);
	temVar.clear();
	temVar.resize(3);
	for (int i = NumCom + 1; i < NumCom + 4; i++) temVar[i - NumCom - 1] = dbls[i];
	hydrateVar.push_back(temVar);*/
	//
	//readMESH_Conne();
	//readMESH_Eleme();
	//readSOLVE();
	//readINCON();
	//readTIMES();
	readFLUID();
	//readPVTDA();
	//readMODEL();
	//readRPCAP();
	//readTIMBC();
	//readINDOM();
	//readDIFFU();

	//readHYDRA();

	//readHYDRATE();

	//vector<string>  KeyWordList;
	//string temStrs[] = {
	//		"TITLE","ROCKS","INCON","PARAM","GENER",
	//		"ELEME","CONNE","MULTI","START","TIMES",
	//		"FLUID","PVTDA","SOLVE","ENDCY","MODEL",
	//		"FOFTS","COFTS","GOFTS","ENDFI","RPCAP",
	//		"TIMBC","INDOM","DIFFU","INGEQ","SELEC",
	//		"ADVAN","SVPAR","TRCON","OUTCT","PREME",
	//		"INICO","WELLS","GRIDB","PVT_D","HYDRA" };
	//int a = sizeof(temStrs);
	//for (int i = 0; i < sizeof(temStrs) / sizeof(temStrs[0]); i++) KeyWordList.push_back(temStrs[i]);

	//int b = sizeof(temStrs) / sizeof(temStrs[0]);



	//int i = 2, j = 4;
	//while (true) {
	//	i++;		// the same as i += 1;
	//	cout << "ok4" << endl;
	//	if (i < j) continue;
	//	if (i > j) break;
	//}
	//cout << "ok5" << endl;


	//string Title_tmp;
	//MainInputFile.open("test.log");
	//cout << MainInputFile.is_open() << endl;
	//getline(MainInputFile, Title_tmp);
	//cout << Title_tmp << endl;
	//if (Title_tmp.length() <= 0) {
	//	cout << "ok1" << endl;
	//}
	//else if (Title_tmp.length() > 200) {
	//	//strcpy(Project_Title, Title_tmp.substr(0, 200).c_str());
	//	cout << "ok2" << endl;
	//}
	//else {
	//	//strcpy(Project_Title, Title_tmp.c_str());
	//	cout << "ok3" << endl;
	//}

	//MainInputFile.open("test.log");
	//cout << MainInputFile.is_open() << endl;
	//cout << !MainInputFile << endl;
	//MainInputFile << "ok" << endl;	// set:fstream MainInputFile means 
	//while (!MainInputFile) {
	//	cout << "ok" << endl;
	//}

	return 0;
}