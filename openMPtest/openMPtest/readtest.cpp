#include <iostream>;
#include <fstream>;
#include <string>;		// for getline function
#include <vector>;
#include <algorithm>;
#include <iomanip>; // for setw and precision


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
	cout << LineText.length() << endl;
	cout << LineText.substr(0, 2) << endl;
	if (LineText.length() < 2 || LineText.substr(0, 2).compare("//") == 0) return 1;

	std::string strTmp = LineText;
	transform(strTmp.begin(), strTmp.end(), strTmp.begin(), ::toupper);

	for (int i = 0; i < AryDescribles.size(); i++)
	{
		cout << (strTmp.find(AryDescribles[i]) == -1) << endl;

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
			//timeUnit = stringToNum<int>(strSub);
		}
		break;
		}
	}
}


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
	int NCom;
	int hydrN;
	int N_ThC;
	int N_SpH;
	int N_Rho;
	int EquationOption;

	double moleF;
	std::vector<double> An, Bn, Dn;
	std::string nameG;
	InhibitorParameters inhibp;
}
// ReadHYDRA
void readHYDRA() {
	using namespace Hydrate_Related_Parameters;

	//if (myid == iMaster) {
	//	HYDRA_read = true;
	LogFile << std::endl << "--->>> Start reading HYDRA" << std::endl;
	//}

	// Initialization of all hydrate related parameters
	//int NCom, hydrN, N_ThC, N_SpH, N_Rho;
	//double moleF, Max_TShift, Y_atMax_TShift, InhibitorMW, InhibitorDens, InhibitorEnthSol, InhibitorCpCoeff, EquationOption;
	//int inhibitor_flag;
	//std::vector<double> An, Bn, Dn;
	//std::string nameG;


	std::vector<std::string> AryDescribles;
	AryDescribles.push_back("NCom");
	//????? When Ncom > 1, need i cycle? 
	AryDescribles.push_back("nameG");
	AryDescribles.push_back("hydrN");
	AryDescribles.push_back("moleF");
	//AryDescribles.push_back("MolWt");
	//AryDescribles.push_back("GasMF");
	//AryDescribles.push_back("H2OMF");

	AryDescribles.push_back("N_ThC");
	AryDescribles.push_back("An");
	AryDescribles.push_back("N_SpH");
	AryDescribles.push_back("Bn");
	AryDescribles.push_back("N_Rho");
	AryDescribles.push_back("Dn");

	AryDescribles.push_back("inhibitor_flag");
	AryDescribles.push_back("Max_TShift");
	AryDescribles.push_back("Y_atMax_TShift");
	AryDescribles.push_back("InhibitorMW");
	AryDescribles.push_back("InhibitorDens");
	AryDescribles.push_back("InhibitorEnthSol");
	AryDescribles.push_back("InhibitorCpCoeff");
	AryDescribles.push_back("EquationOption");
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
				NCom = atoi(strSub.c_str());
				////range check!
				//Max_NumNRIterations = rangeCheck(iTem, "Maximum number of Newtonian iterations", 1, 100, 8);
			}
		}
		break;
		case 1:
		{
			if (strSub.length() > 0) {
				nameG = strSub.c_str();
				//if (Max_NumTimeSteps < 1) Max_NumTimeSteps = 1;
			}
		}
		break;

		case 2:
		{
			if (strSub.length() > 0) hydrN = atoi(strSub.c_str());
		}
		break;
		case 3:
		{
			if (strSub.length() > 0) {
				moleF = atof(strSub.c_str());
				//if (CPU_MaxTime < 0.1) CPU_MaxTime = 0.1;
				//CPU_MaxTime = CPU_MaxTime * 60;    //input in minutes
			}
		}
		break;
		case 4:
		{
			if (strSub.length() > 0) N_ThC = atoi(strSub.c_str());
			//if (OutputControl::output_frequency < 1) OutputControl::output_frequency = 100000;
		}
		break;
		case 5:
		{
			//if (strSub.length() > 0) An = atof(strSub.c_str());
		}
		break;
		case 6:
		{
			if (strSub.length() > 0) N_SpH = atoi(strSub.c_str());
		}
		break;
		case 7:
		{
			//if (strSub.length() > 0) Bn = atoi(strSub.c_str());
		}
		break;
		case 8:
		{
			if (strSub.length() > 0) N_Rho = atoi(strSub.c_str());
		}
		break;
		case 9:
		{
			//if (strSub.length() > 0) Dn = atoi(strSub.c_str());
		}
		break;
		case 10:
		{
			if (strSub.length() > 0) inhibp.inhibitor_flag = atoi(strSub.c_str());
		}
		break;
		case 11:
		{
			if (strSub.length() > 0) inhibp.Max_TShift = atof(strSub.c_str());
		}
		break;
		case 12:
		{
			if (strSub.length() > 0) inhibp.Y_atMax_TShift = atof(strSub.c_str());
		}
		break;
		case 13:
		{
			if (strSub.length() > 0) inhibp.InhibMW = atof(strSub.c_str());
		}
		break;
		case 14:
		{
			if (strSub.length() > 0) inhibp.InhibDens = atof(strSub.c_str());
		}
		break;
		case 15:
		{
			if (strSub.length() > 0) inhibp.InhibEnthSol = atof(strSub.c_str());
		}
		break;
		case 16:
		{
			if (strSub.length() > 0) inhibp.InhibCpCoeff = atof(strSub.c_str());
		}
		break;
		case 17:
		{
			if (strSub.length() > 0) EquationOption = atoi(strSub.c_str());
		}
		break;
		default:
			break;
		}
	}
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


// Test of INCON reading
void readINCON()
{
	//using namespace General_External_File_Units;
	//using namespace Basic_Parameters;
	std::string str_tmp = " ";
	short nPrimary = 4;
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
			if (nTmp >= (nPrimary + 4)) break;
		}

		strOut = strName;
		if (strName.length() == 0) { getline(MainInputFile, LineText); Line_Num += 1; continue; }
		for (int i = strName.length(); i < 15; i++) strOut += " ";
		incon_out << strOut;

		for (int i = nPrimary; i < nPrimary + 4; i++)
		{
			if (xn[i] != -999.666) incon_out << std::right << std::setw(15) << std::setprecision(7) << std::scientific << xn[i];
			//cout << xn[i] << endl;
		}
		incon_out << std::endl;

		for (int i = 0; i < nPrimary; i++)
		{
			if (xn[i] != -999.666) incon_out << std::right << std::setw(20) << std::setprecision(10) << std::scientific << xn[i];
		}
		incon_out << std::endl;

		getline(MainInputFile, LineText);
		Line_Num += 1;
	}
	incon_out.close();

	LogFile << std::endl << "--->>> End reading INCON" << std::endl;
	LogFile << std::endl << "===============================================================" << std::endl << std::endl;
}

// main program
int main() {

	//
	LogFile.open("test.log");
	MainInputFile.open("test.dat");
	//
	//readSOLVE();
	readINCON();
	readTIMES();
	readHYDRA();

	vector<string>  KeyWordList;
	string temStrs[] = {
			"TITLE","ROCKS","INCON","PARAM","GENER",
			"ELEME","CONNE","MULTI","START","TIMES",
			"FLUID","PVTDA","SOLVE","ENDCY","MODEL",
			"FOFTS","COFTS","GOFTS","ENDFI","RPCAP",
			"TIMBC","INDOM","DIFFU","INGEQ","SELEC",
			"ADVAN","SVPAR","TRCON","OUTCT","PREME",
			"INICO","WELLS","GRIDB","PVT_D","HYDRA" };
	int a = sizeof(temStrs);
	for (int i = 0; i < sizeof(temStrs) / sizeof(temStrs[0]); i++) KeyWordList.push_back(temStrs[i]);

	int b = sizeof(temStrs) / sizeof(temStrs[0]);



	int i = 2, j = 4;
	while (true) {
		i++;		// the same as i += 1;
		cout << "ok4" << endl;
		if (i < j) continue;
		if (i > j) break;
	}
	cout << "ok5" << endl;


	string Title_tmp;
	MainInputFile.open("test.log");
	cout << MainInputFile.is_open() << endl;
	getline(MainInputFile, Title_tmp);
	cout << Title_tmp << endl;
	if (Title_tmp.length() <= 0) {
		cout << "ok1" << endl;
	}
	else if (Title_tmp.length() > 200) {
		//strcpy(Project_Title, Title_tmp.substr(0, 200).c_str());
		cout << "ok2" << endl;
	}
	else {
		//strcpy(Project_Title, Title_tmp.c_str());
		cout << "ok3" << endl;
	}





	//MainInputFile.open("test.log");
	//cout << MainInputFile.is_open() << endl;
	//cout << !MainInputFile << endl;
	//MainInputFile << "ok" << endl;	// set:fstream MainInputFile means 
	//while (!MainInputFile) {
	//	cout << "ok" << endl;
	//}

	return 0;
}