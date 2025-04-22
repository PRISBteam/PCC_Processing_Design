///================================ DCC Edges types statistics, indices and configuration entropy =============================================================///
///==============================================================================================================================///
/** This subfile calculates the statistical characteristics of Edges incuding their Face indices and configurational entropy  **/
///==============================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random> // Require C++ 11 and above
// #include <execution> // Require C++ 17 and above

// external libraries
#include <Eigen/Core>
#include <Eigen/SparseCore>

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'


vector<double> EdgesStat( std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> const& CellNumbs, Eigen::SparseMatrix<double> const& FES, char* output_dir, double &Face_Entropy_Median, double &Face_Entropy_Skrew, double &informativeness, vector<double> &j_types_fractions)
{
    std::vector<double> TJsTypes(CellNumbs.at(1),0);
    map<unsigned int, unsigned int> res; // Here 100 is an arbitrary number of Edge types
    map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    Face_Entropy_Median = 0; Face_Entropy_Skrew = 0;

    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
//    for (auto sfe : TJsTypes) cout << sfe << "\t" ; cout << endl;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
    double j0s = j0, j1s = j1, j2s = j2, j3s = j3;
    /// using values with pow(10,-10) instead of 0s!
//    if (j0s == 0) j0s = pow(10,-30); if (j1s == 0) j1s = pow(10,-30); if (j2s == 0) j2s = pow(10,-30); if (j3s == 0) j3s = pow(10,-30); //Gives 0 in entropy!
    if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!
    /// Configuration Entropy related with Faces
//    Configurational_Face_Entropy = - (j0s* log2(j0s) + j1s* log2(j1s) + j2s* log2(j2s) + j3s* log2(j3s));
    Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

    /// Median part in the entropy decomposition
    j_types_fractions = {j0, j1, j2, j3}; /// using values with pow(10,-10) instead of 0s!
    if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
        Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(j0 * j1 * j2 * j3);
    } else Face_Entropy_Median = 0.0;

    /// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
    for (int j = 0; j < j_types_fractions.size(); j++)
        for (int i = 0; i < j; i++)
            if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
                Face_Entropy_Skrew +=
                        -(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
                        log2(j_types_fractions[i] / j_types_fractions[j]);
            } else Face_Entropy_Skrew += 0.0;
/*
    /// Complexity/ Informativeness
    // Srand and Smax reading from files
    string input_filename_SstrRand = "Random_Entropy_100.txt"s, input_filename_SstrMAX = "Maximum_Entropy_100.txt"s;
    string input_RandomEntropy_dir = output_dir + input_filename_SstrRand, input_MAXEntropy_dir = output_dir + input_filename_SstrMAX;
    char* RandomEntropy_dir = const_cast<char*>(input_RandomEntropy_dir.c_str()); char* MAXEntropy_dir = const_cast<char*>(input_MAXEntropy_dir.c_str()); // From string to char for the passing folder path to a function
    // Format of the elements:: special Face fraction _ Conf entropy value _ Mean entropy value (log(p1*p2*..*pn))
    vector<vector<double>>  RandomEntropy = VectorVectors4Reader(RandomEntropy_dir);
    vector<vector<double>>  MAXEntropy = VectorVectors4Reader(MAXEntropy_dir);
    //for ( auto tpl : MAXEntropy) cout << tpl[0] << " " << tpl[1] << " " << tpl[2] << endl;

    /// Index of complexity:
    double RandomEntropy_p = 0, MAXEntropy_p = 0;
    vector<double> Delta_MAX;
    double sff =(double) s_faces_sequence.size()/ CellNumbs.at(2); // special face fraction
    for ( auto tpl : MAXEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Smax = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    for ( auto tpl : RandomEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Srand = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    /// Informativeness parameter
//REPAIR    cout << numb_Smax << " " << numb_Srand << endl;
//REPAIR    cout << Configurational_Face_Entropy << " " << RandomEntropy[numb_Srand][1] << " " << MAXEntropy[numb_Smax][1] << endl;
    informativeness =  ( Configurational_Face_Entropy - RandomEntropy[numb_Srand][1])/ ( MAXEntropy[numb_Smax][1] - RandomEntropy[numb_Srand][1]);
    if (informativeness > 1) informativeness = 1;
*/
    /// ====== Data output =====================>
    /// Opening of the output streams
    string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    double special_faces_fraction =(double) s_faces_sequence.size()/ CellNumbs.at(2);
//    cout << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t"<< Face_Entropy_Median << endl;

    OutTJsFile << special_faces_fraction << "\t" << j0 << "\t" << j1 << "\t" << j2 << "\t" << j3 << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t\t" << Face_Entropy_Skrew << endl;
    OutTJsFile.close();

//    int b = 0;
//    for (unsigned int it : TJsTypes) res[b++] = it;

    return TJsTypes;
}
/**
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////// Обработка матриц соседей: начльная статистика связей и матрицы смежности /////////////////////////////////////////////////////
	long lt=0; knn=1;
// Для каждой грани мы строим матрицу JEdgeNeigh[i][j] где на первом месте ее собственный J-тип а дальше список J-типов ее соседей
long knumb=0, J0Ncount=0, J00Ncount=0; 	knn=0;
	int **JEdgeNeigh;
	JEdgeNeigh = new int* [Edgenumb];
	for (int i = 0; i < Edgenumb; i++) JEdgeNeigh[i] = new int [100];	//снова закладываем до 100 гипотетических соседей
	int **JFaceNeigh;
	JFaceNeigh = new int* [Facenumb];
	for (int i = 0; i < Facenumb; i++) JFaceNeigh[i] = new int [100];	//снова закладываем до 100 гипотетических соседей

		for (int i = 0; i < Edgenumb; i++)
			for (int j = 0; j < 100; j++)
				JEdgeNeigh[i][j] = -6;
		for (int i = 0; i < Facenumb; i++)
			for (int j = 0; j < 100; j++)
				JFaceNeigh[i][j] = -1; //Все LAGBs!!!

//Сначала мы выясняем тип самой грани
for(int i = 0; i < Edgenumb; i++)	{
	for (int lk = 0; lk < Facenumb; lk++) if(MFE1[lk][i] == 1) J00Ncount++; //then ==2
		JEdgeNeigh[i][0] = J00Ncount;
			J00Ncount=0;
	lt=0; knn=1;


//Затем проходим всех ее соседей
	do{
		if(EdgeNeighbours[i][lt]>=0) knumb = EdgeNeighbours[i][lt]; else knumb=-1;
//------------------------------------------------------

		if(knumb>=0) for (int lk = 0; lk < Facenumb; lk++) if(MFE1[lk][knumb] == 1) J0Ncount++;  //then ==2

//-----------------------------------------------------		//cout <<i<< "   "<< knn << "   "<< JEdgeNeigh[0][1] << endl;
		JEdgeNeigh[i][knn++] = J0Ncount;

		J0Ncount=0;	lt++;
	}while(knumb>=0);
}

//Вывод в файл JEdgeNeigh.txt
			JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEdgeNeigh.txt", ios::trunc);
			for (int i = 0; i < Edgenumb; i++) {
				for (int j = 0; j < 100; j++) {
					if(JEdgeNeigh[i][j]>=0) JEdgeN << JEdgeNeigh[i][j] << "\t";
				}
			JEdgeN <<"Edge number="<<i<<"\n";
		}
        JEdgeN.close();

//cout<<"initial (HAGBs)JEdgeNeigh.txt has been created"<<endl;

//Чистка файла под мощности тройных стыков с учетом соседей
	JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)TJspow.txt", ios::trunc);
	JEdgeN.close();


	int **JEN;
	JEN = new int* [30];
	for (int i = 0; i < 30; i++) JEN[i] = new int [30];
			for (int i = 0; i < 5; i++)
				for (int j = 0; j < 5; j++)
					JEN[i][j] = -10;

//Анализ матрицы JEdgeNeigh[i][j]
for (int i = 0; i < Edgenumb; i++) //для каждой грани
	for (int l = 0; l < 5; l++) //мы перебираем все варианты какой она может быть (с запасом до 30, хотя реально до 4-5 видов)
		if(JEdgeNeigh[i][0] == l) for (int j = 1; j < 100; j++) //и если она оказалась определенного типа, то мы перебираем всех ее соседей
										for (int k = 0; k < 5; k++)
											if(JEdgeNeigh[i][j] == k) JEN[l][k]++; //так что если сосед оказывается также определенного типа, то мы заносим их связь в матрицу JEN
//очевидно, что при таком алгоритме диагональные элементы учитываются дважды, то есть
for (int l = 0; l < 5; l++)
	for (int k = 0; k < 5; k++)
		if(k == l)  JEN[l][k]= 0.5*JEN[l][k];

//Вывод в файл начального распределения узлов по количеству соседей разного типа JEN.txt
//(сколько 1-узлов соединено с 1-узлами, 2 - узлами..., сколько 3-узлов соединено с 0 - узлами, 1 - узлами и тд)
			JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEN00.txt", ios::trunc);
			for (int i = 0; i < 5; i++) {
				for (int j = 0; j < 5; j++) {
					if(JEN[i][j]>=0) JENStream << JEN[i][j] << "\t";
//					if(JEN[i][j]>=0) JENStream << JEN[i][j]/Edgenumb << "\t"; //В долях /Edgenumb
				}
				JENStream <<"\n";
		}
        JENStream.close();

		JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Jlk.txt", ios::trunc);  //			 JENStream << JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + JEN[0][2] + JEN[2][0] + JEN[0][3] + JEN[3][0] + JEN[1][2] + JEN[2][1] + JEN[1][3] + JEN[3][1] + JEN[2][3] + JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< JEN[0][2] + JEN[2][0]<< "\t"<< JEN[0][3] + JEN[0][3]<< "\t"<< JEN[1][2] + JEN[2][1]<< "\t"<< JEN[1][3] + JEN[3][1]<< "\t"<< JEN[2][3] + JEN[3][2];
		JENStream.close();

//*-------------------------------------------------------------------------------------
//*------------------------------QUADRUPLE POINTS----------------------------
//*-------------------------------------------------------------------------------------
//*-------------------------------------------------------------------------------------
//Ji-Jj analyser
//--------------------------------------------------------------------------------------
for (int in = 0; in < Edgenumb; in++)
for (int j = 0; j < 100; j++)
JEdgeNeigh[in][j] = -10;

for(int ih = 0; ih < Edgenumb; ih++)	{
knn=0;
//Сначала мы выясняем тип самой грани
J1count=0; 	JN0Count=0;
for (int lk = 0; lk < Facenumb; lk++) { if(MFE1[lk][ih] == 2) JN0Count++; if(MFE1[lk][ih] >= 1) J1count++; }
if(J1count > 1) JEdgeNeigh[ih][0] = JN0Count;
//			{cout<<"  il= "<<ih<<"  JEdgeNeigh0= "<<JEdgeNeigh[ih][0]<<endl; system("pause");}
//Сначала в 0-ячейку собственный тип ребра            //
}

for(int ih = 0; ih < Edgenumb; ih++)	{
knn=1;
//Затем проходим всех ее соседей
for (int lmn = 0; lmn < 100; lmn++) {
if(EdgeNeighbours[ih][lmn]>=0)  knumb = EdgeNeighbours[ih][lmn];
else knumb=-1;
//------------------------------------------------------
J1count=0; 		JNCount=0;
if(knumb>=0) { for (int lk = 0; lk < Facenumb; lk++) { if(MFE1[lk][knumb] == 2) JNCount++; if(MFE1[lk][knumb] >= 1) J1count++; }
//-----------------------------------------------------
if(J1count > 1) JEdgeNeigh[ih][knn++] = JNCount;
else JEdgeNeigh[ih][knn+1]=-1; //knn++

}			}		}

//Вывод в файл JEdgeNeigh.txt
if(i==10) {
JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEdgeNeigh10.txt", ios::trunc);
for (int il = 0; il < Edgenumb; il++) {
for (int jl = 0; jl < 100; jl++) {
if(JEdgeNeigh[il][jl]>=0) JEdgeN << JEdgeNeigh[il][jl] << "\t";
}
JEdgeN <<"Edge number="<<il<<"\n";
}
JEdgeN.close();
}

for (int ip = 0; ip < 10; ip++)
for (int j = 0; j < 10; j++)
JEN[ip][j] = 0;

//Анализ матрицы JEdgeNeigh[i][j]
for (int il = 0; il < Edgenumb; il++) //для каждой грани
for (int nl = 0; nl < 10; nl++) //мы перебираем все варианты какой она может быть
if(JEdgeNeigh[il][0] == nl)  for (int mj = 1; mj < 100; mj++) //НАЧАЛО С 1, ЧТОБЫ НЕ УЧИТЫВАТЬ САМУ ГРАНЬ КАК СОСЕДА// и если она оказалась определенного типа, то мы перебираем всех ее соседей
for (int nk = 0; nk < 10; nk++) if(JEdgeNeigh[il][mj] == nk) JEN[nl][nk]++; //так что если сосед оказывается также определенного типа, то мы заносим их связь в матрицу JEN


//Полсчет "мощности" каждого тройного стыка с учетом соседей
for (int il = 0; il < Edgenumb; il++) {
TJpow[il]=0; NEneigh=0; for (int jl = 0; jl < 100; jl++)	if(JEdgeNeigh[il][jl]>=0) {TJpow[il]+= JEdgeNeigh[il][jl]; NEneigh++; };
TJpow[il] = TJpow[il]/NEneigh;
//Только ненулевые стыки
TJpow2[il]=0; NEneigh2=0; for (int jl = 0; jl < 100; jl++)	if(JEdgeNeigh[il][jl]>=0 && JEdgeNeigh[il][0]>0) {TJpow2[il]+= JEdgeNeigh[il][jl]; NEneigh2++; };
TJpow2[il] = TJpow2[il]/NEneigh2;
}


//Подсчет среднего по комплексу и дисперсии DISPERSION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NEneighSum=0.0; NEneighAv=0; for (int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) {NEneighSum += TJpow[il];  NEneighAv++;}
NEneighSum = NEneighSum/NEneighAv;

//Только ненулевые стыки
NEneighSum2=0.0; NEneighAv2=0; for (int il = 0; il < Edgenumb; il++) if(TJpow2[il]>0) {NEneighSum2 += TJpow2[il];  NEneighAv2++;}
NEneighSum2 = NEneighSum2/NEneighAv2;

//Dispersion
SNEneighSum=0.0; NEneighAv=0; for (int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) {SNEneighSum += (TJpow[il]*TJpow[il]); NEneighAv++;}
SNEneighSum = SNEneighSum/(NEneighAv);
//Только ненулевые стыки
SNEneighSum2=0.0; NEneighAv2=0; for (int il = 0; il < Edgenumb; il++) if(TJpow2[il]>0) {SNEneighSum2 += (TJpow2[il]*TJpow2[il]); NEneighAv2++;}
SNEneighSum2 = SNEneighSum2/(NEneighAv2);

if((SNEneighSum - NEneighSum*NEneighSum) >0) DispNES1 = sqrtl(SNEneighSum - NEneighSum*NEneighSum);
else DispNES1 = sqrtl(-SNEneighSum +NEneighSum*NEneighSum);
if((SNEneighSum2 - NEneighSum2*NEneighSum2) >0) DispNES2 = sqrtl(SNEneighSum2 - NEneighSum2*NEneighSum2);
else  DispNES2 = sqrtl(-SNEneighSum2 +NEneighSum2*NEneighSum2);

JEdgeN.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)TJspow.txt", ios::app);
JEdgeN << HAGBsFunc[i][0]<<"\t"<<25.0*powf(10,6)/(sqrtf(RHOD))<<"\t"<<30.0*powf(10,6)/(HAGBsFunc[i][1]*sqrtf(RHOD))<<"\t"<< 30.0*powf(10,6)*NEneighSum/(NEneighSum2*HAGBsFunc[i][1]*sqrtf(RHOD))<<"\t"<< 3.0*27.0*powf(10,6)/(NEneighSum2*HAGBsFunc[i][1]*sqrtf(RHOD));

//				if(NEneighSum2>0)  JEdgeN <<HAGBsFunc[i][0]<<"\t"<<NEneighSum<<"\t"<<NEneighSum2<<"\t"<<DispNES1<<"\t"<<DispNES2<<"\t"<<2.2*powl(NEneighSum,(-1.0/1.0))<<"\t"<<2.2*powl(NEneighSum2,(-1.0/1.0)); //<<"\t"<<SNEneighSum<<"\t"<<SNEneighSum2;
//				else JEdgeN <<HAGBsFunc[i][0]<<"\t"<<NEneighSum<<"\t"<<0.0;
//				 JEdgeN <<"Accumulated Strain =  "<<HAGBsFunc[i][0]<<"\n";
//				 JEdgeN<< NEneighSum << "\t";
//Мощность всех узлов
//				for(int il = 0; il < Edgenumb; il++) if(TJpow[il]>=0) JEdgeN<<TJpow[il] << "\t";
JEdgeN <<"\n";
JEdgeN.close();

//--------------------------------------------------------------------------------------

//Вывод в файл JEN.txt
JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)JEN.txt", ios::app);
JENStream <<AcStrain<<endl;
for (int i = 0; i < 10; i++) {
for (int j = 0; j < 10; j++) if(JEN[i][j]>=0)  JENStream << JEN[i][j] << "\t";
//			JENStream << JEN[i][j] << "\t";
JENStream <<"\n";
}
// Очевидно, что при таком алгоритме диагональные элементы учитываются дважды, а также считаем число всех элементов и делим потом на него, то есть
SumJEN=0; Jenii=0; Jenij=0;
for (int ic = 0; ic < 10; ic++)
for (int jc = 0; jc < 10; jc++) {SumJEN += 0.5*JEN[ic][jc]; if(ic==jc) Jenii+=0.5*JEN[ic][jc]; else Jenij+=0.5*JEN[ic][jc]; }
//cout<< SumJEN <<endl;

JENStream <<"\n"<<"\n"<<"\n";
JENStream.close();
//cout<<"(HAGBs)JEN.txt has been created"<<endl;

JENStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Jlk.txt", ios::app);
//				if(i==0)	JENStream << "STRAIN" << "\t"<<"ii"<<"\t"<<"ij"<<"\t"<<"00"<<"\t"<<"11"<<"\t"<<"22"<<"\t"<<"33"<<"\t"<<"01"<<"\t"<<"02"<<"\t"<<"03"<<"\t"<<"12"<<"\t"<<"13"<<"\t"<<"23"<<endl;
if(i==0)	JENStream << "STRAIN" << "\t"<<"ii"<<"\t"<<"ij"<<"\t"<<endl;
JENStream << HAGBsFunc[i][0]<<"\t"<<Jenii*100.0/SumJEN<<"\t"<< Jenij*100.0/SumJEN;
//				JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3] + JEN[4][4]<<"\t"<< JEN[0][1] + JEN[0][2] + JEN[0][3] + JEN[1][2] + JEN[1][3] + JEN[2][3] + JEN[2][4] + JEN[3][4] << "\t"<< JEN[0][0] << "\t"<< JEN[1][1] << "\t"<< JEN[2][2] << "\t"<< JEN[3][3] << "\t"<< JEN[0][1] << "\t"<< JEN[0][2] << "\t"<< JEN[0][3] << "\t"<< JEN[1][2] << "\t"<< JEN[1][3] << "\t"<< JEN[2][3];     //			 JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + JEN[0][2] + JEN[2][0] + JEN[0][3] + JEN[3][0] + JEN[1][2] + JEN[2][1] + JEN[1][3] + JEN[3][1] + JEN[2][3] + JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< JEN[0][2] + JEN[2][0]<< "\t"<< JEN[0][3] + JEN[0][3]<< "\t"<< JEN[1][2] + JEN[2][1]<< "\t"<< JEN[1][3] + JEN[3][1]<< "\t"<< JEN[2][3] + JEN[3][2];
//			 JENStream << HAGBsFunc[i][0]<<"\t"<<JEN[0][0] + JEN[1][1] + JEN[2][2] + JEN[3][3]<<"\t"<< JEN[0][1] + JEN[0][2] + JEN[0][3] + JEN[1][2] + JEN[1][3] + JEN[2][3] << "\t"<< JEN[0][1] << "\t"<< JEN[0][2] << "\t"<< JEN[0][3] << "\t"<< JEN[1][2] << "\t"<< JEN[1][3] << "\t"<< JEN[2][3]; 	// 			 JENStream << 0*JEN[0][0] + 2*JEN[1][1] + 4*JEN[2][2] + 6*JEN[3][3]<<"\t"<< JEN[0][1] + JEN[1][0] + 2*JEN[0][2] + 2*JEN[2][0] + 3*JEN[0][3] + 3*JEN[3][0] + 3*JEN[1][2] + 3*JEN[2][1] + 4*JEN[1][3] + 4*JEN[3][1] + 5*JEN[2][3] + 5*JEN[3][2]<< "\t"<< JEN[0][1] + JEN[1][0]<< "\t"<< 2*JEN[0][2] + 2*JEN[2][0]<< "\t"<< 3*JEN[0][3] + 3*JEN[0][3]<< "\t"<< 3*JEN[1][2] + 3*JEN[2][1]<< "\t"<< 4*JEN[1][3] + 4*JEN[3][1]<< "\t"<< 5*JEN[2][3] + 5*JEN[3][2];
JENStream <<"\n";
JENStream.close();
**/
