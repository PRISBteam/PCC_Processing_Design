///================================ DCC Edges types statistics, indices and configuration entropy =============================================================///
///==============================================================================================================================///
/** This subfile calculates the statistical characteristics of Edges incuding their Face indices and configurational entropy  **/
///==============================================================================================================================///

#include <iostream>
#include <fstream>
#include <vector>

using namespace std; ///Standard namespace

void EdgesStat(std::vector<unsigned int> const& CellNumbs, unsigned long numerator, map<unsigned int,unsigned int> SpecialCellMap, Eigen::SparseMatrix<double> const& FES, char* output_dir, double special_faces_fraction)
{
    vector<int> TJsTypes(CellNumbs.at(1)+1,0);
    map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;

    for (auto sit : SpecialCellMap)
        for(int k = 1; k < CellNumbs.at(1); k++) // Loop ovel all the
            if( FES.coeff(sit.first,k) == 1) TJsTypes.at(k)++;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
    /// Configurational Entropy related with Faces
    Configurational_Face_Entropy = - (j0* log2(j0) + j1* log2(j1) + j2* log2(j2) + j3* log2(j3));

    /// Data output
    /// Opening of the output streams
    string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::app);

    OutTJsFile << special_faces_fraction << "\t\t" << j0 << "\t\t" << j1 << "\t" << j2 << "\t" << j3 << endl;
    OutSFile << special_faces_fraction << "\t\t" << Configurational_Face_Entropy << endl;

    OutTJsFile.close();
    OutSFile.close();

    return;
}