///================================ DCC Edges types statistics, indices and configuration entropy =============================================================///
///==============================================================================================================================///
/** This subfile calculates the statistical characteristics of Edges incuding their Face indices and configurational entropy  **/
///==============================================================================================================================///

#include <iostream>
#include <fstream>
#include <vector>

using namespace std; ///Standard namespace

void EdgesStat(std::vector<unsigned int> const& CellNumbs, unsigned long numerator, Eigen::SparseMatrix<double> const& FacesGraph, Eigen::SparseMatrix<double> const& FES, char* output_dir)
{
    vector<int> TJsTypes(CellNumbs.at(1)+1,0);
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    for(int i = 1; i < numerator; i++)
        for(int j = 1; j < numerator; j++)
            if( i != j && FacesGraph.coeff(i,j) == 1) for(int k = 1; k < CellNumbs.at(1); k++) if( FES.coeff(j,k) == 1) TJsTypes.at(k)++;

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
    string TJs_output_filename = "TJsTypes.txt"s, Entropy_output_filename = "ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::app);

    OutTJsFile << j0 << "\t\t" << j1 << "\t" << j2 << "\t" << j3 << endl;
    OutSFile << Configurational_Face_Entropy << endl;

    OutTJsFile.close();
    OutSFile.close();

    return;
}