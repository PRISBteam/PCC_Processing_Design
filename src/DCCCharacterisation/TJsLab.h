///================================ DCC Edges types statistics, indices and configuration entropy =============================================================///
///==============================================================================================================================///
/** This subfile calculates the statistical characteristics of Edges incuding their Face indices and configurational entropy  **/
///==============================================================================================================================///

using namespace std; ///Standard namespace
map<unsigned int, unsigned int> EdgesStat(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> const& CellNumbs, Eigen::SparseMatrix<double> const& FES, char* output_dir)
{
    vector<int> TJsTypes(CellNumbs.at(1),0);
    map<unsigned int, unsigned int> res; // Here 100 is an arbitrary number of Edge types
    map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0, Face_Entropy_Median = 0 , Face_Entropy_Skrew = 0;

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
    if (j0s == 0) j0s = pow(10,-30); if (j1s == 0) j1s = pow(10,-30); if (j2s == 0) j2s = pow(10,-30); if (j3s == 0) j3s = pow(10,-30); //Gives 0 in entropy!

    /// Configuration Entropy related with Faces
    Configurational_Face_Entropy = - (j0s* log2(j0s) + j1s* log2(j1s) + j2s* log2(j2s) + j3s* log2(j3s));

    /// Median part in the entropy decomposition
    vector<double> j_types_fractions = {j0s, j1s, j2s, j3s}; /// using values with pow(10,-10) instead of 0s!
    Face_Entropy_Median = -(1.0/j_types_fractions.size())*log2(j0s*j1s*j2s*j3s);

    /// Skrew part (divergence from the uniform distribution -> S_max) in the entropy decomposition
    for (int j = 0; j < j_types_fractions.size(); j++)
        for (int i = 0; i < j; i++)
            Face_Entropy_Skrew += -(1.0/j_types_fractions.size())*(j_types_fractions[i]-j_types_fractions[j])*log2(j_types_fractions[i]/j_types_fractions[j]);

    /// Data output
    /// Opening of the output streams
    string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::app);
    double special_faces_fraction =(double) s_faces_sequence.size()/ CellNumbs.at(2);
//    cout << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t"<< Face_Entropy_Median << endl;

    OutTJsFile << special_faces_fraction << "\t\t" << j0 << "\t\t" << j1 << "\t" << j2 << "\t" << j3 << endl;
    OutSFile << special_faces_fraction << "\t\t" << Configurational_Face_Entropy << "\t\t" << Face_Entropy_Median << "\t\t" << Face_Entropy_Skrew << endl;

    OutTJsFile.close();
    OutSFile.close();

    int b = 0;
    for (unsigned int it : TJsTypes) res[b++] = it;

    return res;
}
