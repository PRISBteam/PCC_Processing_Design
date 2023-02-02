/// Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols) {

    Eigen::SparseMatrix<double> res(Rows, Cols);
    typedef Eigen::Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
            // cout << i << "\t" << j << "\t" << value << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong
//Sparse AB matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements and free the memory from the probe vector
    inAN.close();
    tripletList.clear();
    tripletList.shrink_to_fit();

    return res;
}

/// Function for reading triplet list from file
vector<Triplet<double>> TripletsReader(char* SMpath) {
    typedef Eigen::Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    double i = 0, j = 0, value = 0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i , j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return tripletList;
}

/// Function for reading tuples list from file
vector<vector<double>> VectorVectors4Reader(char* SMpath) {
    vector<vector<double>> res;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            double i = 0.0, j = 0.0, k = 0.0, l = 0;
            inAN >> i >> j >> k >> l;
            res.push_back({i,j,k,l});
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return res;
}

vector<tuple<double, double, double>> TuplesReader(char* SMpath) {
    typedef tuple<double, double, double> Tup; // Eigen library class
    std::vector<Tup> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(make_tuple(i, j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    return tripletList;
}
vector<tuple<double, double, double>> dTuplesReader(char* SMpath, unsigned int &size) {
    typedef tuple<double, double, double> Tup; // Eigen library class
    std::vector<Tup> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(make_tuple(i, j, value));
            //   cout << i << "\t" << j << "\t"<< value << "\t" << endl;
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong

    size = tripletList.size();
    return tripletList;
}

/// Creation int std::Vector from file
std::vector<unsigned int> VectorReader(char* FilePath) {
    std::vector<unsigned int> res;
    unsigned int i=0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else {std:cout << "The file " << FilePath << " cannot be read !" << endl; exit(1); } //If something goes wrong

    return res;
}

/// Creation double std::Vector from file
std::vector<double> dVectorReader(char* FilePath) {
    std::vector<double> res;
    double i = 0.0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else std:cout << "The file " << FilePath << " cannot be read" << endl; //If something goes wrong

    return res;
}

/// * The function count the number of Edges possessing with types J0, J1, J2, J3 and J4 for every junction * ///
/// * Calculation Face-Edge index ::                                                                                                     * ///

vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, vector<double> const& TJsTypes) {
    vector<double> res(100,0); /// Up to 100 types of possible TJs types

    for (unsigned int l = 0; l < FES.rows(); l++) // Loop over all Edges
        if (FES.coeff(l,face_number) == 1) res[TJsTypes.at(l)]++;
    // output in the form res[0] = #TJsTypes[0] incident to the face with the number face_number, res[1] = #TJsTypes[1] incident to the face with the number face_number,...

    return res;
}

/// * Function calculates the vector<int> "TJsTypes" of types TJs in the DCC using its FES incidence matrix and special faces sequence (s_faces_sequence) * ///
/// *                                                                                                                                                    * ///
vector<double> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES)
{
    vector<double> TJsTypes(CellNumbs.at(1),0); // CellNumbs.at(1) is the number of Edges
    for (auto vit: s_faces_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(1); k++) // loop over all Edges
            if (FES.coeff(k, vit) == 1) TJsTypes.at(k)++;

    return TJsTypes;
}

/// EntropyIncreaseList
vector<double> Get_EntropyIncreaseList(std::vector<unsigned int> &S_Vector, vector<double> const &TJsTypes, SpMat const &FES) {
//vector<double> Get_EntropyIncreaseList(std::vector<unsigned int> &S_Vector, double> const &TJsTypes, SpMat const &FES) {
    vector<double> EntropyIncreaseList(CellNumbs.at(2), 0); // vector with values of configuration entropy increases at conversion of each Face
    /// amounts of TJs of different types
    double Jall = 0, J0all = 0, JNall = 0;
    double  J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
    double J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
    double  J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
    double  J0 = CellNumbs.at(1) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
    Jall = J0 + J1 + J2 + J3;

    for (unsigned int k = 0; k < CellNumbs.at(2); k++) { /// loop over all Faces in DCC
        double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_Entropy = 0, CFace_EntropyIncrease = 0;

        if (S_Vector.at(k) == 0) { // Loop over each still ORDINARY element neighbours
            J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
            /// j_types_neigh_fractions calculation ///***function GBIndex(k, FES, TJsTypes) from DCC_SupportFunctions.h
            // output in the form j_types_neigh_fractions[0] = #TJsTypes[0] incident to the face with the number face_number, j_types_neigh_fractions[1] = #TJsTypes[1] incident to the face with the number face_number,...
            vector<double> j_types_neigh_fractions = GBIndex(k, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
//REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;

            /// For a particular Face face_number = k (!)
            /// Values before conversion
            J00 = j_types_neigh_fractions.at(0);
            J10 = j_types_neigh_fractions.at(1);
            J20 = j_types_neigh_fractions.at(2);
            J0all = J00 + J10 + J20;

 /*
            Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
            j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
            double j0s = j0, j1s = j1, j2s = j2, j3s = j3;

            if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!
            /// Configuration Entropy related with Faces
            Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

*/


//REPAIR// cout << " J00= " << J00<< " J10= " << J10 << " J20= " << J20 << endl;
            // Values after conversion
            J1N = J00;
            J2N = J10;
            J3N = J20;
            JNall = J1N + J2N + J3N;
            // The entropy increase calculation for a given Face
            // Conversion from numbers to fractions
            double  j0 = 0.0, j1 = 0.0,  j2 = 0.0, j3 = 0.0, l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0, j1n = 0.0,  j2n = 0.0, j3n = 0.0,  j00 = 0.0,  j10 = 0.0,  j20 = 0.0, l2j0_0 = 0.0, l2j1_0 = 0.0, l2j2_0 = 0.0, l2j0_n = 0.0, l2j1_n = 0.0, l2j2_n = 0.0, l2j3_n = 0.0;

            if (J0all > 0) j00 = (double) J00 / J0all;
            if (j0 - j00 > 0)  l2j0_n = log2(j0 - j00) ;
            if (J0all > 0) j10 = (double) J10 / J0all;
            if (j1 - j10 + j1n > 0) l2j1_n = log2(j1 - j10 + j1n);
            if (J0all > 0) j20 = (double) J20 / J0all;
            if (j2 - j20 + j2n > 0) l2j2_n = log2(j2 - j20 + j2n);
            if (j3 + j3n > 0) l2j3_n = log2(j3 + j3n);

            if (JNall > 0) j1n = (double) J1N / JNall;
            if (j1n > 0) l2j1_n = log2(j1n);
            if (JNall > 0) j2n = (double) J2N / JNall;
            if (j2n > 0) l2j2_n = log2(j2n);
            if (JNall > 0) j3n = (double) J3N / JNall;
            if (j3n > 0) l2j3_n = log2(j3n);

            if (Jall > 0) j0 = (double) J0 / Jall;
            if (j0 > 0) l2j0 = log2(j0);
            if (Jall > 0) j1 = (double) J1 / Jall;
            if (j1 > 0) l2j1 = log2(j1);
            if (Jall > 0) j2 = (double) J2 / Jall;
            if (j2 > 0) l2j2 = log2(j2);
            if (Jall > 0) j3 = (double) J3 / Jall;
            if (j3 > 0) l2j3 = log2(j3);

//            CFace_EntropyIncrease = (j0 * l2j0 - j00 * l2j0_0)  + (j1 * l2j1 + j1n * l2j1_n - j10 * l2j1_0)  + (j2 * l2j2 + j2n * l2j2_n - j20 * l2j2_0)  + (j3 * l2j3 + j3n * l2j3_n);

            CFace_EntropyIncrease = - (j0 * l2j0 + j1 * l2j1 + j2 * l2j2 + j3 * l2j3
                                    - (j0 - j00)* l2j0_n
                                    - (j1 - j10 + j1n)* l2j1_n
                                    - (j2 - j20 + j2n)* l2j2_n
                                    - (j3 + j3n)* l2j3_n);

            //REPAIR    cout  << "\t\t j00 * l2j0_0: " <<  (j0 * l2j0 - j00 * l2j0_0) + (j1 * l2j1 + j1n * l2j1_n - j10 * l2j1_0) + (j2 * l2j2 + j2n * l2j2_n - j20 * l2j2_0) << "\t\t" << endl;

/// The result of one iteration (EntropyIncreaseList value for a particular Face face_number = k)
            EntropyIncreaseList.at(k) = CFace_EntropyIncrease;
        } // if OrdinaryCells (S_Vector.at(Face) == 0)
    } // for (..k < CellNumbs.at(2)..)

    return EntropyIncreaseList;
}

/// DDRX support function :: GFS matrix reading and calculation of new seeds at the centres of GBs
tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates) {
    tuple <double, double, double> res; // find two grain neighbour for fnumber
    vector<double> xx, yy, zz;
    //Triplet<double> res;     // find two grain neighbour for fnumber

    SpMat GFS(CellNumbs.at(2),CellNumbs.at(3));
    /// GFS matrix reading
    GFS = SMatrixReader(paths.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains

    vector<int> Grain_neighbours;
    vector<unsigned int> grainIDs;
#pragma omp parallel for // parallel execution by OpenMP
  /*
        grainIDs.clear();
        for (SparseMatrix<double>::InnerIterator it(GFS, facenumb); it; ++it)
            grainIDs.push_back(it.col());
*/
    //normal loop
for (unsigned int j = 0; j < CellNumbs.at(3); ++j) {
        if (GFS.coeff(facenumb, j) == 1) {
            grainIDs.push_back(j);
            if(j<CellNumbs.at(3)-1) {
             for (unsigned int k = j + 1; k < CellNumbs.at(3); ++k)
                if (GFS.coeff(facenumb, k) == 1)
                    grainIDs.push_back(k);
            }
        }
}
//        cout << grainIDs[0] << " " << grainIDs[1] << endl;

    if(grainIDs.size()>0) xx.push_back(get<0>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) xx.push_back(get<0>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>0) yy.push_back(get<1>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) yy.push_back(get<1>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>0) zz.push_back(get<2>(AllSeeds_coordinates.at(grainIDs[0])));
    if(grainIDs.size()>1) zz.push_back(get<2>(AllSeeds_coordinates.at(grainIDs[1])));
    if(grainIDs.size()>1) res = make_tuple(0.5*(xx[0] + xx[1]), 0.5*(yy[0] + yy[1]), 0.5*(zz[0] + zz[1]));
        else if(grainIDs.size()>0) res = make_tuple(xx[0], yy[0], zz[0]);
        else res = make_tuple(0, 0, 0);
//    cout << "facenumb " << facenumb << " " << 0.5*(xx[0] + xx[1]) << "\t" << 0.5*(yy[0] + yy[1]) << "\t"<< 0.5*(zz[0] + zz[1]) << "\t" << endl;
//    cout << "facenumb " << facenumb << " " << get<0>(res) << "\t" << get<1>(res) << "\t"<< get<2>(res) << "\t" << endl;

    return res;

}

bool is_file_exists(const string fileName)
{
    char* charfileName = const_cast<char*>(fileName.c_str());
    std::ifstream infile(charfileName);
    return infile.good();
}

/// Configuration TJs entropy
double Get_TJsEntropy(vector<unsigned int> special_faces_seq) {
double TJsEntropy = 0.0;

    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    vector<double> TJsTypes;

    SpMat FES(CellNumbs.at(1), CellNumbs.at(2));
    FES = SMatrixReader(paths.at(5 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces

    TJsTypes = EdgesTypesCalc(CellNumbs, special_faces_seq, FES);

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = (double) CellNumbs.at(1);

/// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
    double j0s = j0, j1s = j1, j2s = j2, j3s = j3;

    /// using values with pow(10,-10) instead of 0s!
    if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!

    /// Configuration Entropy related with Faces
    TJsEntropy = - (j0s + j1s + j2s + j3s);

return TJsEntropy;
};

/// ARCHIVE ///
//Erase First Occurrence of given  substring from main string
void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

/*
 *     double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0, Face_Entropy_Median = 0 , Face_Entropy_Skrew = 0;

 */

/*
 *     /// Ordinary face sequence
    ordinary_faces_sequence.clear();
    for (auto  itr = State_sVector.begin(); itr != State_sVector.end(); ++itr)
        if (*itr == 0) ordinary_faces_sequence.push_back(distance(State_sVector.begin(), itr));

 */