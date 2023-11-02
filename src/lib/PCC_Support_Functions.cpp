/// Function common for all other modules in the code

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <Eigen/SparseCore>

using namespace std; // standard namespace
using namespace Eigen; // standard namespace

#include "PCC_Support_Functions.h" // It must be here - first in this list (!)
/// # 1 # Checking if file exists in the directory 'fileName'

bool is_file_exists(const string fileName) {
    char* charfileName = const_cast<char*>(fileName.c_str());
    std::ifstream infile(charfileName);
    return infile.good();
} // END of is_file_exists()

std::vector<unsigned int> VectorIReader(const char* FilePath) { // creation int std::Vector from file
std::vector<unsigned int> res; unsigned int i = 0;  // function output
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { // if the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else { cout << "The file " << FilePath << " cannot be read !" << endl; exit(1); } // if something goes wrong
    return res;
} // END of VectorIReader()

std::vector<double> VectorDReader(const char* FilePath) { // creation double std::vector from file
std::vector<double> res; double d = 0.0;  // function output
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { // if the file was successfully open, then
        while(inCellNumbers >> d) res.push_back(d);
    } else std:cout << "The file " << FilePath << " cannot be read" << endl; // if something goes wrong
    return res;
} // END of VectorDReader()

/// Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols) {
    /// function output
    Eigen::SparseMatrix<double> res(Rows, Cols);

    typedef Eigen::Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    double i = 0.0, j = 0.0, value = 0.0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
// REPAIR cout << i << "\t" << j << "\t" << value << endl;
        }
    } else cout << "WARNING: The file " << SMpath << " cannot be read" << endl; //If something goes wrong
//Sparse AB matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements and free the memory from the probe vector
    inAN.close();
    tripletList.clear();
    tripletList.shrink_to_fit();

    return res;
}

/// Function for reading triplet list from file
vector<Eigen::Triplet<double>> TripletsReader(char* SMpath) {
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


/// DDRX support function :: GFS matrix reading and calculation of new seeds at the centres of GBs

std::tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates) {
    tuple <double, double, double> res; // find two grain neighbour for fnumber
    vector<double> xx, yy, zz;
    //Triplet<double> res;     // find two grain neighbour for fnumber

    Eigen::SparseMatrix<double> GFS(CellNumbs.at(2),CellNumbs.at(3));
    /// GFS matrix reading
    GFS = SMatrixReader(paths.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains

    vector<int> Grain_neighbours;
    vector<unsigned int> grainIDs;
// #pragma omp parallel for // parallel execution by OpenMP

//          grainIDs.clear();
//          for (SparseMatrix<double>::InnerIterator it(GFS, facenumb); it; ++it)
//              grainIDs.push_back(it.col());

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

} /// END of std::tuple<double, double, double> find_aGBseed() function



/**
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

/// * The function count the number of Edges possessing with types J0, J1, J2, J3 and J4 for every junction * ///
/// * Calculation Face-Edge index ::                                                                                                     * ///

vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, vector<double> const& TJsTypes) {
    vector<double> res(100,0); /// Up to 100 types of possible TJs types

    for (unsigned int l = 0; l < FES.rows(); l++) // Loop over all Edges
        if (FES.coeff(l,face_number) != 0)
            res[TJsTypes.at(l)]++;
    // output in the form res[0] = #TJsTypes[0] incident to the face with the number face_number, res[1] = #TJsTypes[1] incident to the face with the number face_number,...

    return res;
}

/// * Function calculates the vector<int> "TJsTypes" of types TJs in the PCC using its FES incidence matrix and special faces sequence (s_faces_sequence) * ///
/// *                                                                                                                                                    * ///
vector<double> NodesTypesCalc(std::vector<unsigned int> const &CellNumbs, vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &ENS)
{
vector<double> TJsTypes(CellNumbs.at(0),0); // CellNumbs.at(1) is the number of Edges
for (auto vit: s_faces_sequence) // loop over all Special Faces
for(int k = 0; k < CellNumbs.at(0); k++) // loop over all Edges
if (ENS.coeff(k, vit) != 0)
TJsTypes.at(k)++;

return TJsTypes;
}

vector<double> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES)
{
vector<double> TJsTypes(CellNumbs.at(1),0); // CellNumbs.at(1) is the number of Edges
for (auto vit: s_faces_sequence) // loop over all Special Faces
for(int k = 0; k < CellNumbs.at(1); k++) // loop over all Edges
if (FES.coeff(k, vit) != 0)
TJsTypes.at(k)++;

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

//                       Jall = CellNumbs.at(1);
           // Conversion from numbers to fractions
           // (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
//                       j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
//                       double j0s = j0, j1s = j1, j2s = j2, j3s = j3;

//                       if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!
                       /// Configuration Entropy related with Faces
//                       Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);


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

std::vector<vector<int>> Get_cases_list(std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat const &FES, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, double const &p_index) {
std::vector<vector<int>> cases_list; // in every case its own TJs vector

std::vector<int> NewEdgeTypes = EdgeTypes;

if (p_index == 0) { // cases: direct assigment of special faces one by one
unsigned int iter = 0;
for (unsigned int f = 0; f < CellNumbs.at(2 + (dim - 3)); ++f) { // loop over all Faces in PCC
NewEdgeTypes = EdgeTypes; // on each step it is equal to the initial Edge state defined by S_Vector

/// only ORDINARY f will be taken in the "cases list" (!) as if they change their type to special
if (S_Vector.at(f) == 0) { // Loop over each still ORDINARY cell neighbours
for(unsigned int e = 0; e < CellNumbs.at(1 + (dim - 3)); ++e) // loop over all Edges
if (FES.coeff(e, f) != 0) NewEdgeTypes.at(e)++;
//                cout << " J1: " << std::count(NewEdgeTypes.begin(), NewEdgeTypes.end(), 1)<< " J2: " << std::count(NewEdgeTypes.begin(), NewEdgeTypes.end(), 2) << " J3: " << std::count(NewEdgeTypes.begin(), NewEdgeTypes.end(), 3) << endl; // containing 1 incident special face

cases_list.push_back(NewEdgeTypes);
/// map from the Cases to special Faces set
cases_to_sfaces.insert( std::pair<unsigned int, std::vector<unsigned int>> (iter, {f})); // each cases coresponds to the vector with a single element which is the same number of special face
// new # in the case_list
// REPAIR            cout << " case iterator: " <<  iter << "   " << cases_to_sfaces[iter].at(0) << endl;
++iter;
} // end of if (S_Vector.at(f) == 0)
} // end for (unsigned int f = 0; f < CellNumbs.at(2); ++f)
} else if (p_index == 1) { // cases: crystallographic restrictions with "grain rotations"
///...
}

return cases_list;
} // END of Get_cases_list()

std::vector<double> dq_upleft(std::vector<double> const &grain_q_coord, double dq_step) { // function 1
std::vector<double> grain_upleft_qcoord = grain_q_coord;

double new_grain_q_coord;
// x (Q1) - 2dx
new_grain_q_coord = grain_q_coord.at(0) - 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_upleft_qcoord.at(0) = new_grain_q_coord;
else return grain_q_coord;
// y (Q2) + 2dy
new_grain_q_coord = grain_q_coord.at(1) + 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_upleft_qcoord.at(1) = new_grain_q_coord;
else return grain_q_coord;
// z (Q3) = const

return grain_upleft_qcoord;
}

std::vector<double> dq_upright(vector<double> const &grain_q_coord, double dq_step) { // function 2
std::vector<double> grain_upright_qcoord = grain_q_coord;

// x (Q1) = const
double new_grain_q_coord;

// y + 2dy
new_grain_q_coord = grain_q_coord.at(1) + 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_upright_qcoord.at(1) = new_grain_q_coord;
else return grain_q_coord;
// z - 2dz
new_grain_q_coord = grain_q_coord.at(2) - 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_upright_qcoord.at(2) = new_grain_q_coord;
else return grain_q_coord;

return grain_upright_qcoord;
}

std::vector<double> dq_downleft(vector<double> const &grain_q_coord, double dq_step) { // function 3
std::vector<double> grain_downleft_qcoord = grain_q_coord;

double new_grain_q_coord;
// x (Q1) = const

// y - 2dy
new_grain_q_coord = grain_q_coord.at(1) - 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_downleft_qcoord.at(1) = new_grain_q_coord;
else return grain_q_coord;

// z + 2dz
new_grain_q_coord = grain_q_coord.at(2) + 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_downleft_qcoord.at(2) = new_grain_q_coord;
else return grain_q_coord;

return grain_downleft_qcoord;
}

std::vector<double> dq_downright(vector<double> const &grain_q_coord, double dq_step) { // function 3
std::vector<double> grain_downright_qcoord = grain_q_coord;

double new_grain_q_coord;

// x + 2dx
new_grain_q_coord = grain_q_coord.at(0) + 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_downright_qcoord.at(0) = new_grain_q_coord;
else return grain_q_coord;

// y - 2dy
new_grain_q_coord = grain_q_coord.at(1) - 2.0* dq_step;
if (new_grain_q_coord >= 0 && new_grain_q_coord <= 1)
grain_downright_qcoord.at(1) = new_grain_q_coord;
else return grain_q_coord;

// z (Q3) = const

return grain_downright_qcoord;
}


/// Crystallographic BCC disorientations
double Get_2grains_FCCdisorientation(std::vector<double> &grain_quaternion1, std::vector<double> &grain_quaternion2){ /// in degrees (!)
    double disorienation_GBs_angle; // function output

    //DMat Cubic_symmetry;
    Eigen::Matrix<double, 24, 4> Cubic_symmetry;

    Cubic_symmetry <<   0, 0, 0, 1,
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            1.0/sqrt(2.0), 0, 0, 1.0/sqrt(2.0),
            0, 1.0/sqrt(2.0), 0, 1.0/sqrt(2.0),
            0, 0, 1.0/sqrt(2.0), 1.0/sqrt(2.0),
            -1.0/sqrt(2.0), 0, 0, 1.0/sqrt(2.0),
            0, -1.0/sqrt(2.0), 0, 1.0/sqrt(2.0),
            0, 0, -1.0/sqrt(2.0), 1.0/sqrt(2.0),
            1.0/sqrt(2.0), 1.0/sqrt(2.0), 0, 0,
            -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0, 0,
            0, 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0,
            0, -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0,
            1.0/sqrt(2.0), 0, 1.0/sqrt(2.0), 0,
            -1.0/sqrt(2.0), 0, 1.0/sqrt(2.0), 0,
            0.5, 0.5, 0.5, 0.5,
            -0.5, -0.5, -0.5, 0.5,
            0.5, -0.5, 0.5, 0.5,
            -0.5, 0.5, -0.5, 0.5,
            -0.5, 0.5, 0.5, 0.5,
            0.5, -0.5, -0.5, 0.5,
            -0.5, -0.5, 0.5, 0.5,
            0.5, 0.5, -0.5, 0.5;

    Eigen::Matrix<double, 4, 4> face_misorientations_mat1;
    face_misorientations_mat1 << grain_quaternion1.at(0), grain_quaternion1.at(1),  grain_quaternion1.at(2), grain_quaternion1.at(3),
            -grain_quaternion1.at(1),  grain_quaternion1.at(0),  grain_quaternion1.at(3), -grain_quaternion1.at(2),
            -grain_quaternion1.at(2),  -grain_quaternion1.at(3), grain_quaternion1.at(0), grain_quaternion1.at(1),
            -grain_quaternion1.at(3),  grain_quaternion1.at(2), -grain_quaternion1.at(1), grain_quaternion1.at(0);
//
    Eigen::Matrix<double, 4, 1> qat2; // coloumn vector
    qat2 << grain_quaternion2.at(0),
            grain_quaternion2.at(1),
            grain_quaternion2.at(2),
            grain_quaternion2.at(3);

    Eigen::Matrix<double, 4, 1> face_misorientations_vect1 = face_misorientations_mat1 * qat2;

    Eigen::Matrix<double, 4, 4> face_misorientations_mat2;
    face_misorientations_mat2 << grain_quaternion2.at(0), grain_quaternion2.at(1),  grain_quaternion2.at(2), grain_quaternion2.at(3),
            -grain_quaternion2.at(1),  grain_quaternion2.at(0),  grain_quaternion2.at(3), -grain_quaternion2.at(2),
            -grain_quaternion2.at(2),  -grain_quaternion2.at(3), grain_quaternion2.at(0), grain_quaternion2.at(1),
            -grain_quaternion2.at(3),  grain_quaternion2.at(2), -grain_quaternion2.at(1), grain_quaternion2.at(0);
//
    Eigen::Matrix<double, 4, 1> qat1; // coloumn vector
    qat1 << grain_quaternion1.at(0),
            grain_quaternion1.at(1),
            grain_quaternion1.at(2),
            grain_quaternion1.at(3);

    Eigen::Matrix<double, 4, 1> face_misorientations_vect2 = face_misorientations_mat2 * qat1;

    Eigen::Matrix<double, 4, 4> sym_1m, sym_11m;
    Eigen::Matrix<double, 4, 1> sym_1v, sym_11v, sym_2v, sym_22v;

    disorienation_GBs_angle = 180.0;
/// loop over all BCC symmetries
    for(int i = 0; i < 24; ++i) {
        sym_1m << Cubic_symmetry(i,0), -Cubic_symmetry(i,1), -Cubic_symmetry(i,2), -Cubic_symmetry(i,3),
                Cubic_symmetry(i,1), Cubic_symmetry(i,0), -Cubic_symmetry(i,3), Cubic_symmetry(i,2),
                Cubic_symmetry(i,2), Cubic_symmetry(i,3), Cubic_symmetry(i,0), -Cubic_symmetry(i,1),
                Cubic_symmetry(i,3), -Cubic_symmetry(i,2), Cubic_symmetry(i,1), Cubic_symmetry(i,0);

        sym_1v = sym_1m * face_misorientations_vect1;

        sym_11m <<  Cubic_symmetry(i,0), -Cubic_symmetry(i,1), -Cubic_symmetry(i,2), -Cubic_symmetry(i,3),
                Cubic_symmetry(i,1), Cubic_symmetry(i,0), -Cubic_symmetry(i,3), Cubic_symmetry(i,2),
                Cubic_symmetry(i,2), Cubic_symmetry(i,3), Cubic_symmetry(i,0), -Cubic_symmetry(i,1),
                Cubic_symmetry(i,3), -Cubic_symmetry(i,2), Cubic_symmetry(i,1), Cubic_symmetry(i,0);

        sym_11v = sym_11m * face_misorientations_vect2;

//            sym_1v = sym_1v.transpose();
//            sym_11v = sym_11v.transpose();

///
        Eigen::Matrix<double, 4, 4> sym_2_const;
        sym_2_const << sym_1v(0), -sym_1v(1), -sym_1v(2), -sym_1v(3),
                sym_1v(1), sym_1v(0), -sym_1v(3), sym_1v(2),
                sym_1v(2), sym_1v(3), sym_1v(0), -sym_1v(1),
                sym_1v(3), -sym_1v(2), sym_1v(1), sym_1v(0);

        Eigen::Matrix<double, 4, 4> sym_22_const;
        sym_22_const << sym_11v(0), -sym_11v(1), -sym_11v(2), -sym_11v(3),
                sym_11v(1), sym_11v(0), -sym_11v(3), sym_11v(2),
                sym_11v(2), sym_11v(3), sym_11v(0), -sym_11v(1),
                sym_11v(3), -sym_11v(2), sym_11v(1), sym_11v(0);

/// another loop over all BCC symmetries
        for(int j = 0; j < 24; ++j) {
            Eigen::Matrix<double, 4, 1> sym_csj;
            sym_csj << Cubic_symmetry(j,0),
                    -Cubic_symmetry(j,1),
                    -Cubic_symmetry(j,2),
                    -Cubic_symmetry(j,3);

            sym_2v = sym_2_const * sym_csj;
            sym_22v = sym_22_const * sym_csj;

            double disorientaion = std::min(2.0*acos(sym_2v(0)), 2.0*acos(sym_22v(0)));
            if(disorientaion < disorienation_GBs_angle)
                disorienation_GBs_angle = disorientaion;

        } // end of for(i = 0; i < 24; ++i)
// REPAIR cout << "disorientation  " << disorienation_GBs_angle*57.3 << endl;

    } // end of for(j = 0; j < 24; ++j)

    return disorienation_GBs_angle*57.3; /// in degrees (!)
} // END of double Get_2grains_FCCdisorientation()

/// isHAGB checker
bool isHAGB(std::vector<double> &grain_quaternion, std::vector<double> &new_grain_quaternion, double &threshold){
    bool isHAGB = 0;

// Obtain disorientation
    double angle = Get_2grains_FCCdisorientation(grain_quaternion, new_grain_quaternion);
// Check the threshold
    if(angle > threshold)
        isHAGB = 1;

    return isHAGB;
} // END of isHAGB()

std::vector<vector<int>> Get_crystallographic_cases_list(std::vector<vector<double>> &grain_quaternions_list, std::map<unsigned int, std::vector<unsigned int>> &g_gbs_map, double &dq, double &HAGBs_threshold, std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat &GFS, SpMat &FES, std::map<unsigned int, unsigned int> &cases_to_grains, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, std::map<unsigned int, std::vector<double>> &cases_to_new_quaternions, double const &p_index) {
    std::vector<vector<int>> cryst_cases_list; // function output

/// Assighment NEW orientations and related 3-cases
    std::vector<int> NewEdgeTypes = EdgeTypes;

    std::vector<double> new_grain_triangle_quaternions(4);
    std::vector<unsigned int> gb_special_set_1, gb_special_set_2, gb_special_set_3, gb_special_set_4; // cases for grain rotations
    std::vector<double> grain_q_quaternion(3), new_grain_q_quaternion(3), new_grain_quaternion(4);
    double q0_coord, new_q0_coord;

/// Studying cases - through the whole quaternion list

    unsigned int iter = 0; // count CASES
    for (auto  g_itr = grain_quaternions_list.begin(); g_itr != grain_quaternions_list.end(); ++g_itr) {
        gb_special_set_1.clear(); gb_special_set_2.clear();
        gb_special_set_3.clear(); gb_special_set_4.clear();

        std::vector<double> grain_full_quaternion = *g_itr; // a grain (or case) full quaternion from the list
        q0_coord = grain_full_quaternion.at(0);
        grain_q_quaternion = {pow(grain_full_quaternion.at(1),2)/(1.0 - pow(q0_coord,2)), pow(grain_full_quaternion.at(2),2)/(1.0 - pow(q0_coord,2)), pow(grain_full_quaternion.at(3),2)/(1.0 - pow(q0_coord,2))};

/// GBs set for a grain
        std::vector<unsigned int> gb_set = g_gbs_map[distance(grain_quaternions_list.begin(),g_itr)];

        /// Three cases
        // case 1 dq_down - function 3
        new_grain_q_quaternion = dq_upleft(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_1.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_1)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation1 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;

        // case 2
        new_grain_q_quaternion = dq_upright(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_2.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_2)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation2 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;
        // case 3
        new_grain_q_quaternion = dq_downleft(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_3.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_3)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation3 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;
        // case 4
        new_grain_q_quaternion = dq_downright(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

// REPAIR    cout << distance(grain_quaternions_list.begin(),g_itr) << endl;

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_4.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_4)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation3 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;

// REPAIR if (gb_special_set_1.size() != 0 || gb_special_set_2.size() != 0 || gb_special_set_3.size() != 0 || gb_special_set_4.size() != 0) cout << "q_check " << grain_q_quaternion.at(0) + grain_q_quaternion.at(1) + grain_q_quaternion.at(2) << " ngb " << gb_special_set_1.size() << "  " << gb_special_set_2.size() << "  " << gb_special_set_3.size() << "  " << gb_special_set_4.size() <<endl;
//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;

/// New State Vectors for Edges
        if (gb_special_set_1.size() != 0 || gb_special_set_2.size() != 0 || gb_special_set_3.size() != 0 || gb_special_set_4.size() != 0) {
// case 1
            NewEdgeTypes = EdgeTypes;
            for (unsigned int gb: gb_special_set_1) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    for (unsigned int e = 0; e < CellNumbs.at(1 ); ++e) // loop over all Edges
                        if (FES.coeff(e, gb) != 0) NewEdgeTypes.at(e)++;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewEdgeTypes);

// case 2
            NewEdgeTypes = EdgeTypes;
            for (unsigned int gb: gb_special_set_2) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    for (unsigned int e = 0; e < CellNumbs.at(1 ); ++e) // loop over all Edges
                        if (FES.coeff(e, gb) != 0) NewEdgeTypes.at(e)++;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewEdgeTypes);

// case 3
            NewEdgeTypes = EdgeTypes;
            for (unsigned int gb: gb_special_set_3) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    for (unsigned int e = 0; e < CellNumbs.at(1 ); ++e) // loop over all Edges
                        if (FES.coeff(e, gb) != 0) NewEdgeTypes.at(e)++;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewEdgeTypes);

// case 4
            NewEdgeTypes = EdgeTypes;
            for (unsigned int gb: gb_special_set_4) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    for (unsigned int e = 0; e < CellNumbs.at(1 ); ++e) // loop over all Edges
                        if (FES.coeff(e, gb) != 0) NewEdgeTypes.at(e)++;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewEdgeTypes);

        } // end of if (gb_special_set_1.size() != 0 ||...

    } // end of for (auto  g_itr = grain_quaternions_list.begin(); g_itr != grain_quaternions_list.end(); ++g_itr)

    return cryst_cases_list; // function output
} // END of Get_crystallographic_cases_list()

std::vector<vector<int>> Get_crystallographic_cases_random_list(std::vector<vector<double>> &grain_quaternions_list, std::map<unsigned int, std::vector<unsigned int>> &g_gbs_map, double &dq, double &HAGBs_threshold, std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat &GFS, SpMat &FES, std::map<unsigned int, unsigned int> &cases_to_grains, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, std::map<unsigned int, std::vector<double>> &cases_to_new_quaternions, double const &p_index) {
    std::vector<vector<int>> cryst_cases_list; // function output

/// Assighment NEW orientations and related 3-cases
    std::vector<int> NewSpecialFaces = S_Vector;

    std::vector<double> new_grain_triangle_quaternions(4);
    std::vector<unsigned int> gb_special_set_1, gb_special_set_2, gb_special_set_3, gb_special_set_4; // cases for grain rotations
    std::vector<double> grain_q_quaternion(3), new_grain_q_quaternion(3), new_grain_quaternion(4);
    double q0_coord, new_q0_coord;

/// Studying cases - through the whole quaternion list
// cout << "cases_random_list " << "place 1" << endl;
    unsigned int iter = 0; /// count CASES
    for (auto  g_itr = grain_quaternions_list.begin(); g_itr != grain_quaternions_list.end(); ++g_itr) {

        gb_special_set_1.clear(); gb_special_set_2.clear();
        gb_special_set_3.clear(); gb_special_set_4.clear();

        std::vector<double> grain_full_quaternion = *g_itr; // a grain (or case) full quaternion from the list
        q0_coord = grain_full_quaternion.at(0);
        grain_q_quaternion = {pow(grain_full_quaternion.at(1),2)/(1.0 - pow(q0_coord,2)), pow(grain_full_quaternion.at(2),2)/(1.0 - pow(q0_coord,2)), pow(grain_full_quaternion.at(3),2)/(1.0 - pow(q0_coord,2))};

//cout << "CHECK: " << grain_q_quaternion.at(0) +  grain_q_quaternion.at(1) +  grain_q_quaternion.at(2) << "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;

//cout << "cases_random_list " << "place 2: " << distance(grain_quaternions_list.begin(),g_itr) << endl;

/// GBs set for a grain
        std::vector<unsigned int> gb_set = g_gbs_map[distance(grain_quaternions_list.begin(),g_itr)];
//        cout << "cases_random_list " << "place 3: " << distance(grain_quaternions_list.begin(),g_itr) << endl;

        /// Three cases
        // case 1 dq_down - function 3
        new_grain_q_quaternion = dq_upleft(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

//cout << "CHECK_case_1: " << iter << " " << pow(new_grain_quaternion.at(0),2) +  pow(new_grain_quaternion.at(1),2) +  pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2)<< "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;
        /// two maps from the list of Cases to (1) new "shifted" quaternions and (2) the corresponding grain rotated
        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_1.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_1)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation1 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;

        // case 2
        new_grain_q_quaternion = dq_upright(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));
// cout << "CHECK_case_2: " << iter << " "  << pow(new_grain_quaternion.at(0),2) +  pow(new_grain_quaternion.at(1),2) +  pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2)<< "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;

        /// two maps from the list of Cases to (1) new "shifted" quaternions and (2) the corresponding grain rotated
        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_2.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_2)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation2 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;
        // case 3
        new_grain_q_quaternion = dq_downleft(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

//cout << "CHECK_case_3: " << iter << " " << pow(new_grain_quaternion.at(0),2) +  pow(new_grain_quaternion.at(1),2) +  pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2)<< "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;

        /// two maps from the list of Cases to (1) new "shifted" quaternions and (2) the corresponding grain rotated
        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_3.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_3)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation3 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;
        // case 4
        //cout << "CHECK_case_4: " << iter << " " << pow(new_grain_quaternion.at(0),2) +  pow(new_grain_quaternion.at(1),2) +  pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2)<< "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;
        new_grain_q_quaternion = dq_downright(grain_q_quaternion, dq);
/// new_q0_coord remains the same (!!)
        new_q0_coord = q0_coord; // std::sqrt(1.0 - pow(new_grain_q_quaternion[0],2) - pow(new_grain_q_quaternion[1],2) - pow(new_grain_q_quaternion[2],2));
        new_grain_quaternion.at(0) = new_q0_coord; // angle
        // axis
        new_grain_quaternion.at(1) = std::sqrt(new_grain_q_quaternion.at(0)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(2) = std::sqrt(new_grain_q_quaternion.at(1)*(1.0 - pow(new_q0_coord,2)));
        new_grain_quaternion.at(3) = std::sqrt(new_grain_q_quaternion.at(2)*(1.0 - pow(new_q0_coord,2)));

//cout << "CHECK_case_4: " << iter << " " << pow(new_grain_quaternion.at(0),2) +  pow(new_grain_quaternion.at(1),2) +  pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2)<< "  " << distance(grain_quaternions_list.begin(),g_itr) << endl;

        /// two maps from the list of Cases to (1) new "shifted" quaternions and (2) the corresponding grain rotated
        cases_to_new_quaternions.insert(std::pair<unsigned int, std::vector<double>>(iter,new_grain_quaternion));
        cases_to_grains.insert(std::pair<unsigned int, unsigned int> (iter, distance(grain_quaternions_list.begin(),g_itr)));
// REPAIR    cout << distance(grain_quaternions_list.begin(),g_itr) << endl;

        for(unsigned int ngb : gb_set) // all the special GBs related with the rotated grain
            if(isHAGB(grain_full_quaternion, new_grain_quaternion, HAGBs_threshold))
                gb_special_set_4.push_back(ngb);

        /// map from the list of Cases to special Faces set
        cases_to_sfaces.insert(std::pair<unsigned int, std::vector<unsigned int>>(iter++,gb_special_set_4)); // each cases coresponds to the vector with a single element which is the same number of special face

//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;
//    cout << "disorientation3 " << Get_2grains_FCCdisorientation(grain_full_quaternion, new_grain_quaternion) << endl;

// REPAIR if (gb_special_set_1.size() != 0 || gb_special_set_2.size() != 0 || gb_special_set_3.size() != 0 || gb_special_set_4.size() != 0) cout << "q_check " << grain_q_quaternion.at(0) + grain_q_quaternion.at(1) + grain_q_quaternion.at(2) << " ngb " << gb_special_set_1.size() << "  " << gb_special_set_2.size() << "  " << gb_special_set_3.size() << "  " << gb_special_set_4.size() <<endl;
//    cout << "quaternions old " << pow(grain_full_quaternion.at(0),2) + pow(grain_full_quaternion.at(1),2) + pow(grain_full_quaternion.at(2),2) + pow(grain_full_quaternion.at(3),2) << "  " << grain_full_quaternion.at(0) << "  " << grain_full_quaternion.at(1) << "  " << grain_full_quaternion.at(2) << "  " << grain_full_quaternion.at(3) <<endl;
//    cout << "quaternions new " << pow(new_grain_quaternion.at(0),2) + pow(new_grain_quaternion.at(1),2) + pow(new_grain_quaternion.at(2),2) + pow(new_grain_quaternion.at(3),2) << "  " << new_grain_quaternion.at(0) << "  " << new_grain_quaternion.at(1) << "  " << new_grain_quaternion.at(2) << "  " << new_grain_quaternion.at(3) <<endl;

//        cout << "cases_random_list " << "place 4: " << distance(grain_quaternions_list.begin(),g_itr) << endl;

/// New State Vectors for Faces
        NewSpecialFaces = S_Vector;
        if (gb_special_set_1.size() != 0 || gb_special_set_2.size() != 0 || gb_special_set_3.size() != 0 || gb_special_set_4.size() != 0) {
// case 1
            NewSpecialFaces = S_Vector;
            for (unsigned int gb: gb_special_set_1) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    NewSpecialFaces.at(gb) = 1;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewSpecialFaces);
// case 2
            NewSpecialFaces = S_Vector;
            for (unsigned int gb: gb_special_set_2) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    NewSpecialFaces.at(gb) = 1;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewSpecialFaces);
// case 3
            NewSpecialFaces = S_Vector;
            for (unsigned int gb: gb_special_set_3) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    NewSpecialFaces.at(gb) = 1;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewSpecialFaces);
// case 4
            NewSpecialFaces = S_Vector;
            for (unsigned int gb: gb_special_set_4) // loop over all Edges
                if (S_Vector.at(gb) == 0) { // ORDINARY faces only
                    NewSpecialFaces.at(gb) = 1;
                } // end of if (S_Vector.at(gb) == 0)

            cryst_cases_list.push_back(NewSpecialFaces);

        } // end of if (gb_special_set_1.size() != 0 ||...
//        cout << "cases_random_list " << "place 5: " << distance(grain_quaternions_list.begin(),g_itr) << endl;

    } // end of for (auto  g_itr = grain_quaternions_list.begin(); g_itr != grain_quaternions_list.end(); ++g_itr)
//        cout << "cases_random_list " << "place 6: " << endl;

    return cryst_cases_list; // function output
} // END of Get_crystallographic_cases_random_list()

/// Configuration TJs entropy
double Get_TJsEntropy(vector<unsigned int> special_faces_seq) {
    double TJsEntropy = 0.0;

    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    vector<double> TJsTypes;

    SpMat FES(CellNumbs.at(1), CellNumbs.at(2));
    FES = SMatrixReader(paths.at(5 ), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces

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

std::vector<int> state_vector_by_sequence(std::vector<unsigned int> const &cell_sequence, int cell_type) { // cell_type: 0 -nodes, 1 - edges, 2 - faces, 3 - polyhedrons
std::vector<int> state_vector(CellNumbs.at(cell_type), 0); // including 2D case

for(auto cs : cell_sequence)
state_vector.at(cs) = 1;

return state_vector;
}

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
*/

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