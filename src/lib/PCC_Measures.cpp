///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains functions which calculate various measures such as configuration entropy.                                 **/
///==============================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// external libraries (must be added to the STL directory!)
#include <Eigen/Core>
#include <Eigen/SparseCore>

// local libraries
#include "PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace
using namespace Eigen; // standard namespace

// global variables
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> PCCpaths;
extern int dim;

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

#include "PCC_Measures.h"
///------------------------------------------------------------------

std::vector<int> Edge_types_byFaces(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &special_face_sequence, std::vector<double> &j_fractions, std::vector<double> &d_fractions)
{
// Output of the model
    std::vector<int> TJsTypes(CellNumbs.at(1 + (dim - 3)), 0); // CellNumbs.at(1) is the number of Edges

// Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(PCCpaths.at(5 + (dim - 3)), CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))); // Edges-Faces sparse incidence matrix

    for (auto f : special_face_sequence) // loop over all Special Faces
        for(int e = 0; e < CellNumbs.at(1 + (dim - 3)); ++e) // loop over all Edges
            if (FES.coeff(e, f) != 0) TJsTypes.at(e)++;

    unsigned int J0 = 0, J1 = 0, J2 = 0, J3 = 0;
    double Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
/// amounts of TJs of different types
    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
    J0 = CellNumbs.at(1) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
//    Jall = (double) CellNumbs.at(1 + (dim - 3)); // amount of Edges in DCC
    Jall = (double) TJsTypes.size();

// Conversion from numbers to fractions
    double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
    j0 = (double) J0 / Jall;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    if (j0 > 0) l2j0 = log2(j0);
    j1 = (double) J1 / Jall;
    if (j1 > 0) l2j1 = log2(j1);
    j2 = (double) J2 / Jall;
    if (j2 > 0) l2j2 = log2(j2);
    j3 = (double) J3 / Jall;
    if (j3 > 0) l2j3 = log2(j3);

    j_fractions.at(0) = j0; j_fractions.at(1) = j1; j_fractions.at(2) = j2; j_fractions.at(3) = j3;

    for (int k = 0; k < 3; ++k) {
        if (j_fractions.at(0) != 1.0)
            d_fractions.at(k) = j_fractions.at(k + 1) / (1.0 - j_fractions.at(0));
        else d_fractions.at(k) = 0.0;
    }

    return TJsTypes;
}

double Face_edge_index(std::vector<unsigned int> &special_face_sequence, Eigen::SparseMatrix<double> const& FES, double norm_const) {
    double face_edge_index = 0.0;
    std::vector<double> Face_inclusion_index(CellNumbs.at(2),0), j_fractions(4,0), d_fractions(3,0);
    std::vector<double> TJsTypes; std::vector<int> TJsTypes_int;
    TJsTypes_int = Edge_types_byFaces(CellNumbs, special_face_sequence, j_fractions, d_fractions);

    for (int var : TJsTypes_int)
        TJsTypes.push_back(var);

    for (unsigned int fn = 0; fn < CellNumbs.at(2); ++fn) {
        std::vector<double> j_types_neigh_fractions = GBIndex(fn, FES, TJsTypes);         //std::vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, vector<double> const& TJsTypes) {
        /// Inclusion index computation
        Face_inclusion_index.at(fn) = j_types_neigh_fractions.at(1) + 2.0 * j_types_neigh_fractions.at(2) + 3.0 * j_types_neigh_fractions.at(3);
//REPAIR        if(Face_inclusion_index.at(fn) != 0)  cout << "\tFace_inclusion_index.at(fn)\t" << Face_inclusion_index.at(fn) << endl;
    } // end of for ( fn < CellNumbs.at(cell_type))

    face_edge_index = 0.0;
    std::vector<double> BL_indices;
    unsigned int bli_max= *std::max_element(Face_inclusion_index.begin(),Face_inclusion_index.end());
    for(unsigned int bl_numb = 1; bl_numb <= bli_max; ++bl_numb)
        BL_indices.push_back(std::count(Face_inclusion_index.begin(), Face_inclusion_index.end(), bl_numb));

    for(unsigned int bl_numb = 0; bl_numb < BL_indices.size(); ++bl_numb)
        face_edge_index += BL_indices.at(bl_numb) * (double) bl_numb/ (norm_const*CellNumbs.at(2));

    return face_edge_index;
} // end of Face_edge_index

double Node_edge_index(std::vector<unsigned int> &special_face_sequence, Eigen::SparseMatrix<double> const& ENS, double norm_const) { // based on Edges vector
    double node_edge_index = 0.0;
    std::vector<double> Node_inclusion_index(CellNumbs.at(0),0), j_fractions(4,0), d_fractions(3,0);
    std::vector<double> TJsTypes; std::vector<int> TJsTypes_int;
    TJsTypes_int = Edge_types_byFaces(CellNumbs, special_face_sequence, j_fractions, d_fractions);

    for (int var : TJsTypes_int)
        TJsTypes.push_back(var);

    for (unsigned int fn = 0; fn < CellNumbs.at(0); ++fn) {
        std::vector<double> j_types_neigh_fractions = NodeIndex(fn, ENS, TJsTypes);
        /// Inclusion index computation
        Node_inclusion_index.at(fn) = j_types_neigh_fractions.at(1) + 2.0 * j_types_neigh_fractions.at(2) + 3.0 * j_types_neigh_fractions.at(3);
    } // end of for ( fn < CellNumbs.at(cell_type))

    node_edge_index = 0.0;
    std::vector<double> QL_indices;

    unsigned int bli_max= *std::max_element(Node_inclusion_index.begin(),Node_inclusion_index.end());
    for(unsigned int bl_numb = 1; bl_numb < bli_max; ++bl_numb)
        QL_indices.push_back(std::count(Node_inclusion_index.begin(),Node_inclusion_index.end(), bl_numb));

//    for(unsigned int bl_numb = 1; bl_numb <= bli_max; ++bl_numb)
    for(unsigned int bl_numb = 0; bl_numb < QL_indices.size(); ++bl_numb)
        node_edge_index += QL_indices.at(bl_numb) * std::exp((double) bl_numb)/ (norm_const*CellNumbs.at(0));

    return node_edge_index;
} // end of Node_edge_index

/// * Function calculates the vector<int> "EdgeTypes" of types Edges in the PCC using its FES incidence matrix and special faces sequence (special_faces_sequence) * ///
std::vector<double> j_fractions_vector(std::vector<double> const &TJsTypes){ // based on Edges vector
std::vector<double> j_fractions_vector(4); // Function output: TJs fractions

unsigned int J0 = 0, J1 = 0, J2 = 0, J3 = 0;
double j0 = 0, j1 = 0, j2 = 0, j3 = 0, Jall = 0;

/// amounts of TJs of different types
J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
J0 = CellNumbs.at(1) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
Jall = (double) CellNumbs.at(1); // amount of Edges in DCC

// Conversion from numbers to fractions
j_fractions_vector.at(0) = (double) J0 / Jall;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
j_fractions_vector.at(1) = (double) J1 / Jall;
j_fractions_vector.at(2) = (double) J2 / Jall;
j_fractions_vector.at(3) = (double) J3 / Jall;

return j_fractions_vector;
} // end of j_fractions()

std::vector<double> d_fractions_vector(vector<double> const &TJsTypes){ // based on Edges vector
std::vector<double> j_fractions(4), d_fractions_vector(4); // Function output: TJs fractions
/// j_fractions calculation
j_fractions = j_fractions_vector(TJsTypes);

// Conversion from numbers to fractions
/// (!) d[0] - > 1 ; d[1] - > 2 ; d[2] - > 3
d_fractions_vector.at(0) = j_fractions.at(1) / (1.0 - j_fractions.at(0));
d_fractions_vector.at(1) = j_fractions.at(2) / (1.0 - j_fractions.at(0));
d_fractions_vector.at(2) = j_fractions.at(3) / (1.0 - j_fractions.at(0));

return d_fractions_vector;
} // end of d_fractions()

double Configuration_Entropy(vector<double> &TJsTypes){ // based on Edges vector
double Configuration_Face_Entropy = 0.0; // Function output: Configuration Entropy related with Faces

unsigned int J0 = 0, J1 = 0, J2 = 0, J3 = 0;
double j0 = 0, j1 = 0, j2 = 0, j3 = 0, Jall = 0;
/// amounts of TJs of different types
J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
J0 = CellNumbs.at(1) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
Jall = (double) CellNumbs.at(1); // amount of Edges in DCC

// Conversion from numbers to fractions
double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
j0 = (double) J0 / Jall;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j0 > 0) l2j0 = log2(j0);
j1 = (double) J1 / Jall;
if (j1 > 0) l2j1 = log2(j1);
j2 = (double) J2 / Jall;
if (j2 > 0) l2j2 = log2(j2);
j3 = (double) J3 / Jall;
if (j3 > 0) l2j3 = log2(j3);

return Configuration_Face_Entropy = -(j0 * l2j0 + j1 * l2j1 + j2 * l2j2 + j3 * l2j3);
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;
} // end of Configuration_Entropy()

double Configuration_Entropy(std::vector<double> const &j_fractions){ // based on Edges fraction vector
double Configuration_Face_Entropy = 0.0; // Function output: Configuration Entropy related with Faces

double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j_fractions.at(0) > 0)
l2j0 = log2(j_fractions.at(0));
if (j_fractions.at(1) > 0)
l2j1 = log2(j_fractions.at(1));
if (j_fractions.at(2) > 0)
l2j2 = log2(j_fractions.at(2));
if (j_fractions.at(3) > 0)
l2j3 = log2(j_fractions.at(3));

Configuration_Face_Entropy = -(j_fractions.at(0) * l2j0 + j_fractions.at(1) * l2j1 + j_fractions.at(2) * l2j2 + j_fractions.at(3) * l2j3);
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;

return Configuration_Face_Entropy;
} // end of Configuration_Entropy()

std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions){ // based on Edges fraction vector
//double Configuration_Face_Entropy = 0.0; // Function output: Configuration Entropy related with Faces
double Face_Entropy_Median = 0, Face_Entropy_Skrew = 0;

double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j_fractions.at(0) > 0)
l2j0 = log2(j_fractions.at(0));
if (j_fractions.at(1) > 0)
l2j1 = log2(j_fractions.at(1));
if (j_fractions.at(2) > 0)
l2j2 = log2(j_fractions.at(2));
if (j_fractions.at(3) > 0)
l2j3 = log2(j_fractions.at(3));

double j0s = j_fractions.at(0), j1s = j_fractions.at(1), j2s = j_fractions.at(2), j3s = j_fractions.at(3);
std::vector<double> j_types_fractions = {j_fractions.at(0), j_fractions.at(1), j_fractions.at(2), j_fractions.at(3)};

/// Median part in the entropy decomposition
if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(j_fractions.at(0) * j_fractions.at(1) * j_fractions.at(2) * j_fractions.at(3));
} else Face_Entropy_Median = 0.0;

/// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
for (int j = 0; j < j_types_fractions.size(); ++j)
for (int i = 0; i < j; ++i)
if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
Face_Entropy_Skrew +=
-(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
log2(j_types_fractions[i] / j_types_fractions[j]);
} else Face_Entropy_Skrew += 0.0;
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;

return std::make_tuple(Face_Entropy_Median, Face_Entropy_Skrew);
} // end of Configuration_Entropy()

/// Configuration Entropy Differences

double Configuration_Entropy_change(vector<int> const &TJsTypes, vector<int> const &TJsTypesNew){ // based on Edges vector
double Configuration_Face_Entropy = 0, Configuration_Face_Entropy_new = 0, Configuration_Face_Entropy_increase = 0.0; // Function output: Configuration Entropy related with Faces

unsigned int J0 = 0, J1 = 0, J2 = 0, J3 = 0;
double j0 = 0, j1 = 0, j2 = 0, j3 = 0;
/// amounts of TJs of different types
J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
J0 = CellNumbs.at(1) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
double Jall = (double) CellNumbs.at(1); // amount of Edges in DCC

// Conversion from numbers to fractions
double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
j0 = (double) J0 / Jall;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j0 > 0) l2j0 = log2(j0);
j1 = (double) J1 / Jall;
if (j1 > 0) l2j1 = log2(j1);
j2 = (double) J2 / Jall;
if (j2 > 0) l2j2 = log2(j2);
j3 = (double) J3 / Jall;
if (j3 > 0) l2j3 = log2(j3);

// initial configuration entropy
Configuration_Face_Entropy = -(j0 * l2j0 + j1 * l2j1 + j2 * l2j2 + j3 * l2j3);

unsigned int J0n = 0, J1n = 0, J2n = 0, J3n = 0;
double j0n = 0, j1n = 0, j2n = 0, j3n = 0;
/// amounts of TJs of different types
J1n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 1); // containing 1 incident special face
J2n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 2); // containing 2 incident special face
J3n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 3); // containing 3 incident special face
J0n = CellNumbs.at(1) - J1n - J2n - J3n; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))

// Conversion from numbers to fractions
double l2j0n = 0.0, l2j1n = 0.0, l2j2n = 0.0, l2j3n = 0.0;
j0n = (double) J0n / Jall;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j0n > 0) l2j0n = log2(j0n);
j1n = (double) J1n / Jall;
if (j1n > 0) l2j1n = log2(j1n);
j2n = (double) J2n / Jall;
if (j2n > 0) l2j2n = log2(j2n);
j3n = (double) J3n / Jall;
if (j3n > 0) l2j3n = log2(j3n);

// updated configuration entropy
Configuration_Face_Entropy_new = -(j0n * l2j0n + j1n * l2j1n + j2n * l2j2n + j3n * l2j3n);

return Configuration_Face_Entropy_increase = Configuration_Face_Entropy_new - Configuration_Face_Entropy;
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;
} // end of Configuration_Entropy_change()

double Configuration_Entropy_change(std::vector<double> const &j_fractions, vector<int> const &TJsTypesNew){ // based on Edges vector
double Configuration_Face_Entropy = 0, Configuration_Face_Entropy_new = 0, Configuration_Face_Entropy_increase = 0.0; // Function output: Configuration Entropy related with Faces

double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j_fractions.at(0) > 0)
l2j0 = log2(j_fractions.at(0));
if (j_fractions.at(1) > 0)
l2j1 = log2(j_fractions.at(1));
if (j_fractions.at(2) > 0)
l2j2 = log2(j_fractions.at(2));
if (j_fractions.at(3) > 0)
l2j3 = log2(j_fractions.at(3));

// initial configuration entropy
Configuration_Face_Entropy = -(j_fractions.at(0) * l2j0 + j_fractions.at(1) * l2j1 + j_fractions.at(2) * l2j2 + j_fractions.at(3) * l2j3);

unsigned int J0n = 0, J1n = 0, J2n = 0, J3n = 0;
double j0n = 0, j1n = 0, j2n = 0, j3n = 0;
/// amounts of TJs of different types
J1n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 1); // containing 1 incident special face
J2n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 2); // containing 2 incident special face
J3n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 3); // containing 3 incident special face
J0n = CellNumbs.at(1) - J1n - J2n - J3n; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))

// Conversion from numbers to fractions
double l2j0n = 0.0, l2j1n = 0.0, l2j2n = 0.0, l2j3n = 0.0;
double Jall = (double) CellNumbs.at(1); // amount of Edges in DCC

// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
j0n = (double) J0n / Jall;
if (j0n > 0) l2j0n = log2(j0n);
j1n = (double) J1n / Jall;
if (j1n > 0) l2j1n = log2(j1n);
j2n = (double) J2n / Jall;
if (j2n > 0) l2j2n = log2(j2n);
j3n = (double) J3n / Jall;
if (j3n > 0) l2j3n = log2(j3n);

// updated configuration entropy
Configuration_Face_Entropy_new = -(j0n * l2j0n + j1n * l2j1n + j2n * l2j2n + j3n * l2j3n);

Configuration_Face_Entropy_increase = Configuration_Face_Entropy_new - Configuration_Face_Entropy;

return Configuration_Face_Entropy_increase;
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;
} // end of Configuration_Entropy()

double Configuration_DevEntropy_change(std::vector<double> const &j_fractions, vector<int> const &TJsTypesNew){ // based on Edges vector
double Configuration_Edge_DevEntropy_increase = 0.0; // Function output: Configuration Entropy related with Faces

double l2j0 = 0.0, l2j1 = 0.0, l2j2 = 0.0, l2j3 = 0.0;
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
if (j_fractions.at(0) > 0)
l2j0 = log2(j_fractions.at(0));
if (j_fractions.at(1) > 0)
l2j1 = log2(j_fractions.at(1));
if (j_fractions.at(2) > 0)
l2j2 = log2(j_fractions.at(2));
if (j_fractions.at(3) > 0)
l2j3 = log2(j_fractions.at(3));

double Edge_Entropy_Skew = 0;
std::vector<double> j_types_fractions = {j_fractions.at(0), j_fractions.at(1), j_fractions.at(2), j_fractions.at(3)};
/// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
for (int j = 0; j < j_types_fractions.size(); ++j)
for (int i = 0; i < j; ++i)
if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
Edge_Entropy_Skew +=
-(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
log2(j_types_fractions[i] / j_types_fractions[j]);
} else Edge_Entropy_Skew += 0.0;

unsigned int J0n = 0, J1n = 0, J2n = 0, J3n = 0;
double j0n = 0, j1n = 0, j2n = 0, j3n = 0;
/// amounts of TJs of different types
J1n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 1); // containing 1 incident special face
J2n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 2); // containing 2 incident special face
J3n = std::count(TJsTypesNew.begin(), TJsTypesNew.end(), 3); // containing 3 incident special face
J0n = CellNumbs.at(1) - J1n - J2n - J3n; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))

// Conversion from numbers to fractions
double l2j0n = 0.0, l2j1n = 0.0, l2j2n = 0.0, l2j3n = 0.0;
double Jall = (double) CellNumbs.at(1); // amount of Edges in DCC

// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
j0n = (double) J0n / Jall;
j1n = (double) J1n / Jall;
j2n = (double) J2n / Jall;
j3n = (double) J3n / Jall;

double Edge_Entropy_Skew_new = 0;
std::vector<double> j_types_fractions_new = {j0n, j1n, j2n, j3n};
/// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
for (int j = 0; j < j_types_fractions_new.size(); ++j)
for (int i = 0; i < j; ++i)
if (j_types_fractions_new[i]!=0 && j_types_fractions_new[j]!=0) {
Edge_Entropy_Skew_new +=
-(1.0 / j_types_fractions_new.size()) * (j_types_fractions_new[i] - j_types_fractions_new[j]) *
log2(j_types_fractions_new[i] / j_types_fractions_new[j]);
} else Edge_Entropy_Skew_new += 0.0;

Configuration_Edge_DevEntropy_increase = Edge_Entropy_Skew - Edge_Entropy_Skew_new;

return Configuration_Edge_DevEntropy_increase;
// REPAIR cout << "Configuration_Face_Entropy: " << Configuration_Face_Entropy << endl;
} // end of Configuration_DeviatoricEntropy()
