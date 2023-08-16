///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains functions for generation of new identifiers for m-Cells, m={0,1,2} based on the already assigned types    **/
/**  of the k-Cells with k > m, in the PCC Processing module. It created an "imposed" structure of m-Cells in a PCC.              **/
///==============================================================================================================================///

//#include "Processing_Imposed_functions.h"

vector<int> Edge_types_byFaces(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &special_face_sequence, std::vector<double> &j_fractions, std::vector<double> &d_fractions)
{
// Output of the model
    vector<int> TJsTypes(CellNumbs.at(1 + (dim - 3)), 0); // CellNumbs.at(1) is the number of Edges

// Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))); // Edges-Faces sparse incidence matrix

    for (auto f: special_face_sequence) // loop over all Special Faces
        for(int e = 0; e < CellNumbs.at(1 + (dim - 3)); ++e) // loop over all Edges
            if (FES.coeff(e, f) != 0) TJsTypes.at(e)++;

    unsigned int J0 = 0, J1 = 0, J2 = 0, J3 = 0;
    double Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
/// amounts of TJs of different types
    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1); // containing 1 incident special face
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2); // containing 2 incident special face
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3); // containing 3 incident special face
    J0 = CellNumbs.at(1 + (dim0 - 3)) - J1 - J2 - J3; // containing no incident special face (the total amount of TJs = total amount of Edges in DCC = CellNumbs.at(1))
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