///================================ DCC Processing module =============================================================///
///=======================================================================================================================///
/** The function in this library generate quasi-random process of changes in the elements of the pre-constructed          ///
*   discrete sell complex (DCC) with the whole set of incidence and adjacency martices                                   **/
///=======================================================================================================================///

/// Standard (STL) C++ libraries:
///------------------------------
/// Attached user defined C++ libraries:
///-------------------------------------
#include "Processing_Functions.h"
///-------------------------------------

using namespace std;
using namespace Eigen;
//using namespace Spectra;

/// Function defenitions
std::vector<unsigned int> VectorReader(char* FilePath);

int DCC_Processing3D(std::vector<unsigned int> &State_Vector, std::vector<unsigned int>  &special_faces_sequence, char stype, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs, std::vector<char*> const paths) {
// State_Vector in the form : [Element number] - > [Type]
// CellNumbs :: vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
// Maximal fraction (max_sFaces_fraction) for simulation loop max_sFaces_fraction = [0,1]

//    double ordinary_faces_fraction = 1.0, special_faces_fraction = 0.0;
//    std::vector<int> SpecialCells(CellNumbs->at(2), 0), OrdinaryCells(CellNumbs->at(2),1); // New vectors initialised with 0 for SpecialCells and 1 for OrdinaryCells
//    vector<map<unsigned int, unsigned int>> Vector_SCellMaps(number_of_types);
//    vector<unsigned int>::iterator vit;// Special iterator for such vector

    if (stype == 'R') { //  Random generation case
        Processing_Random(State_Vector, special_faces_sequence, max_sFaces_fraction, number_of_types, CellNumbs);
    } ///End of 'R' type simulations

    else if (stype == 'S') { // Maximum entropy production
         Processing_maxEntropy(State_Vector, special_faces_sequence, max_sFaces_fraction, number_of_types, CellNumbs, paths);
    } ///End of 'S' type simulations

    else if (stype == 'D') { // DDRX recrystalisation process
        Processing_DDRX(State_Vector, special_faces_sequence, max_sFaces_fraction, paths, number_of_types, CellNumbs);
    } ///End of 'D' type simulations

    else if (stype == 'I') {

    }

    else if (stype == 'E') {

    }
    else { cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl; return 888;}

// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces /// Sparse adjacency matrix from the corresponding triplet list of special faces
//        SpMat Id(numerator, numerator), SAM_FacesGraph(numerator, numerator), Face_Laplacian(numerator, numerator), Sym_Face_Laplacian(numerator, numerator), RW_Face_Laplacian(numerator, numerator);
//        SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

//*  /// Key point: creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
//*            for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
//*                if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
//*                    SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

/*            long unsigned int sit = SpecialCellMap[New2CellNumb]; // Just a useful variable for the number of newly converted Face (= numerator--)
            /// Creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
            for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
                if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
                    SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}
*/
/*
              ///=============================================================================================================================================////
            ///=================================== Characterisation module =============================>>
            if( OCellAmount % output_step == 0 && ordinary_faces_fraction > 0.05) { // Periodicity of characterisation output
                unsigned int SAM_size = Structure_Characterisation(numerator, CellNumbs, SFaces_Triplet_list, SpecialCellMap, special_faces_sequence, FES, odir, special_faces_fraction, configuration);

                cout << "Fraction of special faces:" << 1.0 - ordinary_faces_fraction << ";\t Size of S-Face Adjacency matrix \t" << SAM_size << endl;
            } // End of analysis and output iterator ( IF: iterator % X == 0 )
 */

    return 0;
} /// The end of HAGBsProbability3D()

/// 2D version (DCC with max 2-Cells without 3-Cells and higher dimensions) of the code
void HAGBsProbability2D(char* AN, char* AE, char* AF, char* MEN, char* MFE, char* CNumbs, char* config, char* odir, char stype) {
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// odir - source path
// stype - simulation type ('R', 'S', 'I',....)
} /// The end of HAGBsProbability2D()


/// ================================== Related functions ==================================///

/// Heap
// (1)    cout << std::numeric_limits<long unsigned int>::max() << endl;

// (2) ofstream NONStructGlCoordOut; NONStructGlCoordOut.open("C:VAnRes\\(MFE21)GNC_NONStructured.txt", ios::trunc);

//(3) ARCHIVE testing output  //   cout << "The number of columns in the created matrices:\t" << "\tANS --\t" << ANS.cols() << "\tAES --\t" << AES.cols()<< "\tAFS --\t" << AFS.cols() << "\tAGS --\t" << AGS.cols() << "\tENS --\t" << ENS.cols() << "\tFES --\t" << FES.cols() << "\tGFS --\t" << GFS.cols() << endl;
