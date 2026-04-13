///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains top-down (based on the (k+1)-cell types) and bottom-up (based on the (k-1)-cell types) indexing of        **/
/**  k-cells, k={0,1,2,3} using the 'special' and 'induced' cell types assigned in the PCC Processing module.                     **/
///==============================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random> // Require C++ 11 and above
// #include <execution> // Require C++ 17 and above

// external libraries
#include "../../../external/Eigen/Core"
#include "../../../external/Eigen/SparseCore"

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

extern std::vector<unsigned int> CellNumbs;
extern std::vector<std::string> paths_to_PCC_matrices;
extern int PCC_dimension;

#include "processing_indexing.h"
///------------------------------------------------------------------
/*!
 * @details
 * @param k_cell_type
 * @param Configuration_State
 * @return unsigned int Vector 'TopDownTypes(CellNumbs.at(k_cell_type))' of the (k+1) labels induced directly by the number of incident (k+1)-cell types
 */
std::vector<unsigned int> TopDown_cell_indexing(int k_cell_type, std::vector<std::vector<unsigned int>> &Configuration_State) {
// Output vector of the (k+1) labels induced directly by the number of incident (k+1)-cell types
    std::vector<unsigned int> TopDownTypes(CellNumbs.at(k_cell_type), 0); // CellNumbs.at(k_cell_type) is the number of k-cells

// Obtaining (k+1)-cells (coloumns) - k-cells (rows) Incidence matrix using the file paths.at(5 + (dim - 3))
    SpMat FES;
    if(k_cell_type != 3) {
        FES = SMatrixReader(paths_to_PCC_matrices.at((k_cell_type + 4) + (PCC_dimension - 3)), CellNumbs.at(k_cell_type + (PCC_dimension - 3)), CellNumbs.at((k_cell_type + 1) + (PCC_dimension - 3))); // k-(k+1)) sparse incidence matrix; (k_cell_type + 4) is the corresponding paths to B_k incidence matrices
    }
    else{
        cout << "ERROR: TopDownTypes() function cannot be applied to 3-cells (!)" << endl;
        exit(1);
    }
// Creating sequence from the corresponding State Vector (configuration)
    std::vector<unsigned int> special_kp1_local_sequence;
    for (auto it = Configuration_State[k_cell_type+1].begin(); it != Configuration_State[k_cell_type+1].end(); ++it)
        if(*it > 0) {
            special_kp1_local_sequence.push_back(distance(Configuration_State[k_cell_type+1].begin(), it)); // add new element to the s_cells_sequence
        } // end if(it)

    for (auto kp1: special_kp1_local_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(k_cell_type + (PCC_dimension - 3)); ++k) // loop over all Edges
            if (FES.coeff(k, kp1) != 0) TopDownTypes.at(k)++;

    return TopDownTypes;
}

// RESTRICTION: up to 3 types of grain phases (!)
std::vector<unsigned int> TopDown_cell_indexing(int k_cell_type, std::vector<unsigned int> &higher_order_cells_state_vector) {
// Output vector of the (k+1) labels induced directly by the number of incident (k+1)-cell types
    std::vector<unsigned int> TopDownTypes(CellNumbs.at(k_cell_type), 0); // CellNumbs.at(k_cell_type) is the number of k-cells

// Obtaining (k+1)-cells (coloumns) - k-cells (rows) Incidence matrix using the file paths.at(5 + (dim - 3))
    SpMat FES;
    if(k_cell_type != 3) {
        FES = SMatrixReader(paths_to_PCC_matrices.at((k_cell_type + 4) + (PCC_dimension - 3)), CellNumbs.at(k_cell_type + (PCC_dimension - 3)), CellNumbs.at((k_cell_type + 1) + (PCC_dimension - 3))); // k-(k+1)) sparse incidence matrix; (k_cell_type + 4) is the corresponding paths to B_k incidence matrices
    }
    else{
        cout << "ERROR: TopDownTypes() function cannot be applied to 3-cells (!)" << endl;
        exit(1);
    }
// Creating sequence from the corresponding State Vector (configuration)
    std::vector<unsigned int> special_kp1_local_sequence, special_kp2_local_sequence, special_kp3_local_sequence;
    for (auto it = higher_order_cells_state_vector.begin(); it != higher_order_cells_state_vector.end(); ++it)
        if(*it == 1) {
            special_kp1_local_sequence.push_back(distance(higher_order_cells_state_vector.begin(), it)); // add new element to the s_cells_sequence
        } // end if(it)
        else if(*it == 2) {
            special_kp2_local_sequence.push_back(distance(higher_order_cells_state_vector.begin(), it)); // add new element to the s_cells_sequence
        }
        else if(*it == 3) {
            special_kp3_local_sequence.push_back(distance(higher_order_cells_state_vector.begin(), it)); // add new element to the s_cells_sequence
        }

    for (auto kp1: special_kp1_local_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(k_cell_type + (PCC_dimension - 3)); ++k) // loop over all Edges
            if (FES.coeff(k, kp1) != 0) TopDownTypes.at(k)++;

    for (auto kp2: special_kp2_local_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(k_cell_type + (PCC_dimension - 3)); ++k) // loop over all Edges
            if (FES.coeff(k, kp2) != 0) TopDownTypes.at(k) += 3;

    for (auto kp2: special_kp3_local_sequence) // loop over all Special Faces
        for(int k = 0; k < CellNumbs.at(k_cell_type + (PCC_dimension - 3)); ++k) // loop over all Edges
            if (FES.coeff(k, kp2) != 0) TopDownTypes.at(k) += 5;

    return TopDownTypes;
}