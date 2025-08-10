///================================ PCC Design module ===================================================================================///
///=========================================================================================================================================///
///* Finds an optimal special and induced types of the State Vectors. *///
///* -----------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design code (CPD code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years *///
/// (https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1)                                                            *///
///==================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

/// Attached user-defined C++ libraries:
// External
#include "../../../src/lib/external/Eigen/SparseCore"

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
//#include "functions/..h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; //number of cells in a PCC defined globally
extern ofstream Out_logfile_stream;
extern std::string source_path;
extern std::string output_dir;
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally

#include "PCC_Design.h"
///* ========================================================= PCC SUBCOMPLEX FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Output sequences of their evolution towards the optimum defined by a 'goal' function.
 * @param configuration
 * @return std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> PCC_Design(Config &configuration){
/// Main output of the module
    std::vector<std::vector<int>> design_list_of_vectors;

    return design_list_of_vectors;
}
