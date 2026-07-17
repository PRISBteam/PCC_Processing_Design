///================================ A part of the PCC Subcomplex module =============================================================///
///=================================================================================================================================///
/** The library contains functions providing k-Cells IDs for a given set of subcomplexes made as clusters of PCC cells               **/
///================================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <random> // Require C++ 11 and above
// #include <execution> // Require C++ 17 and above

// external libraries
#include "../../../external/Eigen-5.0/Core"
#include "../../../external/Eigen-5.0/SparseCore"

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::ofstream subcomplex_logfile_stream;

#include "subcomplex_clusters.h"

/// ======# 1 #================= Subcomplex PCC_Subcomplex_k_order_grain_neighbours() function ==============================================================///
std::set<unsigned int> PCC_Subcomplex_k_order_grain_neighbours(unsigned int grain_id, int k_neighbours_order){
    std::set<unsigned int> k_order_grain_neighbours_set; // function output

    unsigned int grain_number; // 'gb' is for Grain Boundary or interfaces
    if (PCC_dimension == 3)
        grain_number = CellNumbs.at(3);
    else if (PCC_dimension == 2)
        grain_number = CellNumbs.at(2);

    //    Eigen::SparseMatrix<double> adjacency grain matrix
    Eigen::SparseMatrix<double> AGS = SMatrixReader(paths_to_PCC_matrices.at(3 + (PCC_dimension - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose())); // Full symmetric AGS matrix instead of triagonal

std::set<unsigned int> higher_order_grain_neighbours_set, new_grain_neighbours_set;
// k = 0 - grain itself is its own 0-eighbour
    higher_order_grain_neighbours_set.insert(grain_id);
// neighbours k > 0
for (int k_order = 1; k_order <= k_neighbours_order; ++k_order) {    // looking for neighbours for each grain in the PCC
    new_grain_neighbours_set.clear();

        for (unsigned int gnn: higher_order_grain_neighbours_set) {
            for (int i = 1; i < grain_number; ++i) {
                if (AGS.coeff(gnn, i) != 0)
                    new_grain_neighbours_set.insert(i);
            } // i
        } // gnn
    higher_order_grain_neighbours_set = new_grain_neighbours_set;
} // for (int k_order = 1; k_order < k_neighbours_order; ++k_order)

    k_order_grain_neighbours_set = higher_order_grain_neighbours_set;

// ==========================
//    half_plane_cut.Set_sub_polytope_set(half_sub_grain_set); // all grains before cut
//    half_plane_cut.Set_internal_sub_faces_set(half_internal_faces_set); // common faces
//    std::vector<unsigned int> half_sfaces_seq = SetToVector(half_sfaces_set);
//    half_plane_cut.Set_sub_sfaces_sequence(half_sfaces_seq); //special faces
//    std::vector<std::tuple<double, double, double>> half_plane_sfaces_coord = kSequence_barycentre_coordinates(2, half_sfaces_seq);
//    cout << " half_sfaces_set SIZE " << half_sfaces_set.size() << endl;
//    subcomplex_logfile_stream << " half_sfaces_set SIZE " << half_sfaces_set.size() << endl;
//    cout << " half_sfaces_seq SIZE " << half_sfaces_seq.size() << endl;
//    subcomplex_logfile_stream << " half_sfaces_seq SIZE " << half_sfaces_seq.size() << endl;
//    cout << " half_plane_sfaces_coord SIZE " << half_plane_sfaces_coord.size() << endl;
//    subcomplex_logfile_stream << " half_plane_sfaces_coord SIZE " << half_plane_sfaces_coord.size() << endl;

//    half_plane_cut.Set_sub_sfaces_coord(half_plane_sfaces_coord);

    // ==========================


    return k_order_grain_neighbours_set;
}