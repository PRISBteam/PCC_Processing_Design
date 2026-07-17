///================================ PCC Subcomplex module ===================================================================================///
///=========================================================================================================================================///
///* Creates a set of subcomplexes of a PCC *///
///* -----------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2026 years as a module of PCC Processing Design code (CPD code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years *///
/// (https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1)                                                            *///
///==================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>

// External

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
#include "functions/subcomplex_cross_section.h"
#include "functions/subcomplex_clusters.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; //number of cells in a PCC defined globally
extern std::string source_path;
extern std::string output_dir;
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::ofstream subcomplex_logfile_stream;

#include "PCC_Subcomplex.h"
///* ========================================================= PCC SUBCOMPLEX FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Create a vector of PCC complexes with their special and induced labels taken from the read initial PCC. These subcomplexes can be either related to a specific PCC k-cell, or
 * covering a geometrical object in the initial tessellation, such as line, plane or volume. Currently, the subsection types are:
 * 'H' for the half-plane cut like a crack, 'N' for a k-order neighbouring grain set
 * @param configuration
 * @return std::vector<Subcomplex>
 */
std::vector<Subcomplex> PCC_Subcomplex(Config &configuration) {
    std::vector<Subcomplex> subcomplexes_vector; // function output

    cout << endl << "==============================================================================================================================================================" << endl;
    subcomplex_logfile_stream << endl << "==============================================================================================================================================================" << endl;

    // Reads simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen.
    // Change the source directory and simulation type from file ..\config\subcomplex.ini, if needed.
    config_reader_subcomplex(configuration);

    std::string subcomplex_mode = configuration.Get_subcomplex_mode(); // 'H' or 'N' :: This char define the subsection type: 'H' for the half-plane cut like a crack, 'N' for a k-order neighbouring grain set
    cout << endl << "Subcomplex\t" << subcomplex_mode << "\tmode in operation" << endl << "--------------------------------" << endl;
    subcomplex_logfile_stream << endl << "Subcomplex\t" << subcomplex_mode << "\tmode in operation" << endl << "--------------------------------" << endl;

    std::vector<double> plane_orientation(4); // orientation of a geometric plane in a tessellation for the 'P' subcomplex mode
    plane_orientation = configuration.Get_subcomplex_plane(); // 'H' mode only - four numbers {a,b,c,D} specified plane orientation read from the 'config/subcomplex.ini' file
    double cut_length = configuration.Get_subcomplex_cut_length(); // set initial cut length, before its reading from the corresponding ini-file
    unsigned int grain_neighbour_orders = configuration.Get_subcomplex_grain_neighbour_orders(); // the order of k-neighbours of a grain for the 'N' subcomplex generation mode
    bool is_log_file_output = configuration.Get_is_subcomplex_log_file(); // optional output to the subcomplex module log-file

    // TERMINOLOGY:  sub_polytope_set - all grains in the subcomplex, doubled_sub_faces_sequence - all faces in the subcomplex repeated twice (as boundaries of 2 grains each), internal_sub_faces_set - all faces common for two grains within the subcomplex in a set, sub_sfaces_sequence - special faces, sub_cfaces_sequence - induced (fractured, for instance) faces
    std::set<unsigned int> sub_polytope_set, sub_faces_set, internal_sub_faces_set;
    std::vector <unsigned int> doubled_sub_faces_sequence, sub_sfaces_sequence, sub_cfaces_sequence;
    //TODO: check the terminology -- cfaces or ffaces everywhere
    std::vector<tuple<double, double, double>> subcomplex_polytope_coordinates, subcomplex_face_coordinates, internal_faces_coordinates;

    Eigen::SparseMatrix<double> AGS = SMatrixReader(paths_to_PCC_matrices.at(3 + (PCC_dimension - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes (grains\3-cells)
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose()));  //  Full symmetric AGS matrix instead of triagonal
    Eigen::SparseMatrix<double> GFS = SMatrixReader(paths_to_PCC_matrices.at(6 + (PCC_dimension - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Volumes
    Eigen::SparseMatrix<double> FES = SMatrixReader(paths_to_PCC_matrices.at(5 + (PCC_dimension - 3)), (CellNumbs.at(1 + (PCC_dimension - 3))),(CellNumbs.at(2 + (PCC_dimension - 3)))); //all Edges-Faces
/// Vertex coordinates reader from file into triplet double vector
    polytope_coordinates_vector = Tuple3Reader(paths_to_PCC_matrices.at(9)); // grain seeds reader
    node_coordinates_vector = Tuple3Reader(paths_to_PCC_matrices.at(10)); // vertex seeds reader

/// All subcomplex grains (subcomplex_grain_sequence) for the plane cut
    if(subcomplex_mode == "H") {
        Subcomplex new_plane_subcomplex; // initialisation of a new subcomplex

        sub_polytope_set = PCC_Subcomplex_plane_cut_grains(plane_orientation);

        //TEST for (auto u : sub_polytope_set) cout << "sub_polytope_set_grains: " << u << endl; //    cout << "polytope_coordinates_vector.size(): " << polytope_coordinates_vector.size() << endl; // REPAIR   sub_polytope_set.clear(); for (unsigned int k = 0; k < CellNumbs.at(3); ++k) sub_polytope_set.insert(k);
        cout << "Subcomplex polytope set size\t=\t" << sub_polytope_set.size() << endl;
        subcomplex_logfile_stream << "Subcomplex polytope set size\t=\t" << sub_polytope_set.size() << endl;

/// Common grain coordinates for the PCC's 'internal' faces
///-------------------------------------------------
        if (sub_polytope_set.size() > 0) {
            for (auto subgc: sub_polytope_set)
                subcomplex_polytope_coordinates.push_back(polytope_coordinates_vector.at(subgc));
        }
        else {
            cout << "Caution! sub_polytope_set.size() = 0 in PCC_Subcomplex.cpp" << endl;
            subcomplex_logfile_stream << "Caution! sub_polytope_set.size() = 0 in PCC_Subcomplex.cpp" << endl;
        }

/// All subcomplex faces (doubled_sub_faces_sequence)
        doubled_sub_faces_sequence.clear();
        internal_sub_faces_set.clear();

        sub_faces_set = new_plane_subcomplex.Get_unique_subfaces_set(sub_polytope_set, GFS, doubled_sub_faces_sequence);

        cout << "Subcomplex unique faces set size:\t\t=\t" << sub_faces_set.size() << endl;
        subcomplex_logfile_stream << "Subcomplex unique faces set size:\t\t=\t" << sub_faces_set.size() << endl;


        internal_sub_faces_set = new_plane_subcomplex.Get_internal_subfaces_set(doubled_sub_faces_sequence);
        cout << "Internal Subcomplex unique faces sequence size:\t=\t" << internal_sub_faces_set.size() << endl << endl;
        subcomplex_logfile_stream << "Internal Subcomplex unique faces sequence size:\t=\t" << internal_sub_faces_set.size() << endl << endl;


//        cout << "Internal Subcomplex unique edges sequence size:\t=\t" << new_plane_subcomplex.Get_internal_subedges_set(FES).size() << endl << endl;

/// Setting all quantities to the subcomplex new_subPCC with id = 0
        new_plane_subcomplex.Set_internal_sub_faces_set(internal_sub_faces_set);
        new_plane_subcomplex.Set_sub_polytope_set(sub_polytope_set);
        new_plane_subcomplex.Set_sub_faces_set(sub_faces_set);
        new_plane_subcomplex.Set_sub_internal_face_coordinates(internal_faces_coordinates);
        new_plane_subcomplex.Set_sub_polytope_coordinates(subcomplex_polytope_coordinates);

        subcomplexes_vector.push_back(new_plane_subcomplex);

    } // end of the 'H' mode

    // -------------------------------------------------------- // --------------------------------------------------------------------- //
    // ========================================== *** anther computation mode *** ====================================================== //
    // -------------------------------------------------------- // --------------------------------------------------------------------- //

    else if(subcomplex_mode == "N") { // grain clusters M_k(p)
        bool subtype_reading = false;

        std::ifstream subcomplex_instream;
        std::string file_path_name = configuration.Get_source_dir() + "subcomplexes/grain_neighbours_"s + std::to_string(grain_neighbour_orders) + "_orders.txt"s;
        subcomplex_instream.open(file_path_name);

        std::vector<std::vector<int>> vector_of_sub_polytope_vectors;
        if (subcomplex_instream.is_open()) {
            subtype_reading = true;
            vector_of_sub_polytope_vectors = IntListReader(subcomplex_instream);
         }

        std::vector<std::set<unsigned int>> vector_of_sub_polytope_sets;
        std::set<unsigned int> sub_polytope_sets_probe;
        for( std::vector<int> probe_vec : vector_of_sub_polytope_vectors) {
            for (int val: probe_vec) {
                sub_polytope_sets_probe.insert(val);
            }
            vector_of_sub_polytope_sets.push_back(sub_polytope_sets_probe);
        }

        unsigned int grain_number; // 'gb' is for Grain Boundary or interfaces
        if (PCC_dimension == 3)
            grain_number = CellNumbs.at(3);
        else if (PCC_dimension == 2)
            grain_number = CellNumbs.at(2);

        for (unsigned int grain_id = 0; grain_id < grain_number; ++grain_id) {

            if(subtype_reading && vector_of_sub_polytope_sets.size() > 0)
                sub_polytope_set = vector_of_sub_polytope_sets.at(grain_id);
            else
                sub_polytope_set = PCC_Subcomplex_k_order_grain_neighbours(grain_id, grain_neighbour_orders);
        //REPAIR for (auto u : sub_polytope_set) cout << "sub_polytope_set_grains: " << u << endl; //    cout << "polytope_coordinates_vector.size(): " << polytope_coordinates_vector.size() << endl; // REPAIR   sub_polytope_set.clear(); for (unsigned int k = 0; k < CellNumbs.at(3); ++k) sub_polytope_set.insert(k);

/// Common grain coordinates for 'internal' faces
///-------------------------------------------------
        if (sub_polytope_set.size() > 0) {
            for (auto subgc: sub_polytope_set)
                subcomplex_polytope_coordinates.push_back(polytope_coordinates_vector.at(subgc));
        }
        else {
            cout << "Caution! sub_polytope_set.size() = 0 in DCC_Subcomplex.h" << endl;
            subcomplex_logfile_stream << "Caution! sub_polytope_set.size() = 0 in DCC_Subcomplex.h" << endl;
        }

/// All subcomplex faces (doubled_sub_faces_sequence)
        doubled_sub_faces_sequence.clear();
        internal_sub_faces_set.clear();

        if (sub_polytope_set.size() > 0) {
            for (auto grain_id : sub_polytope_set) { // for each grain in a subcomplex
                for (unsigned int l = 0; l < CellNumbs.at(2); ++l) { // for each face
                    if (GFS.coeff(l, grain_id) != 0)
                        doubled_sub_faces_sequence.push_back(l);
                } // end for (auto grain_id : sub_polytope_set)
            } // end of for (unsigned int l = 0; l < CellNumbs.at(2); l++)

/// Full subcomplex faces set
            for (auto unique_sub_faces: doubled_sub_faces_sequence)
                sub_faces_set.insert(unique_sub_faces); // set automatically remove all repetitions

///            cout << "Subcomplex faces set size:\t\t=\t" << sub_faces_set.size() << endl;
///            subcomplex_logfile_stream << "Subcomplex faces set size:\t\t=\t" << sub_faces_set.size() << endl;

            for (auto face_id: doubled_sub_faces_sequence) {
// REPAIR        cout << count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) << endl;
                if (count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) > 1) {
                    internal_sub_faces_set.insert(face_id);
                }
            }
            cout << grain_id << "\tSubcomplex polytope set size =\t" << sub_polytope_set.size() << "\tInternal Subcomplex faces sequence size =\t" << internal_sub_faces_set.size() << endl;
            subcomplex_logfile_stream << grain_id << "\tSubcomplex polytope set size =\t" << sub_polytope_set.size() << "\tInternal Subcomplex faces sequence size =\t" << internal_sub_faces_set.size() << endl;

        }  //if (sub_polytope_set.size() > 0) {

        Subcomplex new_subcomplex; // new subcomplex with its ID
/// Setting all quantities to the subcomplex new_subPCC with id = 0
        new_subcomplex.Set_internal_sub_faces_set(internal_sub_faces_set);
        new_subcomplex.Set_sub_polytope_set(sub_polytope_set);
        new_subcomplex.Set_sub_faces_set(sub_faces_set);
        new_subcomplex.Set_sub_internal_face_coordinates(internal_faces_coordinates);
        new_subcomplex.Set_sub_polytope_coordinates(subcomplex_polytope_coordinates);

        subcomplexes_vector.push_back(new_subcomplex);
        } // end for (unsigned int grain_id = 0; grain_id < grain_number; ++grain_id)

        cout << "-----------------" << endl << endl;
        subcomplex_logfile_stream << "-----------------" << endl << endl;
        for (auto sc : subcomplexes_vector) {
            cout << endl;
            subcomplex_logfile_stream << endl;

            auto pol_set = sc.Get_sub_polytope_set();
            for (auto ps : sc.Get_sub_polytope_set()) {
                cout << ps << " ";
                subcomplex_logfile_stream << ps << " ";
            }
        }

    } // end if(subcomplex_mode == "N")

    return subcomplexes_vector;

} /// END of the Subcomplex PCC_Subcomplex()
