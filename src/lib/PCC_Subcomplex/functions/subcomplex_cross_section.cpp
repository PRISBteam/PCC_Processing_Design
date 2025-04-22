///================================ A part of the PCC Subcomplex module =============================================================///
///=================================================================================================================================///
/** The library contains functions providing k-Cells IDs for a given subcomplex cut rules               **/
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
#include <Eigen/Core>
#include <Eigen/SparseCore>

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> PCCpaths; // PCCpaths to PCC files
extern int dim; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally

#include "subcomplex_cross_section.h"
///* ========================================================= PCC Plane Cut functions ======================================================= *///
///* ========================================================================================================================================= *///

/// ======# 1 #================= std::vector<unsigned int> PCC_Plane_cut () function ==============================================================///
std::set<unsigned int> PCC_Plane_cut_grains (std::vector<double> &plane_orientation) {
/// The plane parameters: a_coeff*X + b_coeff*Y + c_coeff*Z = D
    double a_coeff = plane_orientation.at(0);
    double b_coeff = plane_orientation.at(1);
    double c_coeff = plane_orientation.at(2);
    double D_coeff = plane_orientation.at(3);

    std::set<unsigned int> planecut_grains;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file PCCpaths.at(5 + (dim - 3))
    Eigen::SparseMatrix<double> AGS = SMatrixReader(PCCpaths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose())); // Full symmetric AGS matrix instead of triagonal

    Eigen::SparseMatrix<double> ENS = SMatrixReader(PCCpaths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    Eigen::SparseMatrix<double> FES = SMatrixReader(PCCpaths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    Eigen::SparseMatrix<double> GFS = SMatrixReader(PCCpaths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains

//REPAIR    for (auto  itr = node_coordinates_vector.begin(); itr != node_coordinates_vector.end(); ++itr) cout << "\tx\t" << get<0>(*itr) << "\ty\t" << get<1>(*itr) << "\tz\t" << get<2>(*itr) << endl;  exit(0);

/// Creation the list of all polytopes
    std:vector<Polytope> grains_list; // creates vector of Polytopes of size CellNumbs.at(3) and fill it out with '0's
    for(unsigned int m = 0; m < CellNumbs.at(3); ++m) { // each 3-cell
        Polytope new_polytope(m);
        grains_list.push_back(new_polytope);
    }

/// #pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); ++m) { // each 3-cell
        grains_list.at(m).Set_node_ids(GFS, FES, ENS);
//REPAIR  cout << node_coordinates_vector.size() <<  "\tgrain ID\t" << grains_list.at(m).grain_id << endl;         cout << grains_list.at(m).Get_node_coordinates().size() <<  "\tgrain ID\t" << grains_list.at(m).grain_id << endl;
        grains_list.at(m).Set_node_coordinates(node_coordinates_vector);
//REPAIR cout << grains_list.size() << " grain id:  " << grains_list.at(m).grain_id << endl; //for (auto did : grains_list.at(m).Get_node_coordinates()) cout << get<0>(did) << "  "<< get<1>(did) << "  "<< get<2>(did) << "  ";  cout << endl;
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) - Grains
    // REPAIR    for(unsigned int m = 0; m < CellNumbs.at(3); m++) // each Grain cout << " grain_id = " << grains_list.at(m).grain_id << endl;

/// For each grain minmax_coord vector grain_coordinate_extremums of two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}
/// #pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        std::vector<tuple<double, double, double>> grain_coordinate_extremums = grains_list.at(m).Get_minmax_node_coordinates(); // get<0>(grain_coordinate_extremums.at(0)) -> MIN X coordinate, get<1>(grain_coordinate_extremums.at(0)) -> MIN Y coordinate, get<0>(grain_coordinate_extremums.at(1)) -> MAX X coordinate, etc..
//REPAIR cout <<" gain # : " << m << " " << grains_list.at(m).grain_id << endl; //REPAIR cout << "Xmin " << get<0>(grain_coordinate_extremums.at(0)) << " Ymin " <<get<1>(grain_coordinate_extremums.at(0)) << " Zmin " << get<2>(grain_coordinate_extremums.at(0)) << endl; //REPAIR cout << "Xmax " << get<0>(grain_coordinate_extremums.at(1)) << " Ymax " <<get<1>(grain_coordinate_extremums.at(1)) << " Zmax " << get<2>(grain_coordinate_extremums.at(1)) << endl;
        ///MinMax condition: / a_coeff*X + b_coeff*Y + c_coeff*Z = D
        if(a_coeff*get<0>(grain_coordinate_extremums.at(0)) + b_coeff*get<1>(grain_coordinate_extremums.at(0)) + c_coeff*get<2>(grain_coordinate_extremums.at(0)) < D_coeff  &&
           a_coeff*get<0>(grain_coordinate_extremums.at(1)) + b_coeff*get<1>(grain_coordinate_extremums.at(1)) + c_coeff*get<2>(grain_coordinate_extremums.at(1)) > D_coeff ) // simultaneously: z_min < D < z_max
            planecut_grains.insert(m);
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
    //REPAIR for (auto gid : planecut_grains) cout << gid << endl;

return planecut_grains;
} /// END of the std::vector<unsigned int> PCC_Plane_cut () function

/// Overloaded function
std::set<unsigned int> PCC_Plane_cut_grains (double a_coeff, double b_coeff, double c_coeff, double D_coeff) {
/// The plane parameters: a_coeff*X + b_coeff*Y + c_coeff*Z = D

    std::set<unsigned int> planecut_grains;
    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file PCCpaths.at(5 + (dim - 3))
    Eigen::SparseMatrix<double> ENS = SMatrixReader(PCCpaths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    //SpMat AES = SMatrixReader(PCCpaths.at(1 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
    ///  Full symmetric AES matrix instead of triagonal
    //AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose()));

    Eigen::SparseMatrix<double> FES = SMatrixReader(PCCpaths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    //SpMat AFS = SMatrixReader(PCCpaths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    ///  Full symmetric AFS matrix instead of triagonal
    //AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose()));
    Eigen::SparseMatrix<double> GFS = SMatrixReader(PCCpaths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    Eigen::SparseMatrix<double> AGS = SMatrixReader(PCCpaths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose()));

    /// Grains
    vector<Polytope> grains_list(CellNumbs.at(3),0); // vector of grains (class grain3D)
/// #pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        Polytope new_grain = Polytope(m); //        cout << new_grain.grain_id << endl;
        grains_list.at(m) = new_grain; //        cout <<  "grains_list: " << grains_list.size() << endl; //        cout << GFS.cols() << "  " << CellNumbs.at(3) << endl;
        grains_list.at(m).Set_node_ids(GFS, FES, ENS); //        cout << grains_list.at(m).grain_id << endl; //        cout << "m = " << m << " grain list: " << grains_list.at(m).Get_node_ids(m).size() << endl;
        grains_list.at(m).Set_node_coordinates(node_coordinates_vector); //                cout << grains_list.at(m).grain_id << endl;
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) - Grains

    /// For each grain minmax_coord vector grain_coordinate_extremums of two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}
/// #pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        std::vector<tuple<double, double, double>> grain_coordinate_extremums = grains_list.at(m).Get_minmax_node_coordinates(); // get<0>(grain_coordinate_extremums.at(0)) -> MIN X coordinate, get<1>(grain_coordinate_extremums.at(0)) -> MIN Y coordinate, get<0>(grain_coordinate_extremums.at(1)) -> MAX X coordinate, etc..
//REPAIR cout <<" gain # : " << m << " " << grains_list.at(m).grain_id << endl; //REPAIR cout << "Xmin " << get<0>(grain_coordinate_extremums.at(0)) << " Ymin " <<get<1>(grain_coordinate_extremums.at(0)) << " Zmin " << get<2>(grain_coordinate_extremums.at(0)) << endl; //REPAIR cout << "Xmax " << get<0>(grain_coordinate_extremums.at(1)) << " Ymax " <<get<1>(grain_coordinate_extremums.at(1)) << " Zmax " << get<2>(grain_coordinate_extremums.at(1)) << endl;
        ///MinMax condition: / a_coeff*X + b_coeff*Y + c_coeff*Z = D
        if(a_coeff*get<0>(grain_coordinate_extremums.at(0)) + b_coeff*get<1>(grain_coordinate_extremums.at(0)) + c_coeff*get<2>(grain_coordinate_extremums.at(0)) < D_coeff  &&
           a_coeff*get<0>(grain_coordinate_extremums.at(1)) + b_coeff*get<1>(grain_coordinate_extremums.at(1)) + c_coeff*get<2>(grain_coordinate_extremums.at(1)) > D_coeff ) // simultaneously: z_min < D < z_max
            planecut_grains.insert(m);
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain //a,b,c,d
    //REPAIR    for (auto gid : planecut_grains)         cout << gid << endl;

    return planecut_grains;
} /// END of the std::vector<unsigned int> PCC_Plane_cut () function


/// ======# 2 #================= Subcomplex Get_half_plane() function ==============================================================///
Subcomplex Get_half_plane(Subcomplex &plane_subcomplex, double crack_length){
    Subcomplex half_plane_cut;

    int direction = 0;
    if (plane_subcomplex.a_n != 0) { direction = 0; } // x
    else if (plane_subcomplex.b_n != 0) { direction = 1; } // y
    else direction = 2; // z

    //        vector<tuple<double, double, double>> grain_coordinates = polytope_coordinates_vector;
//REPAIR        cout << " polytope_coordinates_vector " << grain_coordinates.size() << endl;
    //cout << "Xmin " << get<0>(minmax_tuple.at(0)) << " Ymin " <<get<1>(minmax_tuple.at(0)) << " Zmin " << get<2>(minmax_tuple.at(0)) << endl;
//    std::vector<tuple<double, double, double>> all_face_coordinates = face_coordinates_vector;
//REPAIR        cout << " face_coordinates_vector " << face_coordinates.size() << endl;

//        for (auto  itr = grain_coordinates.begin(); itr != grain_coordinates.end(); ++itr)
//            if (std::find(plane_subcomplex.Get_sfaces_sequence().begin(), plane_subcomplex.Get_sfaces_sequence().end(), distance(plane_subcomplex.Get_sfaces_sequence().begin(),itr)) != plane_subcomplex.Get_sfaces_sequence().end() && get(direction, *itr) < crack_length)
//                    half_sub_grains_set.push_back(distance(grain_coordinates.begin(),itr));
//REPAIR cout << "half_sub.Get_grains_sequence(0) " << half_sub.Get_grains_sequence(0).size() << endl;

    std::set<unsigned int> half_sub_grain_set, sub_common_face_set;

    for (auto half_grains : plane_subcomplex.Get_sub_polytope_set()) {
        if (get_i(direction,polytope_coordinates_vector.at(half_grains)) < crack_length) {
//REPAIR                cout << " half_grains " << half_grains << " half_grain_coordinates " << get<0>(polytope_coordinates_vector.at(half_grains)) << " crack_length " << crack_length << endl;
            half_sub_grain_set.insert(half_grains);
        }
    }

    std::set<unsigned int> half_internal_faces_set, half_sfaces_set;

    half_internal_faces_set.clear();
        for (auto  itr2 = face_coordinates_vector.begin(); itr2 != face_coordinates_vector.end(); ++itr2)
            if (std::find(sub_common_face_set.begin(), sub_common_face_set.end(), distance(face_coordinates_vector.begin(),itr2)) != sub_common_face_set.end() && get_i(direction, *itr2) < crack_length)
                half_internal_faces_set.insert(distance(face_coordinates_vector.begin(),itr2));

    cout << "\tcrack_length\t\t" << crack_length << "\t" << endl;

    half_sfaces_set.clear();
    std::vector<unsigned int> local_sfaces_sequence = plane_subcomplex.Get_sub_sfaces_sequence();
    cout << "\tlocal_sfaces_sequence\t\t" << local_sfaces_sequence.size() << "\t" << endl;
    cout << "\tface_coordinates_vector\t\t" << face_coordinates_vector.size() << "\t" << endl;

    for (auto  itr3 = face_coordinates_vector.begin(); itr3 != face_coordinates_vector.end(); ++itr3) {
/// REPAIR        cout << "\tdistance(face_coordinates_vector.begin(), itr3)\t\t" << distance(face_coordinates_vector.begin(), itr3) << "\t" << "crack_length\t" << crack_length << "\t" << endl;
        if (std::find(local_sfaces_sequence.begin(), local_sfaces_sequence.end(), distance(face_coordinates_vector.begin(), itr3)) != local_sfaces_sequence.end() && get_i(direction, *itr3) < crack_length) {
                half_sfaces_set.insert(distance(face_coordinates_vector.begin(), itr3));
            }
        }

//    exit(0);

    half_plane_cut.Set_sub_polytope_set(half_sub_grain_set); // all grains before cut
    half_plane_cut.Set_internal_faces_set(half_internal_faces_set); // common faces
    std::vector<unsigned int> half_sfaces_seq = SetToVector(half_sfaces_set);
    half_plane_cut.Set_sub_sfaces_sequence(half_sfaces_seq); //special faces
    std::vector<std::tuple<double, double, double>> half_plane_sfaces_coord = kSequence_barycentre_coordinates(2, half_sfaces_seq);
    cout << " half_sfaces_set SIZE " << half_sfaces_set.size() << endl;
    cout << " half_sfaces_seq SIZE " << half_sfaces_seq.size() << endl;
    cout << " half_plane_sfaces_coord SIZE " << half_plane_sfaces_coord.size() << endl;

    half_plane_cut.Set_sub_sfaces_coord(half_plane_sfaces_coord);
    //half_plane_cut.Set_common_faces_coordinates(common_faces_coordinates);
    // half_plane_cut.Set_sub_grain_coordinates(subcomplex_grain_coordinates);
    //half_plane_cut.Set_cfaces_sequence(c_sub_faces_sequence); //cracked (induced) faces
    half_plane_cut.sub_length = crack_length;

    return half_plane_cut;

} // end of Get_half_plane() function


/// Subcomplex Get_half_plane(Subcomplex new_sub, double crack_length, std::vector<unsigned int> const &half_sub_sfaces_sequence);
