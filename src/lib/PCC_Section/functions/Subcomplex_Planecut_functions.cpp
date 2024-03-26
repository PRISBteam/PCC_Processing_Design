///================================ A part of the PCC Subcomplex module =============================================================///
///=================================================================================================================================///
/** The library contains functions providing k-Cells IDs for a given subcomplex cut rules               **/
///================================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
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
extern std::vector<unsigned int> CellNumbs; //number of cells in a PCC defined globally
extern std::vector<char*> PCCpaths; //PCCpaths to PCC files
extern int dim; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // coordinate vectors defined globally

#include "Subcomplex_Planecut_functions.h"
///* ========================================================= PCC Plane Cut functions ======================================================= *///
///* ========================================================================================================================================= *///

/// ======# 1 #================= std::vector<unsigned int> PCC_Plane_cut () function ==============================================================///
std::vector<unsigned int> PCC_Plane_cut (double a_coeff, double b_coeff, double c_coeff, double D_coeff) {
/// The plane parameters: a_coeff*X + b_coeff*Y + c_coeff*Z = D

    vector<unsigned int> planecut_grains;

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
        Polytope new_grain = Polytope(m);
//        cout << new_grain.grain_id << endl;
        grains_list.at(m) = new_grain;
//        cout <<  "grains_list: " << grains_list.size() << endl;
//        cout << GFS.cols() << "  " << CellNumbs.at(3) << endl;
        grains_list.at(m).Set_node_ids(m, GFS, FES, ENS);
//        cout << grains_list.at(m).grain_id << endl;
//        cout << "m = " << m << " grain list: " << grains_list.at(m).Get_node_ids(m).size() << endl;
        grains_list.at(m).Set_node_coordinates(m);
//                cout << grains_list.at(m).grain_id << endl;
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) - Grains

    /// For each grain minmax_coord vector grain_coordinate_extremums of two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}
/// #pragma omp parallel for // parallel execution by OpenMP
    for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
        vector<tuple<double, double, double>> grain_coordinate_extremums = grains_list.at(m).Get_minmax_node_coordinates(m); // get<0>(grain_coordinate_extremums.at(0)) -> MIN X coordinate, get<1>(grain_coordinate_extremums.at(0)) -> MIN Y coordinate, get<0>(grain_coordinate_extremums.at(1)) -> MAX X coordinate, etc..
//REPAIR cout <<" gain # : " << m << " " << grains_list.at(m).grain_id << endl;
//REPAIR cout << "Xmin " << get<0>(grain_coordinate_extremums.at(0)) << " Ymin " <<get<1>(grain_coordinate_extremums.at(0)) << " Zmin " << get<2>(grain_coordinate_extremums.at(0)) << endl;
//REPAIR cout << "Xmax " << get<0>(grain_coordinate_extremums.at(1)) << " Ymax " <<get<1>(grain_coordinate_extremums.at(1)) << " Zmax " << get<2>(grain_coordinate_extremums.at(1)) << endl;
        ///MinMax condition: / a_coeff*X + b_coeff*Y + c_coeff*Z = D
        if(a_coeff*get<0>(grain_coordinate_extremums.at(0)) + b_coeff*get<1>(grain_coordinate_extremums.at(0)) + c_coeff*get<2>(grain_coordinate_extremums.at(0)) < D_coeff  &&
           a_coeff*get<0>(grain_coordinate_extremums.at(1)) + b_coeff*get<1>(grain_coordinate_extremums.at(1)) + c_coeff*get<2>(grain_coordinate_extremums.at(1)) > D_coeff ) // simultaneously: z_min < D < z_max
            planecut_grains.push_back(m);
    } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
//a,b,c,d

//REPAIR    for (auto gid : planecut_grains)         cout << gid << endl;

    return planecut_grains;
} /// END of the std::vector<unsigned int> PCC_Plane_cut () function

/// ======# 2 #================= Subcomplex Get_half_plane() function ==============================================================///
Subcomplex Get_half_plane(Subcomplex new_sub, double crack_length, std::vector<unsigned int> const &half_sub_sfaces_sequence){
//        vector<tuple<double, double, double>> grain_coordinates = grain_coordinates_vector;
//REPAIR        cout << " grain_coordinates_vector " << grain_coordinates.size() << endl;
    //cout << "Xmin " << get<0>(minmax_tuple.at(0)) << " Ymin " <<get<1>(minmax_tuple.at(0)) << " Zmin " << get<2>(minmax_tuple.at(0)) << endl;
    std::vector<tuple<double, double, double>> face_coordinates = face_coordinates_vector;
//REPAIR        cout << " face_coordinates_vector " << face_coordinates.size() << endl;

    Subcomplex half_plane_cut(1);
    std::vector<unsigned int> half_sub_grains_sequence, half_common_faces_sequence;
//        for (auto  itr = grain_coordinates.begin(); itr != grain_coordinates.end(); ++itr)
//            if (std::find(sub_grains_sequence.begin(), sub_grains_sequence.end(), distance(grain_coordinates.begin(),itr)) != sub_grains_sequence.end()
//                    && get<0>(*itr) < crack_length) half_sub_grains_sequence.push_back(distance(grain_coordinates.begin(),itr));
//REPAIR cout << "half_sub.Get_grains_sequence(0) " << half_sub.Get_grains_sequence(0).size() << endl;
    for (auto half_grains : new_sub.Get_grains_sequence(0)) {
        if (get<1>(grain_coordinates_vector.at(half_grains)) < crack_length) {
//REPAIR                cout << " half_grains " << half_grains << " half_grain_coordinates " << get<0>(grain_coordinates_vector.at(half_grains)) << " crack_length " << crack_length << endl;
            half_sub_grains_sequence.push_back(half_grains); //// get<0> ->> only along X now!!!
        }
    }
///TEMPORARILY!!!
//        for (auto  itr2 = face_coordinates.begin(); itr2 != face_coordinates.end(); ++itr2)
//            if (std::find(sub_common_face_sequence.begin(), sub_common_face_sequence.end(), distance(face_coordinates.begin(),itr2)) != sub_common_face_sequence.end()
//                    && get<0>(*itr2) < crack_length) half_common_faces_sequence.push_back(distance(face_coordinates.begin(),itr2));
//        for (auto  itr3 = face_coordinates.begin(); itr3 != face_coordinates.end(); ++itr3)
//            if (std::find(sub_sfaces_sequence.begin(), sub_sfaces_sequence.end(), distance(face_coordinates.begin(),itr3)) != sub_sfaces_sequence.end()
//                && get<0>(*itr3) < crack_length) half_sub_sfaces_sequence.push_back(distance(face_coordinates.begin(),itr3));

    half_plane_cut.Set_grains_sequence(half_sub_grains_sequence); // all grains before cut
    ///       half_plane_cut.Set_common_faces_sequence(half_common_faces_sequence); // common faces
    half_plane_cut.Set_sfaces_sequence(half_sub_sfaces_sequence); //special faces
    //half_plane_cut.Set_common_faces_coordinates(common_faces_coordinates);
    //half_plane_cut.Set_sub_grain_coordinates(subcomplex_grain_coordinates);
    //half_plane_cut.Set_cfaces_sequence(c_sub_faces_sequence); //cracked (induced) faces
    half_plane_cut.sub_length = crack_length;

    return half_plane_cut;

} // end of Get_half_plane() function


/// Subcomplex Get_half_plane(Subcomplex new_sub, double crack_length, std::vector<unsigned int> const &half_sub_sfaces_sequence);
