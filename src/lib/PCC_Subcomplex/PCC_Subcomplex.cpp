///================================ PCC Subcomplex module ===================================================================================///
///=========================================================================================================================================///
///* Creates a set of subcomplexes of a PCC     *///
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

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
#include "functions/subcomplex_cross_section.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; //number of cells in a PCC defined globally
extern ofstream Out_logfile_stream;
extern string source_path, output_dir;
extern std::vector<std::string> PCCpaths; // PCCpaths to PCC files
extern int dim; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally

#include "PCC_Subcomplex.h"
///* ========================================================= PCC SUBCOMPLEX FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Create a vector of PCC complexes with their special and induced labels taken from the initial PCC
 * @param configuration
 * @return std::vector<Subcomplex>
 */
std::vector<Subcomplex> PCC_Subcomplex(Config &configuration) { // sub_polytope_set - all grains in the subcomplex, doubled_sub_faces_sequence - all faces in the subcomplex, internal_doubled_sub_faces_sequence - all faces common for two grains in the subcomplex, sub_sfaces_sequence - special faces, sub_cfaces_sequence - induced (fractured, for instance) faces
std::vector<Subcomplex> several_cuts; // function output

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
/// The source directory and simulation type from file ..\config\subcomplex.ini
    std::string S_type; /// 'P', 'H' or 'N' :: This char define the subsection type: 'P' for the whole Plane cut, 'H' for the half-plane cut like a crack, 'N' for a k-order neighbouring grain set
    double cut_length = 0; // initial before reading from the corresponding ini-file
    std::vector<double> plane_orientation(4); /// for 'P' and 'H' modes only - four numbers {a,b,c,D} specified plane orientation read from the 'config/subcomplex.ini' file
    unsigned int grain_neighbour_orders = 0;

    /// Reading of the configuration from the '../config/subcomplex.ini' file
    config_reader_subcomplex(source_path, S_type, plane_orientation, cut_length, grain_neighbour_orders, Out_logfile_stream); // void function

    std::set<unsigned int> sub_polytope_set, internal_sub_faces_set, sub_faces_set;
    std::vector <unsigned int> doubled_sub_faces_sequence, sub_sfaces_sequence, sub_cfaces_sequence;
    std::vector<tuple<double, double, double>> subcomplex_polytope_coordinates, subcomplex_face_coordinates, internal_faces_coordinates;

    Eigen::SparseMatrix<double> AGS = SMatrixReader(PCCpaths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose()));  //  Full symmetric AGS matrix instead of triagonal
    Eigen::SparseMatrix<double> GFS = SMatrixReader(PCCpaths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Volumes

/// Vertex coordinates reader from file into triplet double vector
    polytope_coordinates_vector = Tuple3Reader(PCCpaths.at(9)); // grain seeds reader
    node_coordinates_vector = Tuple3Reader(PCCpaths.at(10)); // vertex seeds reader

/// All subcomplex grains (subcomplex_grain_sequence) for the plane cut
    sub_polytope_set = PCC_Plane_cut_grains(plane_orientation);
//REPAIR for (auto u : sub_polytope_set) cout << "sub_polytope_set_grains: " << u << endl; //    cout << "polytope_coordinates_vector.size(): " << polytope_coordinates_vector.size() << endl; // REPAIR   sub_polytope_set.clear(); for (unsigned int k = 0; k < CellNumbs.at(3); ++k) sub_polytope_set.insert(k);
    cout << "Subcomplex polytope set size\t=\t" << sub_polytope_set.size() << endl;

/// Common grain coordinates for 'internal' faces
///-------------------------------------------------
    if (sub_polytope_set.size() > 0) {
        for (auto subgc: sub_polytope_set)
            subcomplex_polytope_coordinates.push_back(polytope_coordinates_vector.at(subgc));
    }
    else cout << "Caution! sub_polytope_set.size() = 0 in DCC_Subcomplex.h" << endl;

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
  //  } // end of if(sub_polytope_set.size() > 0)

//    if (sub_polytope_set.size() > 0) {
//        for (auto grain_id: sub_polytope_set) {
//            Polytope polytope(grain_id);
//            polytope.Set_faces_list(GFS);

//            for (unsigned int l = 0; l < CellNumbs.at(2); l++) {
//                if (GFS.coeff(l, grain_id) == 1) {
//                    doubled_sub_faces_sequence.push_back(l);
//                        } // end of if (GFS.coeff(l, grain_id) == 1)
//                    } // end of for (unsigned int l = 0; l < CellNumbs.at(2); l++)

//                } // end of if (sub_polytope_set.size() > 0)
//            } // end of for (auto grain_id : sub_polytope_set)
//        } // end of for (auto grain_id2 : sub_polytope_set)
// REPAIR        cout << "doubled sub faces sequence size:\t=\t" << doubled_sub_faces_sequence.size() << endl;

            /// Full subcomplex faces set
            for (auto unique_sub_faces: doubled_sub_faces_sequence)
                sub_faces_set.insert(unique_sub_faces); // set automatically remove all repetitions
//        } // end of for (auto grain_id: sub_polytope_set)
        cout << "Subcomplex faces set size:\t\t=\t" << sub_faces_set.size() << endl;

        for (auto face_id: doubled_sub_faces_sequence) {
// REPAIR        cout << count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) << endl;
            if (count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) > 1) {
                internal_sub_faces_set.insert(face_id);
            }
        }
        cout << "Internal Subcomplex faces sequence size:\t=\t" << internal_sub_faces_set.size() << endl;
    }  //if (sub_polytope_set.size() > 0) {

// Common subcomplex faces (internal_doubled_sub_faces_sequence)
/*
if (doubled_sub_faces_sequence.size() > 0) {
    bool is_common = 0;
    unsigned int face_counter = 0;

    for (auto face_id: doubled_sub_faces_sequence) {
        is_common = 0;
        cout << "faces_sequence: " << count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) << endl;

        if (count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) > 1) {
            internal_doubled_sub_faces_sequence.push_back(face_id);
        }

        face_counter = face_id;
        for (auto face_id2: doubled_sub_faces_sequence) {
            if(face_id2 == face_counter) ++face_counter;
            if(face_counter > 2) is_common = 1;
        }
        if(is_common) cout << "Is common face here !" << endl;
    }
}
 */
/// Common face coordinates
/*
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber)
        if(std::find(internal_doubled_sub_faces_sequence.begin(), internal_doubled_sub_faces_sequence.end(), fnumber) != internal_doubled_sub_faces_sequence.end())
            internal_faces_coordinates.push_back(find_aGBseed(fnumber, PCCpaths, CellNumbs, polytope_coordinates_vector)); // tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);

//    for(auto cfs : internal_doubled_sub_faces_sequence) cout << "internal_doubled_sub_faces_sequence : " << cfs << endl;
    cout << "common faces coordinates size: " << internal_faces_coordinates.size() << endl;
*/

Subcomplex new_cut; // new subcomplex with its ID
/// Setting all quantities to the subcomplex new_subPCC with id = 0
    new_cut.Set_sub_polytope_set(sub_polytope_set);
    new_cut.Set_sub_faces_set(sub_faces_set);
    new_cut.Set_internal_face_coordinates(internal_faces_coordinates);
    new_cut.Set_sub_polytope_coordinates(subcomplex_polytope_coordinates);

    several_cuts.push_back(new_cut);
    return several_cuts;

} /// END of the Subcomplex PCC_Subcomplex(Subcomplex &new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> &doubled_sub_faces_sequence, std::vector<unsigned int> const &c_faces_sequence, double a_coeff = 0.0, double b_coeff = 0.0, double c_coeff = 1.0, double D_coeff = 0.6) {




std::vector<Subcomplex> PCC_Subcomplex(Config &configuration, std::vector <unsigned int> &all_sfaces_sequence){
    std::vector<Subcomplex> several_cuts; // function output

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
/// The source directory and simulation type from file ..\config\subcomplex.ini
    std::string S_type; /// 'P', 'H' or 'N' :: This char define the subsection type: 'P' for the whole Plane cut, 'H' for the half-plane cut like a crack, 'N' for a k-order neighbouring grain set
    double cut_length = 0; // initial before reading from the corresponding ini-file
    std::vector<double> plane_orientation(4); /// for 'P' and 'H' modes only - four numbers {a,b,c,D} specified plane orientation read from the 'config/subcomplex.ini' file
    unsigned int grain_neighbour_orders = 0;

    /// Reading of the configuration from the '../config/subcomplex.ini' file
    config_reader_subcomplex(source_path, S_type, plane_orientation, cut_length, grain_neighbour_orders, Out_logfile_stream); // void function

    std::set<unsigned int> sub_polytope_set, internal_sub_faces_set, sub_faces_set;
    std::vector <unsigned int> doubled_sub_faces_sequence, sub_sfaces_sequence, sub_cfaces_sequence;
    std::vector<tuple<double, double, double>> subcomplex_polytope_coordinates, subcomplex_face_coordinates, internal_faces_coordinates;

    Eigen::SparseMatrix<double> AGS = SMatrixReader(PCCpaths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    AGS = 0.5 * (AGS + Eigen::SparseMatrix<double>(AGS.transpose()));  //  Full symmetric AGS matrix instead of triagonal
    Eigen::SparseMatrix<double> GFS = SMatrixReader(PCCpaths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Volumes

/// Vertex coordinates reader from file into triplet double vector
    polytope_coordinates_vector = Tuple3Reader(PCCpaths.at(9)); // grain seeds reader
    node_coordinates_vector = Tuple3Reader(PCCpaths.at(10)); // vertex seeds reader

/// All subcomplex grains (subcomplex_grain_sequence) for the plane cut
    sub_polytope_set = PCC_Plane_cut_grains(plane_orientation);
//REPAIR for (auto u : sub_polytope_set) cout << "sub_polytope_set_grains: " << u << endl; //    cout << "polytope_coordinates_vector.size(): " << polytope_coordinates_vector.size() << endl; // REPAIR   sub_polytope_set.clear(); for (unsigned int k = 0; k < CellNumbs.at(3); ++k) sub_polytope_set.insert(k);
    cout << "Subcomplex polytope set size\t=\t" << sub_polytope_set.size() << endl;

/// Common grain coordinates for 'internal' faces
///-------------------------------------------------
    if (sub_polytope_set.size() > 0) {
        for (auto subgc: sub_polytope_set)
            subcomplex_polytope_coordinates.push_back(polytope_coordinates_vector.at(subgc));
    }
    else cout << "Caution! sub_polytope_set.size() = 0 in DCC_Subcomplex.h" << endl;

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
        //  } // end of if(sub_polytope_set.size() > 0)

//    if (sub_polytope_set.size() > 0) {
//        for (auto grain_id: sub_polytope_set) {
//            Polytope polytope(grain_id);
//            polytope.Set_faces_list(GFS);

//            for (unsigned int l = 0; l < CellNumbs.at(2); l++) {
//                if (GFS.coeff(l, grain_id) == 1) {
//                    doubled_sub_faces_sequence.push_back(l);
//                        } // end of if (GFS.coeff(l, grain_id) == 1)
//                    } // end of for (unsigned int l = 0; l < CellNumbs.at(2); l++)

//                } // end of if (sub_polytope_set.size() > 0)
//            } // end of for (auto grain_id : sub_polytope_set)
//        } // end of for (auto grain_id2 : sub_polytope_set)
// REPAIR        cout << "doubled sub faces sequence size:\t=\t" << doubled_sub_faces_sequence.size() << endl;

        /// Full subcomplex faces set
        for (auto unique_sub_faces: doubled_sub_faces_sequence)
            sub_faces_set.insert(unique_sub_faces); // set automatically remove all repetitions
//        } // end of for (auto grain_id: sub_polytope_set)
        cout << "Subcomplex faces set size:\t\t=\t" << sub_faces_set.size() << endl;

        for (auto face_id: doubled_sub_faces_sequence) {
// REPAIR        cout << count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) << endl;
            if (count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) > 1) {
                internal_sub_faces_set.insert(face_id);
            }
        }
        cout << "Internal Subcomplex faces sequence size:\t=\t" << internal_sub_faces_set.size() << endl;
    }  //if (sub_polytope_set.size() > 0) {

// Common subcomplex faces (internal_doubled_sub_faces_sequence)
/*
if (doubled_sub_faces_sequence.size() > 0) {
    bool is_common = 0;
    unsigned int face_counter = 0;

    for (auto face_id: doubled_sub_faces_sequence) {
        is_common = 0;
        cout << "faces_sequence: " << count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) << endl;

        if (count(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), face_id) > 1) {
            internal_doubled_sub_faces_sequence.push_back(face_id);
        }

        face_counter = face_id;
        for (auto face_id2: doubled_sub_faces_sequence) {
            if(face_id2 == face_counter) ++face_counter;
            if(face_counter > 2) is_common = 1;
        }
        if(is_common) cout << "Is common face here !" << endl;
    }
}
 */
/// Common face coordinates
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber) {
        if (std::find(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), fnumber) != doubled_sub_faces_sequence.end())
            internal_faces_coordinates.push_back(find_aGBseed(fnumber, PCCpaths, CellNumbs, polytope_coordinates_vector)); // tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);
    }

    cout << "common faces coordinates size: " << subcomplex_polytope_coordinates.size() << endl;

/// Special faces
///--------------------
    if (all_sfaces_sequence.size() > 0)
        for (unsigned int sface_number : all_sfaces_sequence)
            if (std::find(sub_faces_set.begin(), sub_faces_set.end(), sface_number) != sub_faces_set.end())
//            if (std::find(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), sface_number) != doubled_sub_faces_sequence.end())
                sub_sfaces_sequence.push_back(sface_number);

    cout << "Subcomplex special faces size:\t\t=\t" << sub_sfaces_sequence.size() << endl;

/// Induced subcomplex faces
    if (sub_cfaces_sequence.size() > 0)
        for (unsigned int cface_number : sub_cfaces_sequence)
            if (std::find(doubled_sub_faces_sequence.begin(), doubled_sub_faces_sequence.end(), cface_number) != doubled_sub_faces_sequence.end())
                sub_cfaces_sequence.push_back(cface_number);

    Subcomplex new_cut; // new subcomplex with its ID
/// Setting all quantities to the subcomplex new_subPCC with id = 0
    new_cut.Set_sub_polytope_set(sub_polytope_set);
    new_cut.Set_sub_faces_set(sub_faces_set);
    new_cut.Set_internal_face_coordinates(internal_faces_coordinates);
    new_cut.Set_sub_polytope_coordinates(subcomplex_polytope_coordinates);
    new_cut.Set_sub_sfaces_sequence(sub_sfaces_sequence); //special faces
    new_cut.Set_sub_cfaces_sequence(sub_cfaces_sequence); //cracked (induced) faces

    several_cuts.push_back(new_cut);
    return several_cuts;
}


/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence

/*
std::ofstream Cracked_pcc_out, Cracked_stat_out;
Cracked_pcc_out.open(output_dir + "Macrocrack_faces_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
for(unsigned int cfs : sub_sfaces_sequence) {
///    for(auto cfs : internal_faces_coordinates) {
//        find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector);
cout << "Unique_sub_faces : " << get<0>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) << "\t" << get<1>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) << "\t" << get<2>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) << endl;
Cracked_pcc_out << get<0>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) * 10.0 << "\t" << get<1>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) * 10.0 << "\t" << get<2>(find_aGBseed(cfs, PCCpaths, CellNumbs, polytope_coordinates_vector)) * 10.0 << "\t" << endl;
///        cout << get<0>(cfs) << "\t" << get<1>(cfs) << "\t" << get<2>(cfs) << endl;
///        Cracked_pcc_out << get<0>(cfs) << "\t" << get<1>(cfs) << "\t" << get<2>(cfs) << endl;
}
Cracked_pcc_out.close();

exit(0);
*/