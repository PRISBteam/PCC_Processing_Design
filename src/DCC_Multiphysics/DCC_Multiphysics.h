///================================ DCC Multiphysics module ===================================================================================///
///========================================================================================================================================///
/** The function in this library assign scalar and vector quantities to PCC k-cells associated with some physical parameters.             ///
*  In particular, elastic energies, temperature, conductivity, etc.                                                          **/         ///
///=====================================================================================================================================///

/// Attached user-defined C++ libraries:
///-------------------------------------
#include "Multiphysics_Functions.h"
///-------------------------------------

using namespace std;
//using namespace Eigen;

/// User-defined function definitions:
/* No one to date */

std::vector <double> DCC_Multiphysics(std::vector <macrocrack> large_cracks_vector, double external_vonMizes_stress) { // std::vector<unsigned int> &s_faces_sequence, string K_type
    std::vector <double> new_face_energies;

    int set_crack_mode = 1;
    double Puasson_coeff = 0.3;

    vector<tuple<double, double, double>> face_coordinates;

    for( auto macro_cracks : large_cracks_vector) {
        double new_crack_length = macro_cracks.Get_crack_length(macro_cracks.crack_id);
        face_coordinates = macro_cracks.Get_common_faces_coordinates(0); /// 0 is not good here !

        /// Call the main function
        crack_modes_stress_field(new_face_energies, set_crack_mode, new_crack_length, external_vonMizes_stress, face_coordinates, Puasson_coeff);
    }

    return new_face_energies;
};


/*
    /// All grain coordinates
    string GCpath_string = input_folder + "grain_seeds.txt"s;
    char* GCpath = const_cast<char*>(GCpath_string.c_str());
    vector<tuple<double, double, double>> grain_coordinates = TuplesReader(GCpath);

    /// All face coordinates
    vector<tuple<double, double, double>> face_coordinates;
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber)
            face_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates));
*/
