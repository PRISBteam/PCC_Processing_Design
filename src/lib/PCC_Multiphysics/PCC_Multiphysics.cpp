///================================ PCC Multi-physics module =================================================================================///
///=========================================================================================================================================///
/** The function in this library assign scalar and vector quantities to PCC k-cells associated with physical energies.                     **/
/**  In particular, elastic energies, thermal (internal), magnetic, etc.                                                                  **/
///* -----------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design code (CPD code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 (https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1) in 2022-2023 years *///
///===================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
/// Attached user-defined C++ libraries:
// External
#include <Eigen/SparseCore>

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
//#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
#include "functions/Multiphysics_functions.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs;
extern ofstream Out_logfile_stream;
extern string source_path;
extern int dim;
extern std::vector<double> face_areas_vector;
extern string source_dir, output_dir, e_mode;

#include "PCC_Multiphysics.h"
///================================ PCC MULTIPHYSICS FUNCTIONS =================================================================================///
///* ========================================================================================================================================= *///
std::vector <double> PCC_Multiphysics(Macrocrack &large_crack, double external_vonMizes_stress) { // std::vector<unsigned int> &s_faces_sequence, string K_type
    std::vector <double> new_face_energies(CellNumbs.at(2),0.0); // initial values 0.0
//    std::vector <double> sample_dimensions;
    //double average_grain_size = 500.0*pow(10,-9); // [metres] | 500 nm grains (!)
    int grains_in_a_row = round(pow(CellNumbs.at(3),0.3333)); /// (!) currently the same for 3 directions

    /// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
// The source directory and simulation type from file config.txt
    string S_type; // 'P' or 'H' :: This char define the subsection type: 'P' for the whole Plane cut, 'H' for the half-plane cut like a crack
    //std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string &main_type, std::string &e_type);
    std::vector<int> ConfigVector = config_reader_main(source_path, source_dir, output_dir, main_type, e_type);

    std::tuple<double, double, double> sample_dimensions;
    config_reader_multiphysics(source_path, sample_dimensions, Out_logfile_stream);

    bool SubcomplexON(char* config, bool time_step_one); // Check the Subcomplex (Section) module status (On/Off) in the config.txt file

    int set_crack_mode = 1;
    double Puasson_coeff = 0.3; /////////// TEMPORARILY (!)

    ///setting crack lengths (!)
    large_crack.Set_real_crack_length(sample_dimensions.at(1)); //Set_real_crack_length(double sample_size)
    double new_crack_length = large_crack.Get_real_crack_length();
//REPAIR
    cout << "grains_in_a_row: " << grains_in_a_row << " sample_dimensions(0): " << sample_dimensions[0] << " macrocrack length: " << new_crack_length << endl;
    /// Call the main function
    crack_modes_stress_field(new_face_energies, set_crack_mode, new_crack_length, external_vonMizes_stress, Puasson_coeff);
    /// Surface energy of a macrocrack
    double sface_energy_matrix = 2.0, adhesion_energy_rGO = 1.0; /// crack real energies  here (!)

//     large_cracks_vector.at(itr).surface_energy = 2.0*new_crack_length*sample_dimensions.at(0)*sface_energy_matrix; // 2 (two surfaces)*L_y*Size_x
    large_crack.surface_energy = 2.0 * new_crack_length * sample_dimensions.at(0) * sface_energy_matrix; // 2 (two surfaces)*L_y*Size_x
    double bridging_coeff = 4.0; /// show how many GBs contributes to the bridging effect

    double total_in_crack_sfaces_area = 0.0;
    for (unsigned int gb : large_crack.Get_sfaces_sequence())
        total_in_crack_sfaces_area  += face_areas_vector.at(gb); // total area of all sfaces
    large_crack.bridging_energy = bridging_coeff * adhesion_energy_rGO * total_in_crack_sfaces_area * sample_dimensions.at(0) * sample_dimensions.at(1);
    //large_cracks_vector.at(itr).bridging_energy = bridging_coeff*adhesion_energy_rGO*total_in_crack_sfaces_area*sample_dimensions.at(0)*sample_dimensions.at(1);

//REPAIR    cout << " macro_crack.bridging_energy: " << large_cracks_vector.at(0).bridging_energy << " macro_crack.surface_energy: " << large_cracks_vector.at(0).surface_energy <<endl; ++itr;  } // end of for( auto macro_crack : large_cracks_vector)

    return new_face_energies;
};