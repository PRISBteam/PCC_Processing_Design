///======================================= PCC Multiphysics module ================================================================================ ///
///=============================================================================================================================================== ///
///* The interface use functions from Multiphysics_<***>_functions.h C++ libraries allowing to set physical scales and the corresponding k-cell  *///
///* elastic, thermal and other type of energies including the stress concentrators internal fields described in the                            *///
///* multiphysics_internal_stresses.cpp. The module also employs Subcomplexes, like a half-planes associated with the presence of macrocracks. *///                                                                                                                   *///
///* -----------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code)   *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years     *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                    *///
///======================================================================================================================================= ///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include <Eigen/SparseCore>

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../ini/ini_readers.h"
#include "../PCC_Objects.h"
#include "../PCC_Measures.h"
#include "../PCC_Subcomplex/functions/subcomplex_cross_section.h"

// Local
///---------------------------------------------------------
#include "functions/multiphysics_internal_stresses.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs;
extern std::vector<double> face_areas_vector;
extern ofstream Out_logfile_stream;
extern string source_path;
extern std::vector<std::string> PCCpaths;
extern int dim;

#include "PCC_Multiphysics.h"
///* ========================================================= PCC MULTIPHYSICS FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
* @details C++ libraries allowing to set physical scales and the corresponding k-cell elastic, thermal and other type of energies including the stress concentrators internal fields described in the multiphysics_internal_stresses.cpp.
 * The module also employs Subcomplexes, like a half-planes associated with the presence of macrocracks. The main output of the module Cell Energies (CE) class object. In particular, it contains (1) elastic external energies (and stress state),
 * (2) elastic internal energies (and stress state), (3) thermal cell energies. All the initial settings are written in 'multiphysics.ini' file, including the Material ID, Cuboid sample dimensions, Time scale, and external Stress tensor components
 * The Material ID corresponding to the materials Database contained locally in the 'CPD_material_database' directory.
 * @param configuration
 * @param pcc_subcomplexes
 * @return
 */
std::vector<CellEnergies> PCC_Multiphysics(Config &configuration, std::vector<Subcomplex> &pcc_subcomplexes, std::vector<Macrocrack> &crack_growth_series) {
/// Main output of the module Cell Energies (CE) class object. In particular, it contains (1) elastic external energies (and stress state), (2) elastic internal energies (and stress state), (3) thermal cell energies
    std::vector<CellEnergies> CE_vector;
    CellEnergies face_energies;

    std::vector<double> self_f_energy_densities, self_f_energies;
    std::vector<double> current_f_energy_densities, current_f_energies;

    std::string Mid_matrix, Mid_inclusion;
    std::tuple<double, double, double> sample_dimensions; // [m]
    double tau; // [seconds]
    Eigen::MatrixXd external_stress_tensor(3,3); // [Pa]
    std::vector<double> macrocrack_ini;

    const char *cfav = PCCpaths.at(7).c_str(); const char *cncv = PCCpaths.at(10).c_str();
    face_areas_vector = VectorDReader(cfav);
// REPAIR for (auto fav : face_areas_vector)  cout << fav << endl;  exit(0);
    config_reader_multiphysics(source_path, Mid_matrix, Mid_inclusion, sample_dimensions, tau, external_stress_tensor, macrocrack_ini, Out_logfile_stream);

    Material matrix_material(Mid_matrix); // set a material from the '../config/CPD_material_database/'
//  Material inclusion_material(Mid_inclusion); // set a material from the '../config/CPD_material_database/'

/// Cohesion energy for each face (2-cell/grain boundary)
    for(unsigned int i = 0; i < CellNumbs.at(2); ++i) // set the SAME cohesion energy
        self_f_energies.push_back(matrix_material.Get_gb_cohesion_energy()*face_areas_vector.at(i));
/// CellEnergy parameter 'vector<double> f_self_energies' for all 2-cells
    face_energies.Set_f_self_energies(self_f_energies);

    /// ======================================================= MACROCRACKS =================================
    double    crack_stress_mode = macrocrack_ini.at(1);
    double    max_crack_lenghts = macrocrack_ini.at(2);
    double    min_crack_lenghts = macrocrack_ini.at(3);
    double    number_of_crack_sizes = 1;
    if(macrocrack_ini.at(4) > 1)
        number_of_crack_sizes = macrocrack_ini.at(4);

    for(int subn = 0; subn < pcc_subcomplexes.size(); ++subn) { /// LOOP over all MACROCRACKS in a PCC

        for (int cl = 1; cl<= number_of_crack_sizes; ++cl) {
            double current_crack_length = (max_crack_lenghts - min_crack_lenghts) * double(cl) / number_of_crack_sizes;

            /// Half-plane individual sub-complex for each of the macrocrack lengths
            Subcomplex half_plane_crack_pcc = Get_half_plane(pcc_subcomplexes.at(subn), current_crack_length);

            cout << " half_plane_crack_pcc.Get_sub_sfaces_sequence().size() " << half_plane_crack_pcc.Get_sub_sfaces_sequence().size() << endl;
            cout << " half_plane_crack_pcc.Get_sub_sfaces_coord().size() " << half_plane_crack_pcc.Get_sub_sfaces_coord().size() << endl;

            /// Creation of the series of macrocracks - each with the corresponding length parameter and a subcomplex as the element of the Macrocrack object
            Macrocrack new_crack((cl - 1), subn, half_plane_crack_pcc, current_crack_length, crack_stress_mode); // (1) macrocrack ID, (2) half-plane subcomplex with length

            cout << " new_crack.Get_sfaces_sequence().size() " << new_crack.Get_sfaces_sequence().size() << endl;
            cout << " new_crack.Get_sfaces_coordinates().size() " << new_crack.Get_sfaces_coordinates().size() << endl;
            crack_growth_series.push_back(new_crack);

            }

            current_f_energy_densities = Multiphysics_crack_stress_field(crack_growth_series.back(), matrix_material, external_stress_tensor, sample_dimensions);

//REPAIR            cout << " Energies: "; for (double cfd : current_f_energy_densities) { std::cout << cfd << "  "; } cout << endl;

            current_f_energies.clear();
        /// Claculation energy density per unit AREA, like [Energy*Vol/Area]
        double av_grain_size = get<0>(sample_dimensions)/std::round(std::pow(CellNumbs.at(3),0.3333));
            for (unsigned int i = 0; i < CellNumbs.at(2); ++i) // set face associated VOLUMETRIC elastic energies (!)
                current_f_energies.push_back(current_f_energy_densities.at(i)*av_grain_size); /// from density to the full GB elastic energy (!
                /// current_f_energies.push_back(current_f_energy_densities.at(i) * face_areas_vector.at(i) * matrix_material.Get_gb_width()); /// from density to the full GB elastic energy (!

//            cout << " Energies: ";
//            for (double cfe : current_f_energies) {
//                std::cout << cfe << "  ";         }
//            cout << endl;

/// CellEnergy parameter 'vector<double> current_f_energies' for all 2-cells
            face_energies.Set_f_elastic_energies(current_f_energies);
            CE_vector.push_back(face_energies);

    } // end of for(pcc_subcomplexes.size())

/*

    int set_crack_mode = 1;
    double Puasson_coeff = 0.3;
    /// Set global sample size vector here (!)
    double average_grain_size = 500.0*pow(10,-9); // [metres] | 500 nm grains (!)
    int grains_in_a_row = round(pow(CellNumbs.at(3),0.3333)); /// (!) currently the same for 3 directions
    sample_dimensions = {grains_in_a_row * average_grain_size, grains_in_a_row * average_grain_size, grains_in_a_row * average_grain_size }; // in meters!

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
*/

/*
    /// All grain coordinates
    string GCpath_string = input_dir + "grain_seeds.txt"s;
    char* GCpath = const_cast<char*>(GCpath_string.c_str());
    vector<tuple<double, double, double>> grain_coordinates = TuplesReader(GCpath);

    /// All face coordinates
    vector<tuple<double, double, double>> face_coordinates;
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber)
            face_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates));
*/

    return CE_vector;
}