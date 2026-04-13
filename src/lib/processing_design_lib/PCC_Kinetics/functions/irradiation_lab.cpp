///================================ A part of the PCC Kinetics module =============================================================///
///=================================================================================================================================///
/** The library contains functions calculating 'time' from the [0,1] interval for various materials subjected to the process         *
/* of degradation via irradiation damage.                                                                                       */
///================================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// external libraries
#include "../../../external/Eigen/Core"
#include "../../../external/Eigen/SparseCore"

// local libraries
#include "../ini/ini_readers.h"
#include "../ini/ini_materials_reader.h"

#include "../../PCC_Objects.h"
#include "../../PCC_Measures.h"
#include "../../PCC_Support_Functions.h"
#include "../../PCC_Processing/functions/processing_assigned_labelling.h"
#include "../../PCC_Processing/functions/processing_indexing.h"

using namespace std; // standard namespace

/// External variables
extern int PCC_dimension;
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::string output_dir;
extern ofstream irradiation_damaged_output, irradiation_damaged_fractions_output;

#include "irradiation_lab.h"

///* ========================================================= PCC Irradiation Damage ======================================================= *///
///* ========================================================================================================================================= *///

/// ======# 1 #================= void interface_irradiation_damage() function ==============================================================///

/*!
 * @details The function implement simulation of a interface/grain boundary damage process due to irradiation via the PCC faces
 * according to the damage rate depending on interface properties, thermodynamic parameters and local combinatorics of affected cells
 * @param p_cells_history
 */
std::vector<std::vector<double>> interface_irradiation_damage(Config &config, CellDesign &processing_cells_design, std::vector<CellEnergies> &gb_energies){
// change p_cells_history.at(2) -- times when faces (2-cells) change their generated types to 'fractured' due to irradiation damage process
    std::vector<std::vector<double>> p_cells_history;

    int face_cell_type = PCC_dimension - 1;
// a single cell corrosion history
    unsigned int gb_number; // 'gb' is for Grain Boundary or interfaces
    if (PCC_dimension == 3)
        gb_number = CellNumbs.at(2);
    else if (PCC_dimension == 2)
        gb_number = CellNumbs.at(1);

    // Corrosion current in grain boundaries
/// creation_state_f_vector
// from PCC_Support_Functions.cpp
    std::vector<unsigned int> polyhedra_state_vector, face_state_vector; // contains only {0,1,2..} ID values
    polyhedra_state_vector = processing_cells_design.Get_p_design();

//    for(auto op : polyhedra_state_vector)  cout << "op\t" << op << "\t";

    /// Indexing of GBs by Grain types
    face_state_vector = TopDown_cell_indexing(2, polyhedra_state_vector);

///    cout << "T H E R E !\t" << polyhedra_state_vector.size() << endl;

// REPAIR:  for (auto polyhedra : face_state_vector) cout << polyhedra; cout << endl; exit(17);

    // SpMat class defined in main.cpp from the Eigen external library
    SpMat FES(CellNumbs.at(1 + (PCC_dimension - 3)), CellNumbs.at(2 + (PCC_dimension - 3))); // adapted for grain boundaries - either faces in 3-PCC or edges in 2-PCC
    FES = SMatrixReader(paths_to_PCC_matrices.at(5 + (PCC_dimension - 3)), (CellNumbs.at(1 + (PCC_dimension - 3))),
                        (CellNumbs.at(2 + (PCC_dimension - 3)))); //all Edges-Faces

/// BL local normalised face indices
    std::vector<double> irradiation_GB_normalised_coefficients(gb_number,0), CL_normalised_coefficients(gb_number,1);
    std::vector<unsigned int> special_f_sequence = processing_cells_design.Get_f_special_sequence();
    std::vector<unsigned int> special_g_sequence(gb_number,0);
//    if (processing_cells_design.Get_f_induced_sequence().size() > 0)
//        std::vector<unsigned int> special_g_sequence = processing_cells_design.Get_f_induced_sequence();
//    corrosion_GB_normalised_coefficients = face_edge_normalised_local_indices(special_f_sequence, FES); // function from Measures.h
/// TODO:    if (special_g_sequence.size() > 0)
///       CL_normalised_coefficients = face_edge_normalised_local_indices(special_g_sequence, FES); // function from Measures.h
    //GB size
    std::vector<unsigned int> all_face_numbers;
    for (unsigned int fn = 0; fn < CellNumbs.at(face_cell_type); ++fn)
        all_face_numbers.push_back(fn);

    std::vector<std::tuple<double, double, double>> face_barycentres_vector; // coordinate vectors defined globally
    face_barycentres_vector = Tuple3Reader(paths_to_PCC_matrices.at(13)); // grain barycentres
    if (face_barycentres_vector.size() == 0)
        face_barycentres_vector = face_sequence_barycentre_coordinates(all_face_numbers);

    std::vector<double> face_areas_vector;
    const char *cfav = paths_to_PCC_matrices.at(7).c_str(); // face areas in 3-PCC
    face_areas_vector = VectorDReader(cfav);
    // REPAIR for (auto fav : face_areas_vector) cout << fav << endl; exit(0);

    /// Projections (scalar multiplication) of areas to the beam direction
    std::vector<std::tuple<double,double,double>> face_normals;
    face_normals = Tuple3Reader(paths_to_PCC_matrices.at(11));


    std::tuple<double, double, double> beam_direction = config.Get_kinetics_beam_direction();

    std::vector<double> projected_face_areas(gb_number,0);
    for (unsigned int fn = 0; fn < CellNumbs.at(face_cell_type); ++fn) {
        projected_face_areas.at(fn) = face_areas_vector.at(fn) *
                                      abs((std::get<0>(face_normals.at(fn)) * std::get<0>(beam_direction) +
                                       std::get<1>(face_normals.at(fn)) * std::get<1>(beam_direction) +
                                       std::get<2>(face_normals.at(fn)) * std::get<2>(beam_direction)));
    }
    /// read gb sizes, von Mizes (equivalent) stresses and temperatures
    config_reader_multiphysics(config);

    std::tuple<double, double, double> sample_dimensions = config.Get_multiphysics_sample_dimensions(); // [m]

    /// Externally applied Stress and Temperature
    std::vector<double> gb_equivalent_stress(CellNumbs.at(face_cell_type)), gb_temperature(CellNumbs.at(face_cell_type));

    Eigen::MatrixXd est = config.Get_multiphysics_external_stress_tensor(); // external stress tensor
    double ambient_temperature = config.Get_multiphysics_temperature(); // external ambient temperature
    std::tuple<double, double, double, double, double, double, double, double, double>
            external_gb_stress = make_tuple(est(0,0), est(0,1),est(1,2),est(1,0),est(1,1),est(1,1),est(2,0),est(2,1),est(2,2));

    gb_equivalent_stress = gb_energies.at(0).Get_von_Mises_stress();
    gb_temperature = gb_energies.at(0).Get_ambient_temperature();

    config_reader_kinetics(config);
    double irradiation_damage_rate_coeff = config.Get_kinetics_irradiation_damage_rate();
    double beam_energy_flux = config.Get_kinetics_beam_energy_flux();
    double beam_current = config.Get_kinetics_beam_current();
    double energy_dissipation_rate = config.Get_kinetics_energy_dissipation_rate();
    double observation_output_time = config.Get_kinetics_observation_time();

    std::string Mid_matrix = config.Get_kinetics_material_id();
    Material material(Mid_matrix);

    /// METRICS assigning for all grain boundaries
    std::vector<double> gb_volumes(gb_number,0);
    double gb_width = material.Get_gb_width()*std::pow(10,-9); // in nanometres in the material database
    for (unsigned int i = 0; i < CellNumbs.at(face_cell_type); ++i) {
        face_areas_vector.at(i) = face_areas_vector.at(i) * (get<0>(sample_dimensions) * get<1>(sample_dimensions)); /// WARNING! Works well only for cubic samples!
        gb_volumes.at(i) = face_areas_vector.at(i) * gb_width;
        projected_face_areas.at(i) = projected_face_areas.at(i) * (get<0>(sample_dimensions) * get<1>(sample_dimensions));
    }

/// irradiation interfaces DAMAGE RATE
    std::vector<double> interface_irradiation_damage_rate(gb_number);
    double Boltzmann_constant = 1.380649*std::pow(10,-23); // [J/K]
    double Burgers_vecor_module = material.Get_Burgers_vector();
    //    double damage_activation_volume = 1.0*std::pow(10,-30);

/// ================================================ //
///    Irradiation process
/// ================================================ //

/// Initial
    std::vector<double> particle_trapping_probability(gb_number,0);
    std::string material_id = config.Get_kinetics_material_id();
    double austenite_interface_param = 0, martensite_interface_param = 0, ferrite_interface_param = 0, austenite_martensite_interface_param = 0, austenite_ferrite_interface_param = 0, martensite_ferrite_interface_param = 0;
    std::vector<double> radiation_damaged_face_area_fractions(1,0);

    /// material-related parameters
    material_irradiation_reader(material_id, austenite_interface_param, martensite_interface_param, ferrite_interface_param, austenite_martensite_interface_param, austenite_ferrite_interface_param, martensite_ferrite_interface_param);
//cout << austenite_interface_param << "\t\t" << martensite_interface_param << "\t\t" << ferrite_interface_param << "\t\t" << austenite_martensite_interface_param << "\t\t" << austenite_ferrite_interface_param << "\t\t" << martensite_ferrite_interface_param << endl;
// exit(22);
    // GB types::
    // 0 - austenite/austenite interface | 1 - martensite/austenite interface | 2 - martensite/martensite interface
    // 3 - austenite/ferrite | 4 - martensite/ferrite | 6 - ferrite/ferrite
    //----------------------------------------------------------------------------------------------------------------
    for (auto  itr = face_state_vector.begin(); itr != face_state_vector.end(); ++itr) {
            if (*itr == 0)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = austenite_interface_param;
            else if (*itr == 1)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = austenite_martensite_interface_param;
            else if (*itr == 2)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = martensite_interface_param;
            else if (*itr == 3)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = austenite_ferrite_interface_param;
            else if (*itr == 4)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = martensite_ferrite_interface_param;
            else if (*itr == 6)
                particle_trapping_probability.at(std::distance(face_state_vector.begin(),itr)) = ferrite_interface_param;
    }

    double irradiation_time = 0, time_step_coeff = 0, time_step = 0;
    std::vector<double> irradiation_time_vector(gb_number,0), irradiation_gb_damage_rate(gb_number,0), irradiation_gb_damage(gb_number,0);
/// ===================== START OF THE IRRADIATION PROCESS LOOP ============================================ ///
    double PCC_total_face_area = 0;
    for (auto it : face_areas_vector)
        PCC_total_face_area += it;

    do {
        // current computation time
        irradiation_time += time_step;

/// corrosion RATE equation
        std::vector<unsigned int> radiation_face_damages(gb_number,0), damaged_face_sequence(gb_number,0);
        for (unsigned int gbn = 0; gbn < gb_number; ++gbn) {
            irradiation_gb_damage_rate.at(gbn) =
                    irradiation_damage_rate_coeff * (beam_current * projected_face_areas.at(gbn)) * particle_trapping_probability.at(gbn) *
                    (std::pow(Burgers_vecor_module, 3.0) / gb_volumes.at(gbn)); /// add stress&temperature effects *exp(gb_equivalent_stress.at(gbn) * CL_normalised_coefficients.at(gbn) * corrosion_activation_volume / (Boltzmann_constant * gb_temperature.at(gbn)));

// REPAIR: cout << "(std::pow(Burgers_vecor_module, 3.0) / gb_volumes.at(gbn))" << "\t\t" << (std::pow(Burgers_vecor_module, 3.0) / gb_volumes.at(gbn)) << "\t\tbeam_current" << "\t\t" << beam_current << endl;
        }

        /// finding corrosion grain boundary DAMAGE time
        std::vector<double> time_vector(gb_number, 1);
        for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
            if (irradiation_gb_damage_rate.at(gbn) > 0)
                time_vector.at(gbn) = 1.0 /irradiation_gb_damage_rate.at(gbn);
        }
// irradiation TIME STEP
        time_step_coeff = config.Get_kinetics_time_scale(); // taken from 'kinetic_time_scale' in the config/Kinetic.ini file
        time_step = time_step_coeff * (*std::min_element(time_vector.begin(), time_vector.end()));

/// irradiation GB DAMAGE
        for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
            irradiation_gb_damage.at(gbn) += time_step * irradiation_gb_damage_rate.at(gbn);

// REPAIR:            cout << "\t" << "irradiation_gb_damage.at(gbn):\t" << time_step * irradiation_gb_damage_rate.at(gbn) << "\t" << endl;

            if (irradiation_gb_damage.at(gbn) > 1 && special_g_sequence.at(gbn) == 0) {
                special_g_sequence.at(gbn) = 8;
                irradiation_time_vector.at(gbn) = irradiation_time;

                radiation_damaged_face_area_fractions.push_back(radiation_damaged_face_area_fractions.back() + (face_areas_vector.at(gbn)/PCC_total_face_area));
                cout << radiation_damaged_face_area_fractions.back() << "\t" << irradiation_time << endl;

                //corrosion_output
                //cout << gbn << "\t" << irradiation_time_vector.at(gbn) << "\t" << get<0>(face_barycentres_vector.at(gbn)) << "\t" << get<1>(face_barycentres_vector.at(gbn)) << "\t" << get<2>(face_barycentres_vector.at(gbn)) << endl;
                irradiation_damaged_output << gbn << "\t" << irradiation_time_vector.at(gbn) << "\t" << get<0>(face_barycentres_vector.at(gbn)) << "\t" << get<1>(face_barycentres_vector.at(gbn)) << "\t" << get<2>(face_barycentres_vector.at(gbn)) << endl;
                irradiation_damaged_fractions_output << radiation_damaged_face_area_fractions.back() << "\t" << irradiation_time << endl;

            }

        } // end of for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {

/* /// 'cout' check
        for (auto cs_itr = 0; cs_itr < corrosion_time_vector.size(); ++cs_itr) {
            if (corrosion_time_vector.at(cs_itr) > 0)
                cout << "corrosion damage time at grain boundary #\t"
                     << std::distance(corrosion_time_vector.begin(), corrosion_time_vector.begin() + cs_itr) << "\tis equal to\t" << corrosion_time_vector.at(cs_itr) << endl;
        }
*/
        cout << "\t" << "Full irradiation process Time:\t" << irradiation_time << "\t" << "with the current Time Step:\t" << time_step << endl;
        cout << "\t" << "Zero elements in the 'irradiation damage time vector':\t" << std::count(irradiation_time_vector.begin(), irradiation_time_vector.end(), 0) << "\t out of\t" << CellNumbs.at(face_cell_type) << endl;

    } while( std::count(irradiation_time_vector.begin(), irradiation_time_vector.end(), 0) > 0.8*irradiation_time_vector.size() ); // END do while();
// } while( irradiation_time < config.Get_kinetics_observation_time()); // END do while();

/// result
    for (unsigned int gbn = 0; gbn < gb_number; ++gbn) {
        p_cells_history.push_back(
                {double(gbn), irradiation_time_vector.at(gbn), radiation_damaged_face_area_fractions.back()});
    }

    return p_cells_history;
} // END of bulk_irradiation()