///================================ A part of the PCC Kinetics module =============================================================///
///=================================================================================================================================///
/** The library contains functions calculating 'time' from the [0,1] interval for various materials subjected to the process         *
/* of degradation via corrosive environments.                                                                                       */
///================================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// external libraries
#include "../../../external/Eigen/Core"
#include "../../../external/Eigen/SparseCore"

// local libraries
#include "../ini/ini_readers.h"
#include "../../PCC_Objects.h"
#include "../../PCC_Measures.h"
#include "../../PCC_Support_Functions.h"
#include "../../PCC_Processing/functions/processing_assigned_labelling.h"

using namespace std; // standard namespace

/// External variables
extern int PCC_dimension;
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::string output_dir;

#include "corrosion_lab.h"

///* ========================================================= PCC Corrosion Degradation ======================================================= *///
///* ========================================================================================================================================= *///

/// ======# 1 #================= void surface_interface_corrosion() function ==============================================================///

/*!
 * @details The function implement simulation of a degradation/corrosion process progressing via the PCC faces (material interfaces)
 * according to the degradation rate depending on interface properties, thermodynamic parameters and local combinatorics of affected   
 * @param p_cells_history
 */
std::vector<std::vector<double>>  surface_interface_corrosion(Config &config, Material &material, CellDesign &processing_cells_design){
// change p_cells_history.at(2) -- times when faces (2-cells) change their generated types to 'fractured' due to corrosion process
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
    std::vector<int> face_state_vector; // contains only {0,1,2..} ID values
    face_state_vector = state_vector_by_sequence(processing_cells_design.Get_f_special_sequence(), face_cell_type);

    std::vector<double> gb_corrosion_current(gb_number);
    for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
        if (face_state_vector.at(gbn) == 0) {
            gb_corrosion_current.at(gbn) = material.Get_lagbs_corrosion_current();
        }        else if (face_state_vector.at(gbn) == 1)
            gb_corrosion_current.at(gbn) = material.Get_hagbs_corrosion_current();
        else if (face_state_vector.at(gbn) == 2)
            gb_corrosion_current.at(gbn) = material.Get_sigma3_corrosion_current();
    }

    // SpMat class defined in main.cpp from the Eigen external library
    SpMat FES(CellNumbs.at(1 + (PCC_dimension - 3)), CellNumbs.at(2 + (PCC_dimension - 3))); // adapted for grain boundaries - either faces in 3-PCC or edges in 2-PCC
    FES = SMatrixReader(paths_to_PCC_matrices.at(5 + (PCC_dimension - 3)), (CellNumbs.at(1 + (PCC_dimension - 3))),
                        (CellNumbs.at(2 + (PCC_dimension - 3)))); //all Edges-Faces

/// BL local normalised face indices
    std::vector<double> corrosion_GB_normalised_coefficients(gb_number,0), CL_normalised_coefficients(gb_number,1);
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
    //    const char *fncv = paths_to_PCC_matrices.at(13).c_str(); // face barycentres
   // REPAIR for (auto fav : face_areas_vector)  cout << fav << endl;  exit(0);

    /// read gb sizes, von Mizes (equivalent) stresses and temperatures
    config_reader_multiphysics(config);

    std::tuple<double, double, double> sample_dimensions = config.Get_multiphysics_sample_dimensions(); // [m]

    std::vector<double> gb_equivalent_stress(CellNumbs.at(face_cell_type)), gb_temperature(CellNumbs.at(face_cell_type));

    std::vector<CellEnergies> gb_energies(gb_number);
    Eigen::MatrixXd est = config.Get_multiphysics_external_stress_tensor(); // external stress tensor
    double ambient_temperature = config.Get_multiphysics_temperature(); // external ambient temperature
    std::tuple<double, double, double, double, double, double, double, double, double>
            external_gb_stress = make_tuple(est(0,0), est(0,1),est(1,2),est(1,0),est(1,1),est(1,1),est(2,0),est(2,1),est(2,2));

    for (unsigned int i = 0; i < CellNumbs.at(face_cell_type); ++i){
        gb_energies.at(i).Set_von_Mises_stress(external_gb_stress);
        gb_energies.at(i).Set_ambient_temperature(ambient_temperature);
    }

    //assigning for all grain boundaries
    for (unsigned int i = 0; i < CellNumbs.at(face_cell_type); ++i) {
        face_areas_vector.at(i) = face_areas_vector.at(i) * (get<0>(sample_dimensions) * get<1>(sample_dimensions)); /// WARNING! Works well only for cubic samples!
        gb_equivalent_stress.at(i) = gb_energies.at(i).Get_von_Mises_stress();
        gb_temperature.at(i) = gb_energies.at(i).Get_ambient_temperature();
    }

/// corrosion rate
    std::vector<double> gb_corrosion_rate(gb_number);
    double Boltzmann_constant = 1.380649*pow(10,-23); // [J/K]

    ///! coefficient under stress exponent
    double corrosion_activation_volume = 1.0*pow(10,-30);

///    Corrosion process
//================================================
    std::vector<unsigned int> special_corr_sequence; // corrosive grain boundaries

// only surface GBs in the set
    std::vector<unsigned int> surface_gb_set;
    for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
        if (get<2>(face_barycentres_vector.at(gbn)) == 0)
            surface_gb_set.push_back(gbn);
    }

// initial
    unsigned int NewCellNumb = 0;
    /// initial corrosion settings
    double corr_gb_initial_fraction = 0.2; /// 20% "hardcoded" initial condition on surface of the sample
    for (unsigned int gbn = 0; gbn < corr_gb_initial_fraction*surface_gb_set.size(); ++gbn) {
        NewCellNumb = NewCellNumb_R(surface_gb_set.size()); /// advanced NewCellNumb_R generator of special cell IDs
        special_corr_sequence.push_back(surface_gb_set.at(NewCellNumb)); // defined in the assigned labelling library; "\functions" subfolder
    }

    double corrosion_time = 0.0;
    std::vector<double> corrosion_gb_damage(gb_number);
    double time_step_coeff = 0, time_step = 0;
    std::vector<double> corrosion_time_vector(gb_number,0);
/// ===================== START OF THE SURFACE CORROSION LOOP ============================================ ///
    do {
        corrosion_time += time_step;

//! The "BL" grain boundary index counting the weighted amount of their local corrosive GB neighbours
        corrosion_GB_normalised_coefficients = face_edge_normalised_local_indices(special_corr_sequence, FES); // function from Measures.h

/// corrosion RATE equation
        for (unsigned int gbn = 0; gbn < gb_number; ++gbn) {
            gb_corrosion_rate.at(gbn) = gb_corrosion_current.at(gbn) * corrosion_GB_normalised_coefficients.at(gbn); /// add stress&temperature effects *exp(gb_equivalent_stress.at(gbn) * CL_normalised_coefficients.at(gbn) * corrosion_activation_volume / (Boltzmann_constant * gb_temperature.at(gbn)));
//                cout << "gb_corrosion_rate.at(gbn)" << "\t\t" << gb_corrosion_rate.at(gbn) << endl;
//                cout << "face_sizes.at(gbn)" << "\t\t" << face_areas_vector.at(gbn) << endl;

        /// update the number of corrosive GBs in the PCC
            if (gb_corrosion_rate.at(gbn) > 0 &&
                std::find(special_corr_sequence.begin(), special_corr_sequence.end(), gbn) == special_corr_sequence.end())
                special_corr_sequence.push_back(gbn);
        }

    /// finding corrosion grain boundary DAMAGE time
        std::vector<double> time_vector(gb_number, 1);
        for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
            if (gb_corrosion_rate.at(gbn) > 0)
                time_vector.at(gbn) = face_areas_vector.at(gbn) / gb_corrosion_rate.at(gbn);
        }
/// corrosion TIME STEP
        time_step_coeff = config.Get_kinetics_time_scale(); // taken from 'kinetic_time_scale' in the config/Kinetic.ini file
        time_step = time_step_coeff * (*std::min_element(time_vector.begin(), time_vector.end()));

/// corrosion GB DAMAGE
        for (unsigned int gbn = 0; gbn < CellNumbs.at(face_cell_type); ++gbn) {
            corrosion_gb_damage.at(gbn) =
                    corrosion_gb_damage.at(gbn) + time_step * gb_corrosion_rate.at(gbn) / face_areas_vector.at(gbn);
/// updating corrosion damage vector
            if (corrosion_gb_damage.at(gbn) > 1 && special_g_sequence.at(gbn) == 0) {
                special_g_sequence.at(gbn) = 5;
                corrosion_time_vector.at(gbn) = corrosion_time;
            }
        }
/// 'cout' check
        for (auto cs_itr = 0; cs_itr < corrosion_time_vector.size(); ++cs_itr) {
            if (corrosion_time_vector.at(cs_itr) > 0)
                cout << "corrosion damage time at grain boundary #\t"
                     << std::distance(corrosion_time_vector.begin(), corrosion_time_vector.begin() + cs_itr) << "\tis equal to\t" << corrosion_time_vector.at(cs_itr) << endl;
        }

        cout << "\t" << "Full corrosion process Time:\t" << corrosion_time << "\t" << "with the current Time Step:\t" << time_step << endl;
        cout << "\t" << "Zero elements in the 'corrosion time vector':\t" << std::count(corrosion_time_vector.begin(), corrosion_time_vector.end(), 0) << "\t out of\t" << CellNumbs.at(face_cell_type) << endl;

//Repair        cout << "special_corr_sequence size\t" << special_corr_sequence.size() << "\t" << "BL 0s number \t" << std::count(corrosion_GB_normalised_coefficients.begin(), corrosion_GB_normalised_coefficients.end(), 0) << endl;
//Repair        cout << "BL size\t" << 1.0 - std::count(corrosion_GB_normalised_coefficients.begin(), corrosion_GB_normalised_coefficients.end(), 0)/ double(gb_number) << endl;

///            } while( corrosion_time < 7.0*pow(10,-12)); // END do while( corrosion_time < 1.0 )
    } while( std::count(corrosion_time_vector.begin(), corrosion_time_vector.end(), 0) != 0 ); // END do while( corrosion_time < 1.0 )

/// result
    for (unsigned int gbn = 0; gbn < gb_number; ++gbn) {
        p_cells_history.push_back(
                {double(gbn), corrosion_time_vector.at(gbn), get<2>(face_barycentres_vector.at(gbn))});
    }

    return p_cells_history;
} // END of surface_interface_corrosion()