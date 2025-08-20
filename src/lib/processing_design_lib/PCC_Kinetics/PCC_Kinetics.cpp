///========================================= PCC Kinetics module ================================================================================ ///
///============================================================================================================================================= ///
///* The interface use functions from PCC_Kinetics/functions C++ libraries to generate for each p-cell in a PCC the moment of 'time'           *///
///* in the [0,1] range when it changed its special 'generated' type as the result of a 'kinetic' process (e.g. corrosion or irradiation)     *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2025 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include "../../../src/lib/external/Eigen/SparseCore"

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../ini/ini_readers.h"
#include "../PCC_Objects.h"
#include "../PCC_Measures.h"

// Local
///---------------------------------------------------------
#include "functions/corrosion_lab.h"
#include "functions/irradiation_lab.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs;
extern ofstream kinetics_logfile_stream;
extern std::string source_path;
extern std::string output_dir;
extern std::vector<std::string> paths_to_PCC_matrices;
extern int PCC_dimension;

#include "PCC_Kinetics.h"
///* ========================================================= PCC KINETICS FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Employ functions from PCC_Kinetics/functions C++ libraries to generate for each p-cell in a PCC the moment of 'time'
 * in the [0,1] range when it changed its special 'generated' type as the result of a 'kinetic' process (e.g. corrosion or irradiation)
 * All the initial settings are written in 'kinetics.ini' file, including the particular material.
 * @param configuration
 * @return std::vector<vector<double>> p_cells_history
 */
///* The interface use functions from PCC_Kinetics/functions C++ libraries to generate for each p-cell in a PCC the moment of 'time'           *///
///* in the [0,1] range when it changed its special 'generated' type as the result of a 'kinetic' process (e.g. corrosion or irradiation)     *///


std::vector<vector<double>> PCC_Kinetics(Config &kinetics_configuration, CellDesign &processing_cells_design) {
/// Main output of the module 'p_cells_history' -- a vector contained the 'time' moment in the [0,1] range when each p-cell special generated p-cell changed its type.
/// First line 0-cells, Second line 1-cells, Third line 2-cells, and Fourth line - 3-cells (for 3-PCCs only)
/// The length of each line corresponds to the number of cells in the PCC.

// Function output
    std::vector<std::vector<double>> p_cells_history;

    cout << "=========================================================================" << endl << endl;
    kinetics_logfile_stream << "==============================================================================================================================================================" << endl << endl;

    config_reader_kinetics(kinetics_configuration);

    std::string material_id = kinetics_configuration.Get_kinetics_material_id();
    Material studied_material(material_id);

    if(kinetics_configuration.Get_kinetics_fk_mode() == "C"s) {
        p_cells_history = surface_interface_corrosion(kinetics_configuration, studied_material, processing_cells_design);

    }
    else if(kinetics_configuration.Get_kinetics_fk_mode() == "I"s){
        p_cells_history = interface_irradiation_damage(kinetics_configuration, studied_material, processing_cells_design);
    }
    return p_cells_history;
}
