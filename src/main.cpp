///******************************************************************************************************************************///
///************************   Polytopal Cell Complex (PCC) Processing Design :: (CPD code) (c)   *******************************///
///****************************************************************************************************************************///
///*                                        Version 5.0 | 12/06/2025                                                         *///
///**************************************************************************************************************************///
///************************************ Dr Elijah Borodin, Manchester, UK **************************************************///
///**************************************** Spring 2022 - Summer 2025  ****************************************************///
///***********************************************************************************************************************///
///*
///*    Code source:    https://github.com/PRISBteam/PCC_Processing_Design/
///*    Documentation:  https://prisbteam.github.io/
///*    PCC sources:    https://materia.team/
///*
///*  The project provides a reliable tool for (1) Generating, (2) Analysing and (3) Optimising of 'design vectors' as the sequences                      *///
///*  of k-cells containing in the k-skeletons, where k = {0,1,2,3}, of a Polytopal Cell Complex (PCC). Such PCCs can be created by external             *///
///*  codes based on the tessellation of 2D or 3D spaces by an agglomeration of polytopes (polygons in the 2D case or polyhedrons in 3D).               *///
///*  Graphs and networks (without loops) are considered as 1-complexes (1-PCCs) and also available for analysis similarly to the 2D and 3D cases.     *///

///* Key terminology:                                                                                                                                                                                                 *///
/// Tessellation's elements     ::   'nodes, 'edges', 'faces', 'polytopes' (with their measures - lengths, areas and volumes - and barycenter coordinates)                                                            ///
/// PCC's elements              ::   'k-cells' containing in 'k-skeletons', where k = {0,1,2,3}, with their types, fractions, and incident (k-1)-cells and (k+1)-cells.                                              ///
/// Material's elements         ::   'quadruple points', 'grain boundary junctions', 'grain boundaries', and 'grains' (with their orientations and barycenter coordinates often get by EBSD or X-ray analysis)      ///

///* ----------------------------------------- *
///* Standard C++ (STL) libraries
///* ----------------------------------------- *
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <set>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries:
///* ------------------------------------------------------------------------------- *

/// Eigen source: https://eigen.tuxfamily.org/ (2024)
/* Alternative way - the libraries must be preliminary copied in the local STL directory (!)
/* #include <Eigen/Core> #include <Eigen/Dense> #include <Eigen/SparseCore> */
#include "../src/lib/external/Eigen/Core"
#include "../src/lib/external/Eigen/Dense"
#include "../src/lib/external/Eigen/SparseCore"

/// Spectra source: https://spectralib.org/ (2024)
/* Alternative way - the libraries must be preliminary copied in the local STL directory (!)
/* #include <Spectra/GenEigsSolver.h> #include <Spectra/SymEigsSolver.h> */
#include "../src/lib/external/Spectra/GenEigsSolver.h"
#include "../src/lib/external/Spectra/SymEigsSolver.h"

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
// Included only in the parallelized version of the code.

///------------------------------------------
using namespace std; // standard/STL namespace
using namespace Eigen; // Eigen library namespace
using namespace Spectra; // Spectra library namespace

/// Eigen library-based classes
typedef Triplet<double> Tr; // <Eigen> library class, which declares a triplet type with the nickname 'Tr' as the objects in the form T = T(i, j, value), where i and j are element's a(i,j) indices in the corresponding dense matrix and the third variable is its value
typedef SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'
typedef MatrixXd DMat; // <Eigen> library class, which declares a dense matrix type of doubles with the nickname 'DMat'

/// * ---------------------------------------------------------------------------------------------------------- *///
/// * ======================================== GLOBAL PROJECT VARIABLES ======================================== *///
/// * ------ Declaration of GLOBAL variables which can be seen in all the project modules and libraries -------- *///

/// Technical variables::
std::string source_path = "../config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // 'source_path' is a path to the directory read from the 'config/main.ini' file
std::string main_type;  // 'mode' parameter read from the config/main.ini file:
                        //| * 'LIST' for the execution one by one all the active (ON in the config file) project modules;
                        //| * 'TUTORIAL' as a specific educational mode;
                        //| * 'PERFORMANCE_TEST' as a special test for a computer's performance and its ability to work with large PCCs, and
                        //| * 'TASK' mode, where user-defined task scripts described in separate 'tasks/*.cpp' files are included with all the necessary modules and functions from the project's libraries.

std::vector<std::string> paths_to_PCC_matrices;   // The vector containing the paths to all the PCC's matrices, measures and other supplementary data files
std::string source_dir, output_dir;     // Input directory [source_dir] to be used if the initial configuration must be read from the file, and output directory [output_dir] for the Writer module and the project Log file as it is written in the 'config/main.ini' file
std::string simulation_tasks_dir;       // Path to the corresponding 'tasks/*.cpp' file containing a 'simulation task' (for 'TASK' execution mode only!) as it is written in the 'config/main.ini' file

// The execution time interval variables for different parts (modulus) of the CPD code:
double Main_execution_time = 0.0, Subcomplex_execution_time = 0.0, Multiphysics_execution_time = 0.0, Processing_execution_time = 0.0, Characterisation_execution_time = 0.0, Design_execution_time = 0.0, Writer_execution_time = 0.0, Kinetics_execution_time = 0.0;

/// Global 'log.txt' file output
std::ofstream main_logfile_stream, subcomplex_logfile_stream, multiphysics_logfile_stream, processing_logfile_stream, characterisation_logfile_stream, design_logfile_stream, writer_logfile_stream;
// 'Processing_Design.log' file output of the entire computation process as a not exact copy of the console output
std::ofstream Out_logfile_stream;
// TODO: DELETE obsolete 'Out_logfile_stream'
std::ofstream corrosion_damaged_output, corrosion_damaged_fractions_output, corrosion_affected_output, corrosion_affected_fractions_output, face_barycentre_coord_outstream, edge_barycentre_coord_outstream;

/// PCC:: (- related variables)
int PCC_dimension;                  // Tessellation dimension corresponding to the maximal value of 'k' in the PCC's k-cell ranks:
    // *  1 -- for graphs and networks,
    // *  2 -- for the 2D plane polygonal tessellations, and
    // *  3 -- for the 3D bulk volumetric tessellations, as it is specified in the 'config/main.ini' file.

/// Combinatorial::
std::vector<unsigned int> CellNumbs;    // vector named CellNumbs containing the numbers of k-cells of different types 'k'. It is read from the 'number_of_cells.txt' file of a PCC (see the corresponding PCC standard for more details).
    //  first line in the CellNumbs.txt file must be the number of nodes (0-cells),
    //  second - edges (1-cells),
    //  third - faces (2-cells) (in the 2D and 3D cases only),
    //  fourth - polyhedra (3-cells) (in the 3D case only).

/// Geometry::
std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // vectors containing barycenter Cartesian coordinates of the corresponding tessellation's elements
// Global vectors of Cartesian coordinates for: (1) vertex coordinates, (2) barycentres of edges, (3) barycentres of faces and (4) barycentres of polyhedrons

/// Metricised PCC-related values:: (metric is needed)
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector; // Global vectors of measures: edge lengths, face areas and polyhedra volumes

/// * ===================== MODULES and LIBRARIES ==============================* ///
///* =========================================================================== *///

/*! Various supplementary useful functions are defined here */
#include "lib/processing_design_lib/PCC_Support_Functions.h" // It must be here - first in this list of libraries (!)

/*! An Objects library contains classes of various objects related to the PCC's substructures and k-cells */
#include "lib/processing_design_lib/PCC_Objects.h"

/*! Various set combinatorial measures are defined here */
#include "lib/processing_design_lib/PCC_Measures.h"

// * Each 'Module' have its own *.ini file for tailored input and *.log for tailored output defined by user */
// * There are only 2 'Principal' modules -- Processing and Design -- creating 'processing_vector' and 'design_vector' output containing 'history' of all the 'pcc_structure' changes  */
// ----------------------------------------------------------------------------------------------------------
// * The numeration order of the Modules is strict throughout the CPD code:
    //|   I.    Main               (technical)
    //|   II.   Subcomplex         (supplementary)
    //|   III.  Multiphysics       (supplementary)
    //|   IV.   Processing         (principal)
    //|   V.    Kinetics           (principal)
    //|   VI.   Characterisation  (supplementary)
    //|   VII.  Design            (principal)
    //|   VIII. Writer            (supplementary)

/*! SUBCOMPLEX module calculates reduced PCC subcomplexes as parts of the initial PCC, inheriting (reduced) sequences of special cells and 'state vectors' of the original PCC */
/* Supplementary project module */
#include "lib/processing_design_lib/PCC_Subcomplex/PCC_Subcomplex.h"

/*! MULTIPHYSICS module set self, elastic and thermal energies with any energy-related values associated with all k-cells in a PCC taking data from the "CPD_material_database" and config/multiphysics.ini files */
/* Supplementary project module */
#include "lib/processing_design_lib/PCC_Multiphysics/PCC_Multiphysics.h"

/*! PROCESSING module assigned special IDs (labels, colours) for various PCC cells -- Nodes, Edges, Faces, Polytopes/Polyhedrons in the corresponding tessellation of space */
/* Principal project module */
#include "lib/processing_design_lib/PCC_Processing/PCC_Processing.h"
/* Output: module generates a 'processed_structure_sequences' as the lists containing the sequences of k-cells possessing 'special' IDs (labels, colours), including
 * (1) ASSIGNED: k-Cells, k={0,1,2,3}, corresponding to different generation principles (random, maximum entropy,.. etc.),
 * (2) INDUCED: m-Cells (where m < k) directly labelled based on the HIGHER-ORDER k-Cell IDs (labels)
 * (3) GENERATED: g-Cells generated as a result of some KINETIC process. They always depend on the ASSIGNED design 's_cell_sequences' vectors of 'special' k-Cells and/or INDUCED vectors of 's_generated_cell_sequences' */

/*! KINETICS module implements the structural changes (affecting state Vectors) appearing in 'time'
/* Principal project module */
#include "lib/processing_design_lib/PCC_Kinetics/PCC_Kinetics.h"
// Provides for each p-cell in a PCC the moment of 'time' in the [0,1] range when it changed its special 'generated' type because of a 'kinetic' process (e.g. corrosion or irradiation)

/*! DESIGN module implements the optimisation of the assigned and induced State Vectors generated by the Processing module according to some 'goal' functions */
/* Principal project module */
#include "lib/processing_design_lib/PCC_Design/PCC_Design.h"
// Provides vector (or list) of the 'designed_structure_sequences' as the lists of structures formed of the state vectors, containing 'special' and/or 'induced' types) in order of their 'evolution' towards an 'optimum' configuration.

/*! CHARACTERISATION module provides vectors with the characteristics (entropic, spectral, etc.) representing evolution of the PCC state vectors (PCC 'STATE') as provided by the output of the PCC_Processing module */
/* Supplementary project module */
#include "lib/processing_design_lib/PCC_Characterisation/PCC_Characterisation.h"
// Calculates various combinatorial measures associated with various PCC structures and networks

/*! WRITER module performs formatted output of various data structures generated by the PCC_Processing (PCC state vectors) and PCC_Design (PCC design vectors) modules */
/* Supplementary project module */
#include "lib/processing_design_lib/PCC_Writer/PCC_Writer.h"
// Provides formatted output of the sequences of the (1) state vectors, (2) structures made of them, and (3) structural characteristics both -- for the 'processed_structure_sequences' and the 'designed_structure_sequences'

/// PCC special structure-related variables :: see class 'Config' in Objects.h and Object.cpp files
Config initial_configuration, configuration; // Configuration is an object of the 'Config' class described in the Objects.cpp project library

/*!
 * @brief TUTORIAL :: An educational course active in the 'TUTORIAL' execution mode of the CPD code; please see 'config/main.ini' file for more details.
 * @param initial_configuration
 */
void tutorial(Config &initial_configuration);

/*!
 * @brief PERFORMANCE_TEST :: A global testing mode of the CPD code; please see 'config/main.ini' file for more details.
 * @param initial_configuration
 */
void performance_test(Config &initial_configuration);

///* ........................................................................................    Main    ................................................................ *///
//* (.h files) * @brief, @param and @return
//* (.cpp files) * @details (detailed descriptions)
/*!
 * @details MAIN module of the CPD code implements the whole program execution according to the specifications written in the 'config/*.ini' files.
 * In particular, the LIST mode calls for the execution of the project modules one by one (if active) and computes the execution time for each of them.
 * The TASK mode is a flexible research tool, allowing the use of tailored user code (by #include < *.h >), using libraries, classes and functions included in the CPD code.
 * The TUTORIAL and PERFORMANCE_TEST modes implement the educational and hardware testing mode essential for the first execution of the code on a new computation cluster.  * @param void
 * @return 0 and the output to console and the *.log files, if successful.
*/
int main() {
// * ------------------ #Print ------------------------------------
    cout << endl << "---------------------------------- *** CPD code execution begins *** --------------------------------------------------------------" << endl << endl;
// * --------------------------------------------------------------
// Output Year/Day/Time of the computation to cout
    time_t timestamp = time(&timestamp);
    struct tm datetime = *localtime(&timestamp);
    cout << " Execution year - " << 1900 + datetime.tm_year << "; Month - " << datetime.tm_mon << "; Day - " << datetime.tm_mday << "; Time - " << datetime.tm_hour << ":" << datetime.tm_min << "." << endl << endl;

/// ========== Elapsing time of the MAIN module =========== ///
    Main_execution_time = (double) clock();

/// Initial configuration reader and information output to the screen and into the 'Processing_Design.log' file
    initial_configuration.Read_config(initial_configuration); // Read_config() method of the class Config defined in Objects.h and described in Objects.cpp

/// Setting values of the global variables:: all the methods below are in the class Config defined in Objects.h and described in Objects.cpp
    std::vector<int> ConfigVector = initial_configuration.Get_ConfVector();
    PCC_dimension = initial_configuration.Get_dim();
    source_dir = initial_configuration.Get_source_dir();
    output_dir = initial_configuration.Get_output_dir();
    paths_to_PCC_matrices = initial_configuration.Get_paths();
    main_type = initial_configuration.Get_main_type();
    simulation_tasks_dir = initial_configuration.Get_sim_task();
 /// ============================================================================== ///

// ------------------ #Print -----------------------
 //   main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::trunc); // the main_logfile_stream.log stream will be closed at the end of the main function
    main_logfile_stream << endl << "---------------------------------- *** CPD code execution begins *** --------------------------------------------------------------" << endl << endl;

    // Output Year/Day/Time of the computation
    main_logfile_stream << " Execution year - " << 1900 + datetime.tm_year << "; Month - " << datetime.tm_mon << "; Day - " << datetime.tm_mday << "; Time - " << datetime.tm_hour << ":" << datetime.tm_min << "." << endl << endl;
    {
        std::string print_to_string = "Main execution time before modules is equal to  "s + std::to_string(Main_execution_time / pow(10.0, 6.0)) + "  seconds"s;
        cout << print_to_string << endl << endl;
        main_logfile_stream << print_to_string << endl << endl;
    }
// -------------------------------------------------

/// ========================================================================================================================================== ///
/// ================================================= PERFORMANCE_TEST MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    if (main_type == "PERFORMANCE_TEST"s) { // testing mode 'PERFORMANCE_TEST' which should be further replaced in the 'config/main.ini' file with the 'PERFORMANCE_TEST' and then the 'TASK' or the 'LIST' modes
        cout << "==================================================================================================================================================" << endl;
        main_logfile_stream << "==============================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        main_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        performance_test(initial_configuration); // output the 'performance_test.txt' file to the 'output_dir' showing the relative code execution times of the present server comparing with some reference execution times and suggest the preferable PCC sizes for various simulation tasks
    } /// END of the SIMULATION MODE "PERFORMANCE_TEST" as specified in the config/main.ini file

/// ========================================================================================================================================== ///
/// ================================================= TUTORIAL MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    else if (main_type == "TUTORIAL"s) { // TUTORIAL feature to facilitate the first acquaintance with the code: only in the 'TUTORIAL' execution type = the 'mode' variable in the config/main.ini file.
        cout << "==================================================================================================================================================" << endl;
        main_logfile_stream << "==============================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        main_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

        tutorial(initial_configuration);
    } /// END of the SIMULATION MODE "TUTORIAL" as specified in the config/main.ini file

/// ========================================================================================================================================== ///
/// ================================================= LIST MODE STARTS HERE ================================================================== ///
/// ========================================================================================================================================== ///
    else if ( main_type == "LIST"s ) { // In the LIST mode all the functions are calling one after another without additional loops and intermediate data output
    /// For all the more complicated simulation cases the TASK mode should be used - see it following next after the 'LIST' module.
        cout << "==================================================================================================================================================" << endl;
        main_logfile_stream << "=======================================================================================================================================================================================================================================" << endl;
        cout << "\t\t\t\t\t\t\t\t\t\t[\tStart of the PCC Processing Design code \t]\t\t\t\t\t\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
        main_logfile_stream << "\t\t\t\t\t[\tStart of the PCC Processing Design code\t]\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

    /// Initialisation of the current_configuration as equal to the initial_configuration
        configuration = initial_configuration;

    /// ====================== I. PCC Subcomplex module ======================
        std::vector<Subcomplex> pcc_subcomplexes; // vector containing all the PCC subcomplexes (cuts, k-order grain neighbours, etc)

        if (ConfigVector.at(1) == 1) { // if the 'PCC_Section' parameter is switched 'ON' in the config/main.ini file
   //         main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << " START of the PCC Subcomplex module " << endl;
            main_logfile_stream << " START of the PCC Subcomplex module " << endl;

///            pcc_subcomplexes = PCC_Subcomplex(configuration);

            cout << " pcc_subcomplexes size =  " << pcc_subcomplexes.size() << endl;

            // ================ Elapsing time for the Subcomplex module ================
            unsigned int Subcomplex_time = clock();
            Subcomplex_execution_time = (double) Subcomplex_time - Main_execution_time;
            cout << "Section time is equal to  " << Subcomplex_execution_time / pow(10.0, 6.0) << "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
            main_logfile_stream << "Section time is equal to  " << Subcomplex_execution_time / pow(10.0, 6.0) << "  seconds" << endl;
            main_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(SectionON)

        /// ====================== II. PCC Multiphysics module ======================
        std::vector<CellEnergies> new_cells_energies; // a class described in PCC_Objects.h contained (1) all the k-cell elastic energies and (2) all the k-cell thermal energies in the PCC
        // Example: vector<CellEnergies> for several crack lengths in a PCC

        if (ConfigVector.at(2) == 1) { // if the 'PCC_Multiphysics' parameter is switched 'ON' in the config/main.ini file
      //      main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Multiphysics module " << endl << endl;
            main_logfile_stream << "START of the PCC Multiphysics module " << endl << endl;

            /// Defects
            std::vector<Macrocrack> crack_growth_series; // series of objects of the class Macrocrack with different lengths simulating a crack growth

///            new_cells_energies = PCC_Multiphysics(configuration, pcc_subcomplexes, crack_growth_series);

            // ================ Elapsing time for the Processing module ================
            unsigned int Multiphysics_time = clock();
            Multiphysics_execution_time = (double) Multiphysics_time - Subcomplex_execution_time - Main_execution_time;
            cout << endl << "Multiphysics time is equal to  " << Multiphysics_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << endl << "Multiphysics time is equal to  " << Multiphysics_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(MultiphysicsON)

        /// ====================== III. PCC Processing module ======================
        CellDesign new_cells_design; // a class described in PCC_Objects.h contained (1) all special k-cell sequences and (2) all the design_<*>_vectors for all k-cells in the PCC

        if (ConfigVector.at(3) == 1) { // if the 'PCC_Processing' parameter is switched 'ON' in the config/main.ini file
      //      main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Processing module " << endl;
            main_logfile_stream << "START of the PCC Processing module " << endl;

            new_cells_design = PCC_Processing(configuration);

        // ================ Elapsing time for the Processing module ================
            unsigned int Processing_time = clock();
            Processing_execution_time = (double) Processing_time - Subcomplex_execution_time - Multiphysics_execution_time - Main_execution_time;
            cout << "Processing time is equal to  " << Processing_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "Processing time is equal to  " << Processing_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(ProcessingON)

        /// ====================== IV. PCC Kinetics module ======================
        std::vector<vector<double>> p_cells_history;
        /// (CellNumbs.at(0),CellNumbs.at(1),CellNumbs.at(2),CellNumbs.at(3));

        if (ConfigVector.at(7) == 1) { // if the 'PCC_Kinetics' parameter is switched 'ON' in the config/main.ini file
       //     main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Kinetics module " << endl;
            main_logfile_stream << "START of the PCC Kinetics module " << endl;

            // TODO: TEMPORARY MODULE OUTPUT
            corrosion_damaged_output.open(output_dir + "surface_corrosion_damaged_output.txt"s, ios::trunc);
            corrosion_affected_output.open(output_dir + "surface_corrosion_affected_output.txt"s, ios::trunc);

            face_barycentre_coord_outstream.open(output_dir + "face_seeds.txt"s, ios::trunc);
            edge_barycentre_coord_outstream.open(output_dir + "edge_seeds.txt"s, ios::trunc);

            corrosion_damaged_fractions_output.open(output_dir + "area_corrosive_damaged_fraction.txt"s, ios::trunc);
            corrosion_affected_fractions_output.open(output_dir + "area_corrosive_affected_fraction.txt"s, ios::trunc);

            p_cells_history = PCC_Kinetics(configuration, new_cells_design);

            corrosion_damaged_output.close();
            corrosion_affected_output.close();

            edge_barycentre_coord_outstream.close();
            face_barycentre_coord_outstream.close();

            corrosion_damaged_output.close();
            corrosion_affected_fractions_output.close();

            // ================ Elapsing time for the Kinetics module ================
            unsigned int Kinetics_time = clock();
            Kinetics_execution_time = (double) Kinetics_time - Processing_execution_time - Subcomplex_execution_time - Multiphysics_execution_time - Main_execution_time;
            cout << "Kinetics time is equal to  " << Kinetics_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "Kinetics time is equal to  " << Kinetics_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(KineticsON)

            /// ====================== III. PCC Characterisation module ======================
        ProcessedComplex pcc_processed;  // a class described in PCC_Objects.h

        if (ConfigVector.at(4) == 1) { // if the 'PCC_Characterisation' parameter is switched 'ON' in the config/main.ini file
     //       main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Characterisation module" << endl; main_logfile_stream << "START of the PCC Characterisation module" << endl;
            cout << "=========================================================================" << endl;
            main_logfile_stream << "==============================================================================================================================================================" << endl;

///            pcc_processed = PCC_StructureCharacterisation(new_cells_design);

        // ===== Elapsing time for the Characterisation module ================
            unsigned int Characterisation_time = clock();
            Characterisation_execution_time = (double) Characterisation_time - Subcomplex_execution_time - Multiphysics_execution_time - Processing_execution_time - Main_execution_time;
            cout << "Characterisation time is equal to  " << Characterisation_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "Characterisation time is equal to  " << Characterisation_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        }// end if(CharacterisationON)

        /// ====================== IV. PCC Design module ======================
        if (ConfigVector.at(5) == 1) { // if the 'PCC_Design' parameter is switched 'ON' in the config/main.ini file
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;

            cout << "START of the PCC Design module" << endl;
            main_logfile_stream << "START of the PCC Design module" << endl;
            cout << "=========================================================================" << endl;
            main_logfile_stream << "==============================================================================================================================================================" << endl;

            std::vector<std::vector<int>> pcc_design;
            pcc_design = PCC_Design(configuration, new_cells_design);

            // ===== Elapsing time for the PCC Design module ================
            unsigned int Design_time = clock();
            Design_execution_time = (double) Design_time - Subcomplex_execution_time - Main_execution_time - Subcomplex_execution_time - Multiphysics_execution_time - Processing_execution_time - Characterisation_execution_time;
            cout << "Design time is equal to  " << Characterisation_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "Design time is equal to  " << Characterisation_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl; //Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(DesignON)

        /// ====================== V. PCC Writer module ======================
        if (ConfigVector.at(6) == 1) { // if the 'PCC_Writer' parameter is switched 'ON' in the config/main.ini file
            main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
            cout << "-------------------------------------------------------------------------" << endl;
            main_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Writer module" << endl;
            main_logfile_stream << "START of the PCC Writer module" << endl;
            cout << "=========================================================================" << endl;
            main_logfile_stream << "==============================================================================================================================================================" << endl;

// alternative:           PCC_Writer(new_cells_energies, new_cells_design, pcc_processed);
///            PCC_Writer(new_cells_design);

        // ================ Elapsing time for the Writer module ================
            unsigned int Writer_time = clock();
            Writer_execution_time = (double) Writer_time - Main_execution_time - Subcomplex_execution_time - Multiphysics_execution_time - Processing_execution_time - Characterisation_execution_time - Design_execution_time;
            cout << "Writer time is equal to  " << Writer_execution_time / pow(10.0, 6.0) << "  seconds" << endl << endl;
        } // end if(WriterON)

    } /// END of the SIMULATION MODE "LIST" as specified in the config/main.ini file


/// ==========================================================================================================================================
/// ================================================= TASK MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
    else if ( main_type == "TASK"s ) { // In the TASK mode any piece of code using the project libraries can be included

/*
#include "tasks/energy_levels.h"
        energy_plasticity();
exit(0);
*/

/**
        string indir1 = "/Users/user/Dropbox/OFFICE/Communications/Siying/Belgorod_data_2022/Experimental_Data/ECAP2022/texture_effect/Neper_Cu_01Cr_01Zr_pass4_rand/rand1/pass4_xy.txt"s;
        string indir2 = "/Users/user/Dropbox/OFFICE/Communications/Siying/Belgorod_data_2022/Experimental_Data/ECAP2022/texture_effect/Neper_Cu_01Cr_01Zr_pass4_rand/rand2/pass4_xy.txt"s;
        string outdir1 = "/Users/user/Dropbox/OFFICE/Communications/Siying/Belgorod_data_2022/Experimental_Data/ECAP2022/texture_effect/Neper_Cu_01Cr_01Zr_pass4_rand/rand1/"s;
        string outdir2 = "/Users/user/Dropbox/OFFICE/Communications/Siying/Belgorod_data_2022/Experimental_Data/ECAP2022/texture_effect/Neper_Cu_01Cr_01Zr_pass4_rand/rand2/"s;
        std::string f1 = outdir1 + "frand_1.txt"s, f2 = outdir2 + "frand_2.txt"s;
        shuffle_coord_vector(indir1, f1, 10283);
        shuffle_coord_vector(indir2, f2, 10283);
**/

/////////////
/// #include "tasks/macrocrack_growth.h"
////////////

        /// Initialisation of the current_configuration as equal to the initial_configuration
        configuration = initial_configuration;

        main_logfile_stream.open(output_dir + "cpdlog_main.log"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function

        vector<unsigned int> node_coordinates_seq, face_coordinates_seq, polytope_coordinates_seq;
        //      for(unsigned int i = 0; i < CellNumbs.at(0); ++i) { node_coordinates_seq.push_back(i); }
        //       node_coordinates_vector = kCell_barycentre_coordinates(0,node_coordinates_seq);

        for (unsigned int i = 0; i < CellNumbs.at(2); ++i) {
            face_coordinates_seq.push_back(i);
        }

///        face_coordinates_vector = kCell_barycentre_coordinates(2, face_coordinates_seq);
        face_coordinates_vector = Tuple3Reader(paths_to_PCC_matrices.at(13));
        if(face_coordinates_vector.size() == 0)
            face_coordinates_vector = kSequence_barycentre_coordinates(2, face_coordinates_seq);

        //     for(unsigned int i = 0; i < CellNumbs.at(3); ++i) { polytope_coordinates_seq.push_back(i); }
        //     polytope_coordinates_vector = kCell_barycentre_coordinates(3, polytope_coordinates_seq);

        /// Function simulating Microcrach Growth
        ProcessedComplex cracked_pcc;
        std::vector<Macrocrack> macrocrack_growth_series; // series of objects of the class Macrocrack with different lengths simulating a crack growth
        std::ofstream agglomeration_stat_out;
        agglomeration_stat_out.open(output_dir + "Agglomeration_stats.txt"s, ios::trunc);
        agglomeration_stat_out << " powder_number " << "\t" << " max_power " << "\t" << " aggls_in_powder.size() " << "\t" << " counter " << "\n"; // for each POWDER
        agglomeration_stat_out.close();

        unsigned int counter = 0;
/// DO        do{ //while loop

        agglomeration_stat_out.open(output_dir + "Agglomeration_stats.txt"s, ios::app);

///        cracked_pcc = Macrocrack_growth(configuration); //, macrocrack_growth_series);

        ///  cout << "SIZE\t " << cracked_pcc.Get_macrocrack_sfaces_series().at(0).size() << endl; exit(0);
        std::ofstream Cracked_pcc_out, Cracked_betti_pcc_out;
        cout << "-------------------------------------------------------------------------" << endl;
        main_logfile_stream << "-------------------------------------------------------------------------" << endl;
        cout << "START of the Writer " << endl;
        main_logfile_stream << "START of the Writer " << endl;
        cout << "=========================================================================" << endl;
        main_logfile_stream
                << "=============================================================================================================================================================="
                << endl;

        //std::vector<std::vector<unsigned int>> macrocrack_faces_set = cracked_pcc.Get_macrocrack_sfaces();
        // std::vector<std::vector<std::tuple<double, double, double>>> msb_coord_vector;

        /// Finding inclusions in the current powder
//
//        current_PCC.Set_face_barycentre_coordinates();

/*
 *         Cracked_pcc_out.open(output_dir + "Macrocrack_sfaces_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function

        PCC current_PCC;
        std::vector<std::tuple<double, double, double>> FBC = current_PCC.Get_face_barycentre_coordinates();
        std::vector<std::vector<std::tuple<double, double, double>>> vector_of_series_of_sface_coord_tuples;
        for ( Macrocrack mcrack : macrocrack_growth_series) {
            std::vector<unsigned int> new_sface_sequence = mcrack.Get_sfaces_sequence();
            vector_of_series_of_sface_coord_tuples.push_back(face_sequence_barycentre_coordinates(new_sface_sequence, FBC));
            cout << "msb SIZE\t" << vector_of_series_of_sface_coord_tuples.back().size() << endl;
        } // end pf for (auto sf_set: microcrack_faces_set)

            for (auto acc : vector_of_series_of_sface_coord_tuples.back()) {
                cout << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << endl;;
                Cracked_pcc_out << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << endl;;
            }
//            cout << endl; Cracked_pcc_out << endl;

        Cracked_pcc_out.close();
        exit(0);
*/
        Cracked_pcc_out.open(output_dir + "Macrocrack_D_sFace_fractions.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
        int powder_iterator = 2;
        for (auto dfrac: cracked_pcc.de_fractions_sface_vector) {
            Cracked_pcc_out << "Powder #\t" << powder_iterator++ << endl;
            for (auto itdf: dfrac) {
// REPAIR                cout << itdf << "\t";
                Cracked_pcc_out << itdf << "\t";
            }
            cout << endl; Cracked_pcc_out << endl;
        }
        Cracked_pcc_out.close();

            Cracked_betti_pcc_out.open(output_dir + "Macrocrack_sFace_Betti.txt"s,ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
            powder_iterator = 0;
            for (int betti_itr = 0; betti_itr < cracked_pcc.Betti_0_sface.size(); ++betti_itr) {
                Cracked_betti_pcc_out << "Powder #\t" << powder_iterator++ << endl;
                Cracked_betti_pcc_out << cracked_pcc.Betti_0_sface.at(betti_itr) << "\t"  << cracked_pcc.Betti_1_sface.at(betti_itr) << "\t" << cracked_pcc.Betti_2_sface.at(betti_itr) << "\t" << cracked_pcc.inverse_connectivity_sface.at(betti_itr) << endl;
            }
            Cracked_betti_pcc_out.close();

            Cracked_betti_pcc_out.open(output_dir + "Macrocrack_iFace_Betti.txt"s,ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
            powder_iterator = 0;
            for (int betti_itr = 0; betti_itr < cracked_pcc.Betti_0_iface.size(); ++betti_itr) {
                Cracked_betti_pcc_out << "Powder #\t" << powder_iterator++ << endl;
                Cracked_betti_pcc_out << cracked_pcc.Betti_0_iface.at(betti_itr) << "\t"  << cracked_pcc.Betti_1_iface.at(betti_itr) << "\t" << cracked_pcc.Betti_2_iface.at(betti_itr) << "\t" << cracked_pcc.inverse_connectivity_iface.at(betti_itr) << endl;
            }
            Cracked_betti_pcc_out.close();

            Cracked_pcc_out.open(output_dir + "Macrocrack_agglomeration_powers.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function

        int numb = 0, npowd = 0;
        double average_power = 0.0;
        for (auto aggls_in_powder: cracked_pcc.agglomerations_in_powders) {
            ++npowd;
            average_power = 0;
            numb = 0;
            for (auto aggl: aggls_in_powder) {
                ++numb;
                Cracked_pcc_out << aggl.Get_agglomeration_power() << "\t";
                average_power += aggl.Get_agglomeration_power();
            }

            cout << " Powder #:\t" << npowd << " agglomeration fraction:\t" << (double) numb / CellNumbs.at(2) << "\t"
                 << " agglomeration relative fraction:\t" << (double) numb / CellNumbs.at(2) << "\t"
                 << "average agglomeration power:\t" << average_power / (double) numb << "\t" << endl;
            //Cracked_pcc_out << aggl.Get_agglomeration_power() << "\t";
            cout << endl;
            Cracked_pcc_out << endl;
        }
        Cracked_pcc_out.close();

        double max_power = 0.0;
        int powder_number = 0;
        for (auto aggls_in_powder : cracked_pcc.agglomerations_in_powders) {
            ++powder_number;

            max_power = 0;

            cout << "powder number\t" << powder_number << "\t\t\t";
            for (auto aggl : aggls_in_powder) {
                cout << aggl.Get_agglomeration_power() << "\t";
                if (aggl.Get_agglomeration_power() > max_power)
                    max_power = aggl.Get_agglomeration_power();
            }
            cout << endl;
            agglomeration_stat_out << powder_number << "\t" << max_power << "\t" << aggls_in_powder.size() << "\t" << counter << endl; // for each POWDER
        }
        agglomeration_stat_out << endl;
        Cracked_pcc_out.open(output_dir + "Macrocrack_agglomeration_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
        std::vector<unsigned int> aggl_cells_sequence;
        aggl_cells_sequence.clear();

        int powder_numb = 0; /// A SPECIFIC number of a powder for the output of the agglomeration coordinates
        ///============================
        for (auto aggl: cracked_pcc.agglomerations_in_powders.at(powder_numb)) {
            aggl_cells_sequence.push_back(aggl.Get_agglomeration_kcell_number());
        }

        std::vector<std::tuple<double, double, double>> aggl_cells_coords;
        aggl_cells_coords.clear();
        aggl_cells_coords = face_sequence_barycentre_coordinates(aggl_cells_sequence);

        unsigned int a_numerator = 0;
        for (auto acc: aggl_cells_coords) {
// REPAIR            cout << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << cracked_pcc.agglomerations_in_powders[powder_numb][a_numerator].Get_agglomeration_power() << endl;
            Cracked_pcc_out << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t"
                            << cracked_pcc.agglomerations_in_powders[powder_numb][a_numerator].Get_agglomeration_power()
                            << endl;
            a_numerator++;
        }

        Cracked_pcc_out.close();
///        main_logfile_stream.close();

        agglomeration_stat_out.close();

/// WHILE } while (++counter < 1000);

    } /// END of the SIMULATION MODE "TASK" as specified in the config/main.ini file
/// ==========================================================================================================================================
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl << "\t\t\t\t\t\t\t\t\t\t[\tThe end of the PCC Processing Design\t]\t\t\t\t\t\t\t\t\t\t" << endl << "==================================================================================================================================================" << endl;
    main_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl << "\t\t\t\t\t\t\t\t\t\t[\tThe end of the PCC Processing Design\t]\t\t\t\t\t\t\t\t\t\t" << endl << "==================================================================================================================================================" << endl;

/// ================ Total CPD code Elapsing time ================ ///
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "Total " << PCC_dimension << "D " << "runtime of the CPD code is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    main_logfile_stream << "Total " << PCC_dimension << "D " << "runtime of the CPD code is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    main_logfile_stream << "-------------------------------------------------------------------------" << endl;

    // closing all project off-streams
    main_logfile_stream.close();
    subcomplex_logfile_stream.close();
    multiphysics_logfile_stream.close();
    processing_logfile_stream.close();
    characterisation_logfile_stream.close();
    design_logfile_stream.close();
    writer_logfile_stream.close();

    return 0; // success!
} /// The END of the Main function

/// ========================================= FUNCTIONS DEFINED IN MAIN MODULE =======================================///

/// ====================# 1 #========================== TUTORIAL ==================================================== ///
/*!
 * @details The TUTORIAL code execution mode is designed to make the first acquaintance with the code easier. It is simply a tour around
 * the typical output of the code finishing with the discussion of its user interface ('config/_.ini' files) and of how to
 * replace the TUTORIAL execution mode itself with the PERFORMANCE_TEST and then LIST or TASK execution modes.
 * @param initial_configuration
 */
void tutorial(Config &initial_configuration){ // teaching function
    cout << "Hello there! This is the TEACHING mode of the program execution (!), \n"
    "where the user passes tutorial helping they to understand how to work with the CPD code. \n"
    "It can be changed by replacing the 'mode' variable with 'list' or 'task' \n"
    "instead of the current 'tutorial' contained in the 'simulation_mode' in ../config/main.ini file (!)" << endl;
    cout << endl;
    cout << "Continue? Please type 'Y' for 'yes' or 'N' for 'no' " << endl;
    char if_continue = 'Y'; cin >> if_continue;
    if (if_continue == 'N') exit(0);
} /// END of the tutorial() function

/// ====================# 2 #====================== PERFORMANCE TEST ================================================ ///
/*!
 * @details The PERFORMANCE TEST code execution mode is aimed to allow user estimate the size (number of cells) of the PCCs which they can employ in calculation for different types of the Processing, Characterisation and Design tasks.
 * It writes the analysis data to the 'performance_test.txt' file contained in the 'output_dir' directory showing the relative code execution times of the present computer comparing to some reference execution times.
 * So the user can decide which computations can be done using the current facility for the reasonable time.
 * @param initial_configuration
 */
void performance_test(Config &initial_configuration) { // execution test function

    std::ofstream Out_performance_test_stream; // performance_test output stream to the output file

    /// Initialisation of the current_configuration = initial_configuration
    configuration = initial_configuration;

    // times
    double processing_execution_time = 0.0, full_processing_time = 0.0;
    unsigned int prev_time, new_time;

    int counter_max = 50000; //number of calculation series

    Out_performance_test_stream.open(output_dir + "cpd_code_performance_test.txt"s, ios::trunc); // creates 'performance_test.txt' file in the 'output_dir'
    Out_performance_test_stream.close();

    Out_performance_test_stream.open(output_dir + "cpd_code_performance_test.txt"s, ios::app); // open for writing 'performance_test.txt' file in the 'output_dir'
    Out_performance_test_stream << "------------------------------------------------------------------------------------------------" << endl;

    CellDesign new_cells_design;
    for (int counter = 0; counter < counter_max; counter++ ) { // TEST LOOP
        prev_time = clock();
        std::vector<int> ConfigVector = configuration.Get_ConfVector();
        PCC_dimension = configuration.Get_dim();
        source_dir = configuration.Get_source_dir();
        output_dir = configuration.Get_output_dir();
        paths_to_PCC_matrices = configuration.Get_paths();
        /// ---------------------------------------------------------------------- ///
///        new_cells_design = PCC_Processing(configuration);

        /// ============= Elapsing time Processing ================ ///
        new_time = clock();
        processing_execution_time = (double) new_time - (double) prev_time;
        full_processing_time += processing_execution_time;
        //prev_time = processing_execution_time;

        // Output to the 'performance_test.txt' file in the 'output_dir
        Out_performance_test_stream << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl << endl;

        cout << endl << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl;
        Out_performance_test_stream << endl << "Processing iteration " << counter + 1 << " tooks  " << processing_execution_time/ pow(10.0,6.0) <<  "  seconds" << endl << "-------------------------------------------------------------------------" << endl;

    } // end of for (int counter = 0; counter < counter_max; counter++ ) loop

    // Output to the 'performance_test.txt' file in the 'output_dir
    Out_performance_test_stream << endl << "Full execution time for " <<  counter_max << " iterations is equal to  " << full_processing_time/ pow(10.0,6.0) <<  "  seconds" << endl;

    cout << "Full execution time for " <<  counter_max << " iterations is equal to  " << full_processing_time/ pow(10.0,6.0) <<  "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl;
    main_logfile_stream << "Processing time is equal to  " << Processing_execution_time / pow(10.0, 6.0) << "  seconds" << endl;
    main_logfile_stream << "-------------------------------------------------------------------------" << endl;

    Out_performance_test_stream.close(); // close output to the 'performance_test.txt' file

} /// END of the performance_test() function



                                            /// *** H E A P *** ///