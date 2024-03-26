///******************************************************************************************************************************///
///************************   Polytopal Cell Complex (PCC) Processing Design :: (CPD code) (c)   *******************************///
///****************************************************************************************************************************///
///*                                        Version 4.0 | 12/02/2024                                                         *///
///**************************************************************************************************************************///
///************************************ Dr Elijah Borodin, Manchester, UK **************************************************///
///**************************************** Spring 2022 - Spring 2024  ****************************************************///
///***********************************************************************************************************************///
///*
///*    Code source:    https://github.com/PRISBteam/PCC_Processing_Design/
///*    Documentation:  https://prisbteam.github.io/
///*    PCC sources:    https://materia.team/
///*
///*  The project provides a reliable tool for 1. Obtaining, 2. Analysing and 3. Improving of 'design vectors' as the sequences                          *///
///*  of k-cells containing in the k-skeletons, where k = {0,1,2,3}, of a Polytopal Cell Complex (PCC). Such PCCs can be created by external             *///
///*  codes based on the tessellation of 2D or 3D spaces by an agglomeration of polytopes (polygons in the 2D case or polyhedrons in 3D).                *///
///*  Graphs and networks (without loops) are considered as 1-complexes (1-PCCs) and also available for analysis similarly to the the 2D and 3D cases.   *///

///* Key terminology:
/// Material's elements         ::   'quadruple points', 'grain boundary junctions', 'grain boundaries', and 'grains' (with their orientations and barycenter coordinates)
/// Tessellation's elements     ::   'nodes, 'edges', 'faces', 'polytopes' (with their measures - lengths, areas and volumes - and barycenter coordinates)
/// PCC's elements              ::   'k-cells' containing in 'k-skeletons', where k = {0,1,2,3}, with their degree fractions, and incident (k-1)-cells and (k+1)-cells.

///* ----------------------------------------- *
///* Standard C++ (STL) libraries
///* ----------------------------------------- *
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include <cmath>
// #include <numeric>
// #include <algorithm>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2024)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

/// Spectra source: https://spectralib.org/ (2024)
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
// Included only in the parallelized version of the code

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project based on the mINI C++ library (2018 Danijel Durakovic http://pulzed.com/, MIT license is included)
#include "lib/ini/ini_readers.h"

///------------------------------------------
using namespace std; // STL namespace
using namespace Eigen; // Eigen library namespace
using namespace Spectra; // Spectra library namespace

/// Eigen library-based class Triplets class containing objects in the form T = T(i, j, value), where i and j are element's a(i, j) indices in the corresponding dense matrix and the third variable is its value
typedef Triplet<double> Tr; // <Eigen> library class, which declares a triplet type with the nickname 'Tr'
typedef SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'
typedef MatrixXd DMat; // <Eigen> library class, which declares a dense matrix type of doubles with the nickname 'DMat'

/// * ---------------------------------------------------------------------------------------------------------- *///
/// * ============================================ GLOBAL VARIABLES ============================================ *///
/// * ----------- Declaration of GLOBAL variables (can be seen in all the project modules and libraries)-------- *///

/// Technical variables
string source_path = "../config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // 'source_path' is a path to the directory reading from the 'config/main.ini' file
string e_mode = "tutorial"; // '[execution_type]' (etype variable) from the config/main.ini file. The "tutorial" mode by default means the 'educational' mode, where the program learn of how to work with its input files and simulation types
std::string main_type; /// SIMULATION MODE(LIST) or 'SIMULATION MODE(TASK) :: This define the global simulation mode: 'LIST' for the "list" of modules implementing one by one (if ON) and 'TASK' for the user-defined task scripts with modules and functions included from the project's libraries

vector<char*> PCCpaths; // The vector containing the PCCpaths to all the PCC's matrices, measures and other supplementary data
string source_dir, output_dir; // Input and output directories as it is written in the 'config/main.ini' file
string sim_task; // path to the corresponding *.cpp file containing a 'simulation task' (for 'TASK' execution mode only, not 'LIST') as it is written in the 'config/main.ini' file

/// Global 'log.txt' file output
std::ofstream Out_logfile_stream; // 'log.txt' file output of the entire computation process as a copy of the console output

/// PCC - related variables
// [general]
int dim; // PCC's dimension: dim = 1 for graphs and networks, dim = 2 for the 2D plane polygonal complexes, and dim = 3 for the 3D bulk polytopal complexes, as it is specified in the 'config/main.ini' file.

// Combinatorial
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs containing the numbers of k-cells of different types 'k'. It is read from the 'number_of_cells.txt' file.
// First line here is the number of nodes (0-cells), second - edges (1-cells), third - faces (2-cells) (in the 2D and 3D cases only), fourth - polyhedra (3-cells) (in the 3D case only)

// Geometry
vector<tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // vectors containing barycenter Cartesian coordinates of the corresponding tessellation's elements
// Global vectors of Cartesian coordinates for: (1) vertex coordinates, (2) barycentres of edges, (3) barycentres of faces and (4) barycentres of polyhedrons

// Measures
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector; // Global vectors of measures: edge lengths, face areas and polyhedra volumes

/// PCC special structure-related variables
// State_Vector in the form : [Element index] - > [Element type], like [0, 0, 2, 1, 1, 0, 2, 4, 3, 3, 2, 0,... ...,2] containing all CellNumb.at(*) element types
std::vector<int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // Normally the State_<*>_vector of special cells can be calculated based on the corresonding special_cell_sequences
std::vector<int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // separate vectors containing the other 'fractured' labels different from the 'special' ones. To be calculated based on the corresonding fractured_cell_sequences
/// Configuration_State = { State_p_vector, State_f_vector, State_e_vector, State_n_vector } is a list of all 'state vectors': from (1) State_p_vector (on top, id = 0) to (4) State_n_vector (bottom, id = 3)
std::vector<std::vector<int>> Configuration_State = { State_p_vector, State_f_vector, State_e_vector, State_n_vector },
    Configuration_cState = { State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector }; //  is the list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_State but for 'cracked' (or induced) network of k-cells
/* where 'n' :: "nodes", 'e' :: "edges", 'f' :: "faces", and 'p' :: "polyhedrons" */

/// 'CPD Tutorial' :: educational course active in the 'TUTORIAL' (please see 'config/main.ini' file) code execution mode
void tutorial();

/// Time interval variables for different parts (modulus) of the CPD code
double Main_time = 0.0, S_time = 0.0, P_time = 0.0, C_time = 0.0, M_time = 0.0, K_time = 0.0, W_time = 0.0;

/// * MODULES and LIBRARIES * ///
// '#include' all the main project's modules and libraries // IMPORTANT! Each module here is just a C++ library

/// * A task for the future: include all the code modules and functions as a single C++ library here placed in the STL directory *//
/// #include <processing_lib>

/* Various useful functions  */
#include "lib/PCC_Support_Functions.h" // It must be here - first in this list (!)

/* Various set measures */
#include "lib/PCC_Measures.h"

/* Objects library contains classes of various objects related to PCC's substructures and k-cells */
 #include "lib/PCC_Objects.h"

/* Section module calculates reduced PCC subcomplexes (including plain cuts) inheriting reduced sequences of special cells and 'state vectors' of the original PCC */
#include "lib/PCC_Section/PCC_Subcomplex.h"

/* Processing module assigned special IDs for the various elements (Nodes, Edges, Faces, Polytopes/Polyhedrons) of the space tessellation */
/* Output: module generates a design_sequences as the lists containing the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k={0,1,2,3}, corresponding to different generation principles (random, maximum entropy,.. etc.),
 * (2) IMPOSED: m-Cells (where m < k) directly labelled based on the HIGH-ORDER k-Cell IDs
 * (3) INDUCED: i-Cells generated as a result of some KINETIC process. They are always DEPEND on the ASSIGNED 'design_sequences' of special k-Cells (and maybe also imposed special m-Cells structures) */
#include "lib/PCC_Processing/PCC_Processing.h"

/* Writer module performs formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

/// Initial configuration as an object of the Config class described in Objects.cpp project library
Config initial_configuration;
// void initial_configuration(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> PCCpaths, std::vector<vector<int>> Configuration_State, std::vector<vector<int>> Configuration_cState); // Read the 'initial configuration' of the problem set in all the relevant '*.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)

///* ........................................................................................    Main    ................................................................ *///
//* (.h files) * @brief, @param and @return
//* (.cpp files) * @details (detailed descriptions)
/*!
* @brief Implement the whole program execution according to the specifications written in 'config/*.ini' files. I the LIST mode call the execution of the project modules one by one and compute the execution time for each of them.
* @param void
* @return 0, if successful
*/
int main() {
    cout << "------------------------------------------------------------------------------------------------" << endl;
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function
    Out_logfile_stream << "------------------------------------------------------------------------------------------------" << endl;

/// Initial configuration reader and information output to the screen and into the '.log' file
    initial_configuration.Read_config();
    std::vector<int> ConfigVector = initial_configuration.Get_ConfVector();
    /// alternatively :: set 'initial_configuration' object "manually"
//    initial_configuration.Set_config(const std::vector<int> &initial_configuration.Get_ConfVector(), const std::string &initial_configuration.Get_source_dir(), int &initial_configuration.Get_dim(), std::vector<char*> &initial_configuration.Get_paths(), std::vector<vector<int>> &initial_configuration.Get_Configuration_State(), std::vector<vector<int>> &initial_configuration.Get_Configuration_cState());

///    face_coordinates_vector = VectorDReader(PCCpaths.at(12?));
///    grain_coordinates_vector = VectorDReader(PCCpaths.at(9));

    // ===== Elapsing time Main ================
    unsigned int Mn_time = clock();
    Main_time = (double) Mn_time;
    cout << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
    cout << "Main execution time before modules is equal to  " << Main_time/ pow(10.0,6.0) <<  "  seconds" << endl;     Out_logfile_stream << "Main execution time before modules is equal to  " << Main_time/ pow(10.0,6.0) <<  "  seconds" << endl;

/// =========== TUTORIAL feature to facilitate the first acquaintance with the code ===================================
    if (e_mode == "tutorial") // only in the 'tutorial' execution type -  'etype' variable in the main.ini config file.
        tutorial();

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
// In the LIST mode all the functions are calling one after another without additional loops and intermediate data output
// For all the more complicated simulation cases the TASK mode (see it following next after the 'LIST' module)
CellsDesign new_cells_design; // A class (described in PCC_Objects.h) containing (1) all special k-cell sequences and (2) all the design_<*>_vectors for all k-cells

if ( main_type == "LIST"s ) { // 'LIST module'
// I: PCC_Section.h module
    if (ConfigVector.at(1) == 1) { // if PCC_Section is SWITCH ON in the main.ini file
            cout << " START of the PCC Subcomplex module " << endl;
            Out_logfile_stream << " START of the PCC Subcomplex module " << endl;

///* module */// Subcomplex();

// ===== Elapsing time Subcomplex ================
            unsigned int Subcomplex_time = clock();
            S_time = (double) Subcomplex_time - Main_time;
            cout << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl; Out_logfile_stream << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
    } // end if(SectionON)

/// II: PCC_Processing.h module
        if (ConfigVector.at(2) == 1) { // if PCC_Processing is SWITCH ON in the main.ini file
            cout << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
            cout << "START of the PCC Processing module " << endl; Out_logfile_stream << "START of the PCC Processing module " << endl;

            new_cells_design = PCC_Processing(Configuration_State, Configuration_cState);

// ===== Elapsing time Processing ================
            unsigned int Processing_time = clock();
            P_time = (double) Processing_time - S_time - Main_time;
            cout << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
        } // end if(ProcessingON)

/// VI: PCC_Writing module
    if (ConfigVector.at(6) == 1) { // simply output ON/OFF for the PCC_Writer module on the screen
            cout << "START of the PCC Writer module" << endl; Out_logfile_stream << "START of the PCC Writer module" << endl;

///            if(pcc_processed) PCC_Writer(new_cells_design, pcc_processed); else
            PCC_Writer(new_cells_design);

// ===== Elapsing time Writer ================
            unsigned int Writer_time = clock();
            double W_time = (double) Writer_time - Main_time- C_time - S_time - P_time - M_time - K_time;
            cout << "Writer time is equal to  " << W_time/ pow(10.0,6.0) <<  "  seconds" << endl;
    } // end if(WriterON)

} /// end of the SIMULATION MODE "LIST" in the main.ini file

/// ==========================================================================================================================================
/// ================================================= THE TASK MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
else if ( main_type == "TASK"s ) {

/// #include.. .cpp

    } /// end of the SIMULATION MODE "TASK" in the main.ini file
/// ==========================================================================================================================================

/// ================ Total Elapsing time ================ ///
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "-------------------------------------------------------------------------" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
    cout << dim << "D " << "runtime is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    Out_logfile_stream << dim << "D " << "runtime is equal to  " << fulltime / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;

/// Closing off-stream *.log file
    Out_logfile_stream.close();

    return 0;
} /// The END of Main function

/// ================================== FUNCTIONS DEFINED IN MAIN MODULE ==================================///

/// ================== # 7 # Initial configuration reader function (for main.cpp module) ==================
// void initial_configuration(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> PCCpaths, std::vector<vector<int>> Configuration_State, std::vector<vector<int>> Configuration_cState) {
//} /// END of the 'initial_configuration' function

/// ====================# 2 #============================= TUTORIAL ================================================= ///

void tutorial(){
    cout << "Hello there! This is the TEACHING mode of the program execution (!), \n"
    "where the user passes tutorial helping they to understand how to work with the CPD code. \n"
    "It can be changed by replacing the 'emode' variable with 'scientific' \n"
    "instead of the current 'tutorial' contained in the 'execution_mode' in ../config/main.ini file (!)" << endl;
    cout << endl;
    cout << "Continue? Please type 'Y' for 'yes' or 'N' for 'no' " << endl;
    char if_continue = 'Y'; cin >> if_continue;
    if (if_continue == 'N') exit(0);
    } /// END of the tutorial() function

                            /// *** H E A P *** ///