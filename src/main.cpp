///******************************************************************************************************************************///
///************************   Polyhedral Cell Complex (PCC) Processing Design :: (CPD code) (c)   ******************************///
///****************************************************************************************************************************///
///*                                        version 4.0 | 16/10/2023                                                          *///
///**************************************************************************************************************************///
///************************************ Dr Elijah Borodin, Manchester, UK **************************************************///
///**************************************** Spring 2022 - Winter 2023  ****************************************************///
///***********************************************************************************************************************///
///*    Code source: https://github.com/PRISBteam/PCC_Processing_Design/                                                     *///
///*    PCC sources: https://materia.team/                                                                                 *///

///* The project provides a reliable tool for 1. Obtaining, 2. Analysing and 3. Improving of the 'design vectors' as the sequences           *///
///*  of k-cells containing in the k-skeletons, k={0,1,2,3} of a Polyhedral Cell Complex (PCC). Such PCCs can be created by external        *///
///* codes based on the tessellation of 2D or 3D spaces by an agglomeration of polytopes (in the 2D case) or polyhedrons (in the 3D case). *///

// Key terminology:
// Material :: 'quadruple points', 'triple junctions', 'grain boundaries', and 'grains' with their orientations and barycenter coordinates
// Tessellation :: 'nodes, 'edges', 'faces', 'polytopes' with their measures (lengths, areas, volumes) and barycenter coordinates
// PCC :: 'k-cells' containing in 'k-skeletons' (k = {0,1,2,3}) with their degree fractions

///* ----------------------------------------- *
///* Standard C++ (STL) libraries
///* ----------------------------------------- *
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include <numeric>
#include <algorithm>
#include <cmath>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2022)
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

/// Spectra source: https://spectralib.org/ (2022)
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
// Included only in the parallelized version of the code

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project
#include "lib/ini/ini_readers.h"

///------------------------------------------
using namespace std; // STL namespace
using namespace Eigen; // Eigen library namespace
using namespace Spectra; // Spectra library namespace

/// Eigen library based Triplets class containing objects in the form T = T(i, j, value), where i and j are element's a(i, j) indices in the corresponding dense matrix and the third variable is its value
typedef Triplet<double> Tr; // <Eigen> library class, which declares a triplet type with the nickname 'Tr'
typedef SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'
typedef MatrixXd DMat; // <Eigen> library class, which declares a dense matrix type of doubles with the nickname 'DMat'

/// * ---------------------------------------------------------------------------------------------------------- *///
/// * ============================================ GLOBAL VARIABLES ============================================ *///
/// * ----------- Declaration of GLOBAL variables (can be used in all the project modules and libraries)-------- *///

/// Technical variables
// 'source_path' is a path to the directory containing *.ini files
string source_path = "../config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // path *.ini files as it is written in the main.ini file (!)
string e_mode = "tutorial"; // 'execution_type' (etype variable) from the config/main.ini file

vector<char*> paths; // The vector with the paths to all the PCC's matrices, measures and other supplementary data
string source_dir, output_dir; // Input and output directories as it is written in the file 'config/main.ini' file
string sim_task; // path to the corresponding *.cpp file with a simulation task (for TASK main mode only) as it is written in the main.ini file

/// Global 'log.txt' file output
ofstream Out_logfile_stream; // log.txt file output of the entire computation process as a copy of the console output

/// PCC - related variables
// General
int dim;// PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.

// Combinatorial
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs with the number of k-cells of different types 'k'. It will be readed from the 'number_of_cells.txt' file.
// First line here is the number of nodes (0-cells), second - edges (1-cells), third - faces (2-cells) (in the 2D and 3D cases only), fourth - polyhedra (3-cells) (in the 3D case only)
// In the 2D case "grains" -> 2-cells; "faces" -> 1-cells; "triple junctions (lines) -> 0-cells".

// Geometry
// Global vectors of Cartesian coordinates for: (1) vertex coordinates, (2) barycentres of edges, (3) barycentres of faces and (4) barycentres of polyhedrons
vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector; // vectors containing barycenter Cartesian coordinates of the corresponding tessellation's elements

// Measures
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector;

/// PCC special structure-related variables
// Configurations and  states

// State_Vector in the form : [Element index] - > [Element type]
// Configuration_State = {State_p_vector, State_f_vector, State_e_vector, State_n_vector} is a list of all 'state vectors': from (1) State_p_vector (on top, id = 0) to (4) State_n_vector (bottom, id = 3)
std::vector<int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // based on each sequences in the list, the State_<*>_vector's of special cells can be calculated, including State_n_vector, State_e_vector, State_f_vector and State_p_vector
std::vector<int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // based on each sequences in the list, the State_<*>_vector's of special cells can be calculated, including State_n_vector, State_e_vector, State_f_vector and State_p_vector
std::vector<vector<int>> Configuration_State = { State_p_vector, State_f_vector, State_e_vector, State_n_vector },
                         Configuration_cState = { State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector }; //  is the list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_State but for 'cracked' (or induced) network of k-cells
/* where n = "nodes", e = "edges", f = "faces" and p = "polyhedrons" */

/// Initial configuration reader
void initial_configuration(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<vector<int>> Configuration_State, std::vector<vector<int>> Configuration_cState);

/// 'CPD Tutorial' :: educational course active in the 'tutorial' (please see 'main.ini file') code execution mode
void tutorial();

/// Time interval variables for different parts (modulus) of the CPD code
double Main_time = 0.0, S_time = 0.0, P_time = 0.0, C_time = 0.0, M_time = 0.0, K_time = 0.0, W_time = 0.0;

/// * MODULES and LIBRARIES * ///
// '#include' all the main project's modules and libraries // IMPORTANT! Each module here is just a C++ library

// * A task for the future: include all the code modules and functions as a single C++ library here placed in the STL directory *//
/// #include <processing_lib>

/* Various useful functions  */
#include "lib/PCC_Support_Functions.h" // It must be here - first in this list (!)

/* Various set measures */
/// #include "lib/PCC_Measures.h"

/* Objects library contains classes of various objects related to PCC's substructures and k-cells */
 #include "lib/PCC_Objects.h"

/* Section module calculates reduced PCC subcomplexes (including plain cuts) inheriting reduced sequences of special cells and 'state vectors' of the original PCC */
/// #include "lib/PCC_Section/PCC_Subcomplex.h"

/* Processing module assigned special IDs for the various elements (Nodes, Edges, Faces, Polytopes/Polyhedrons) of the space tessellation */
/* Output: module generates a design_sequences as the lists containing the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k={0,1,2,3}, corresponding to different generation principles (random, maximum entropy,.. etc.),
 * (2) IMPOSED: m-Cells (where m < k) directly labelled based on the HIGH-ORDER k-Cell IDs
 * (3) INDUCED: i-Cells generated as a result of some KINETIC process. They are always DEPEND on the ASSIGNED 'design_sequences' of special k-Cells (and maybe also imposed special m-Cells structures) */
#include "lib/PCC_Processing/PCC_Processing.h"

/* Writer module performs formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

///* ........................................................................................    Main    ................................................................ *///
int main() {
    cout << "------------------------------------------------------------------------------------------------" << endl;
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function
    Out_logfile_stream << "------------------------------------------------------------------------------------------------" << endl;

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen and *.log file.
    std::string main_type; /// SIMULATION MODE(LIST) or 'SIMULATION MODE(TASK) :: This define the global simulation mode: 'LIST' for the "list" of modules implementing one by one (if ON) and 'TASK' for the user-defined task scripts with modules and functions included from the project's libraries

    std::vector<int> ConfigVector = config_reader_main(source_path, source_dir, output_dir, main_type, e_mode);
// ConfigVector (../config/main.ini) contains ALL the control variables needed for the program execution

/// Initial configuration reader and information output to the screen and into the .log file
    initial_configuration(ConfigVector, source_dir, dim, paths, Configuration_State, Configuration_cState);

///    face_coordinates_vector = VectorDReader(paths.at(12?));
///    grain_coordinates_vector = VectorDReader(paths.at(9));

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

/// Subcomplex();

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
//const std::ofstream &Out_logfile_stream
void initial_configuration(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<vector<int>> Configuration_State, std::vector<vector<int>> Configuration_cState) {

    dim = ConfigVector.at(0); // [0] space dimension of the problem (dim = 1, 2 or 3);
    if (dim != 1 && dim != 2 && dim != 3) {
        cout << "Wrong dimension ERROR! Please change 'dim' parameter in ../config/mail.ini file to 1, 2, or 3"
             << endl;
        Out_logfile_stream
                << "Wrong dimension ERROR! Please change 'dim' parameter in ../config/mail.ini file to 1, 2, or 3"
                << endl;
        exit(1);
    }

/// Several file paths to the sparse PCC's matrices which must already exit in the 'source_dir' and have the same names as below (!)
// They can be obtained by the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or PCC Structure Generator tool (https://github.com/PRISBteam/PCC_Structure_Generator) for plasticity problems
//  PCC matrices, metrices and coordinates
    std::string ssd0, ssd1, ssd2, ssd3, ssd4, ssd5, ssd6, ssd7, ssd8, ssd9, ssd10, ssd11, ssd12, ssd13, ssd14; // PCC matrices

    if (is_file_exists(source_dir + "algebraic/A0.txt"s) && is_file_exists(source_dir + "algebraic/A1.txt"s) &&
        is_file_exists(source_dir + "algebraic/B1.txt"s) && dim >= 0)
        ssd0 = source_dir + "algebraic/A0.txt"s, ssd1 = source_dir + "algebraic/A1.txt"s, ssd4 =
                source_dir + "algebraic/B1.txt"s; // 1D
    if (is_file_exists(source_dir + "/algebraic/A2.txt"s) && is_file_exists(source_dir + "algebraic/B2.txt"s) &&
        dim >= 2)
        ssd2 = source_dir + "algebraic/A2.txt"s, ssd5 = source_dir + "algebraic/B2.txt"s; // 2D
    if (is_file_exists(source_dir + "/algebraic/A3.txt"s) && is_file_exists(source_dir + "algebraic/B3.txt"s) &&
        dim == 3)
        ssd3 = source_dir + "algebraic/A3.txt"s, ssd6 = source_dir + "algebraic/B3.txt"s; // 3D
    if (dim > 3) {
        cout
                << "INPUT DATA ERROR (!) dim > 3 as it specified in the ../config/main.ini file. Please, make it equal to 1, 2 or 3."
                << endl;
        Out_logfile_stream
                << "INPUT DATA ERROR (!) dim > 3 in the ../config/main.ini file. Please, make it equal to 1, 2 or 3."
                << endl;
        exit(1);
    }
    //The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char *>(ssd0.c_str()));
    paths.push_back(const_cast<char *>(ssd1.c_str()));
    paths.push_back(const_cast<char *>(ssd2.c_str()));
    paths.push_back(const_cast<char *>(ssd3.c_str()));
    paths.push_back(const_cast<char *>(ssd4.c_str()));
    paths.push_back(const_cast<char *>(ssd5.c_str()));
    paths.push_back(const_cast<char *>(ssd6.c_str()));

//  PCC measures
/// ---> edge lengths HERE (!) // see next 'ssd14 = source_dir + "edge_lengths.txt"s;' and replace with 'ssd6'
    if (is_file_exists(source_dir + "measures/face_areas.txt"s) && dim >= 2)
        ssd7 = source_dir + "measures/face_areas.txt"s; // 2D
    paths.push_back(const_cast<char *>(ssd7.c_str()));
    if (is_file_exists(source_dir + "measures/polyhedron_volumes.txt"s) && dim == 3)
        ssd8 = source_dir + "measures/polyhedron_volumes.txt"s; // 3D
    paths.push_back(const_cast<char *>(ssd8.c_str()));

//  PCC geometry (barycenter coordinates)
//    For vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;
    if (is_file_exists(source_dir + "coordinates/polyhedron_seeds.txt"s) && dim == 3)
        ssd9 = source_dir + "coordinates/polyhedron_seeds.txt"s; // grain (polyhedron) seeds
    paths.push_back(const_cast<char *>(ssd9.c_str()));
    if (is_file_exists(source_dir + "coordinates/vertex_seeds.txt"s) && dim >= 1)
        ssd10 = source_dir + "coordinates/vertex_seeds.txt"s; // vertex coordinates
    paths.push_back(const_cast<char *>(ssd10.c_str()));
    if (is_file_exists(source_dir + "other/face_normals.txt"s) && dim >= 2)
        ssd11 = source_dir + "other/face_normals.txt"s; // face normal vectors
    paths.push_back(const_cast<char *>(ssd11.c_str()));
    if (is_file_exists(source_dir + "coordinates/edge_seeds.txt"s) && dim >= 1)
        ssd12 = source_dir + "coordinates/edge_seeds.txt"s; // edge barycentres coordinates
    paths.push_back(const_cast<char *>(ssd12.c_str()));
    if (is_file_exists(source_dir + "coordinates/face_seeds.txt"s) && dim >= 2)
        ssd13 = source_dir + "coordinates/face_seeds.txt"s; // face barycentres coordinates
    paths.push_back(const_cast<char *>(ssd13.c_str()));

    if (is_file_exists(source_dir + "measures/edge_lengths.txt"s) && dim >= 1)
        ssd14 = source_dir + "measures/edge_lengths.txt"s; // edge barycentres coordinates
    paths.push_back(const_cast<char *>(ssd14.c_str()));

/// Vector with rhe numbers of PCC k-cells for k\in{0,1,2,3} from file
    if (is_file_exists(source_dir + "/combinatorial/number_of_cells.txt"s)) {
        std::string ncells = source_dir + "/combinatorial/number_of_cells.txt"s;
        char *number_of_cells = const_cast<char *>(ncells.c_str());
// :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Polyhedrons number
        CellNumbs = VectorIReader(
                number_of_cells); // VectorReader is a function from the PCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file

/// CellNumbs output
        cout << "=====================================================================================" << endl;
        Out_logfile_stream
                << "=========================================================================================================================================================================="
                << endl;
        unsigned int t_length = 0;
        for (int cell_numb: CellNumbs) {
            cout << t_length++ << "-cells #\t" << cell_numb << endl;
            Out_logfile_stream << t_length << "-cells #\t" << cell_numb << endl;
        } // end for (int cell_numb : CellNumbs)
    } // if voro_Ncells exists
    else
        cout << "ERROR: The file " << source_dir + "/combinatorial/number_of_cells.txt"s
             << " does not exists (!)" << endl;

/// Initial state::
    if (dim < 3) { // Does not exist in 2D: CellNumbs.at(3)
        State_p_vector.resize(0, 0);
        State_pfracture_vector.resize(0, 0);
        if (dim == 1)
            State_f_vector.resize(0, 0);
        State_ffracture_vector.resize(0, 0);
    } else {
        State_p_vector.resize(CellNumbs.at(3), 0);
        State_pfracture_vector.resize(CellNumbs.at(3), 0);
        State_f_vector.resize(CellNumbs.at(2), 0);
        State_ffracture_vector.resize(CellNumbs.at(2), 0);
    }
    State_e_vector.resize(CellNumbs.at(1), 0);
    State_n_vector.resize(CellNumbs.at(0), 0);
    State_efracture_vector.resize(CellNumbs.at(1), 0);
    State_nfracture_vector.resize(CellNumbs.at(0), 0);

    if (dim == 1) { // 1D case
        Configuration_State.resize(2);
        Configuration_State = {State_n_vector, State_e_vector};
        Configuration_cState.resize(2);
        Configuration_cState = {State_nfracture_vector, State_efracture_vector};
    }
    if (dim == 2) { // 2D case
        Configuration_State.resize(3);
        Configuration_State = {State_n_vector, State_e_vector, State_f_vector};
        Configuration_cState.resize(3);
        Configuration_cState = {State_nfracture_vector, State_efracture_vector, State_ffracture_vector};
    } else { // 3D case
        Configuration_State.resize(4);
        Configuration_State = {State_n_vector, State_e_vector, State_f_vector, State_p_vector};
        Configuration_cState.resize(4);
        Configuration_cState = {State_nfracture_vector, State_efracture_vector, State_ffracture_vector,
                                State_pfracture_vector};
    } // Configuration_State in 3D: [0] -> nodes, [1] -> edges, [2] -> faces, [3] -> polyhedrons

/// Output paths.vector to console and logfile out
    int npath = 0;
    cout << "_____________________________________________________________________________________" << endl;
    Out_logfile_stream
            << "_____________________________________________________________________________________" << endl;
    for (auto path: paths) {
        cout << "[" << npath++ << "]" << " paths:\t" << path << endl;
        Out_logfile_stream << "[" << npath << "]" << " paths:\t" << path << endl;
    }
} /// END of the 'initial_configuration' function

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
            // HERE ARE OLD AND/OR UNUSED BUT USEFUL PIECES OF CODE