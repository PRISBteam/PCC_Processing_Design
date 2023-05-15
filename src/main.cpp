///******************************************************************************************************************************///
///*                                                                                                                           *///
///****************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   ***************************************///
///**************************************************************************************************************************///
///************************   Polyhedral Cell Complex (PCC) Processing Design :: (CPD code) (c)   **************************///
///************************************************************************************************************************///
// Material :: quadruple points, triple junctions, grain boundaries, and grains
// Tessellation :: nodes, edges, faces, polytopes
// PCC :: k-cells with their measures (areas, volumes) and barycenter coordinates

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
#include <thread>
#include <execution>
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
// #include "/usr/local/Cellar/libomp/15.0.7/include/omp.h"
#include <omp.h>

/// Simple reader (library) for .ini files and specific DPD code-related functions for reading its .ini files
#include "lib/ini/ini.h"
#include "lib/ini/ini_readers.h"

// * For the future -- include the whole CPD modules and functions as a single C++ library here placed in the STL directory *//
//include { Processing + SupportFunctions + Objects + Characterisation + Writer + Section + Multiphysics + Kinetic modules }
///#include <CPD_libs>

///------------------------------------------
using namespace std; //Standard namespace
using namespace Eigen; //Eigen namespace
using namespace Spectra; //Spectra namespace

/// Eigen library based Triplets class containing objects in the form T = T(i, j, value), where i and j are element's a(i, j) indices in the corresponding dense matrix and the third variable is its value
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type with name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles with name - SpMat
typedef MatrixXd DMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles with name - SpMat

/// * ------------------------------------------------------------------------------- -------------------------- *
/// * ============================================ GLOBAL VARIABLES ============================================ * ///
/// * ----------- Declaration of GLOBAL variables (can be used in all the project modules and libraries)-------- *
// Technical variables
/// File path to the configuration profile (MainPath is the path to the code's main directory)
string MainPath = "../";
string source_path = MainPath + "config/"s; char* sourcepath = const_cast<char*>(source_path.c_str()); // file path with main.ini

vector<char*> paths; // the vector with paths to all the PCC's matrices, sizes and other elements related to the input and output directories
string source_dir, output_dir; // input and output directories read from the file config_main.txt file
string sim_task; // path to the corresponding .cpp file (with "..") with a simulation task (form TASK main mode) as it is written in the main.ini file

// Global .log file output
ofstream Out_logfile_stream; //log.txt file output for the whole computation process

/// PCC related variables
// General
int dim; // problem dimension (2D or 3D)
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs with the number of k-cells of different types; will be read from file

// PCC Complex Geometry
// Global vectors of Descartes Coordinates for: (1) vertices, (2) barycenters of edges, (3) barycenters of faces and (4) barycenters of polyhedrons
vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;

// Measures
std::vector<double> edge_lengths_vector, face_areas_vector, polyhedron_volumes_vector;

/// PCC special structure related variables
// Sequences
std::vector<unsigned int> special_n_Sequence, special_e_Sequence, special_f_Sequence, special_p_Sequence;
std::vector<vector<unsigned int>> special_cells_design(4); // is a list of special_<cells>_sequences as an output of the Processing module
// First line here is the sequence of special nodes Sequence_n_vector;
// Second line -- sequence of special edges Sequence_e_vector,
// Third -- sequence of special faces Sequence_f_vector;
// Forth -- sequence of special polyhedrons Sequence_p_vector.

// Configurations and  states
// State_Vector in the form : [Element index] - > [Type] = kind of a code related to the microstructure of PCC
// Configuration_State is a list of all state-vectors: from (1) State_p_vector to (4) State_n_vector
std::vector<int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // based on each sequences in the list, the State_<*>_vector's of special cells can be calculated, including State_n_vector, State_e_vector, State_f_vector and State_p_vector
std::vector<vector<int>> Configuration_State; //  is the list of all mentioned below State_<*>_vectors as a possible output of the Processing module
/* where n = "nodes", e = "edges", f = "faces" and p = "polyhedrons" */

/// * MODULES and LIBRARIES * ///
// #include all the main project's modules and libraries
// IMPORTANT! Each module here is just a library (*_lib.h) or simple function(s) definition without any additional code surrounded it

/* Various useful functions (it must be here - first in this list! ) */
#include "lib/PCC_SupportFunctions.h"

/* Objects library contains classes of various objects related to PCC substructures and elements (it must be here - second in this list! ) */
#include "lib/PCC_Objects.h"

/* Section module calculates reduced PCC subcomplexes (including plain cuts) inheriting reduced sequences of special cells and State vectors of the original PCC */
//#include "lib/PCC_Section/PCC_Subcomplex.h"

/* Processing module assigned special IDs for the various elements (Vertices, Edges, Faces, Polyhedrons) of the space tessellation */
/* Output: module generates an s_sequences as the lists with the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k= {0,1,2,3}, corresponding to different generation principles (random, maximum entropy, etc.),
 * (2) IMPOSED m-Cells (where m < k) by HIGH-ORDER k-Cells, and
 * (3) INDUCED i-Cells by some KINETIC process DEPENDENT on the assigned special k-Cells and imposed special m-Cells structures */
#include "lib/PCC_Processing/PCC_Processing.h"

/* Characterisation module calculates various structural characteristics of special substructures defined on the PCC elements */
#include "lib/PCC_Characterisation/PCC_StructureCharacterisation.h"

/* Multiphysics module assign various physical quantities (energies, temperature, electrical conductivity) to the k-Cells of the PCC's subcomplexes */
//#include "lib/PCC_Multiphysics/PCC_Multiphysics.h"

/* Kinetic module assign some new values for the scalar or vector variables defined on the (Vertices, Edges, Faces, Polyhedrons) of the space tessellation (or new types of identifiers (IDs) of PCC's k-cells) */
/* As its output module currently generates one or several sequences of "fractured" cells' ID's; the new "fractured" State_cvector of faces (2-cells) can be generated */
//#include "lib/PCC_Kinetic/PCC_Kinetic.h"

/* Writer module perform formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

/// PCC special "elements" related to Mutiphysics module
// #1# 0-cells "precipitates"
// #2# 1-cells "disclinations"
// #3# 2-cells "macro-cracks"
// std::vector <macrocrack> large_cracks_vector;
// #4# 3-cells "phases"

///* ........................................................................................    Main    ................................................................ *///
int main() {
    cout << "-------------------------------------------------------------------------------------" << endl;

/// Useful ofstreams cleaning for data output
    // # 1 # LogFile ofstream
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::trunc); // will close at the end of the main (!)
    Out_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    string main_type; /// SIMULATION MODE(LIST) or 'SIMULATION MODE(TASK) :: This define the global simulation mode: "LIST" for the "list" of modules implementing one by one (if ON) and "TASK" for the user-defined task scripts with modules and functions from the project's libraries
// ConfigVector (../config/main.ini) contains all the control variables needed for the program
    vector<int> ConfigVector = config_reader_main(source_path, source_dir, output_dir, main_type, Out_logfile_stream);

        dim = ConfigVector.at(0); // space dimension of the problem (dim = 2 or 3);

/// Below the file names with the sparse PCC matrices which must already exit in the input_dir and have the same names (!)
// They can be obtained by the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or PCC Structure Generator tool (https://github.com/PRISBteam/PCC_Structure_Generator) for plasticity problems
//  PCC matrices
    string ssd0 = source_dir + "A0.txt"s, ssd1 = source_dir + "A1.txt"s, ssd2 = source_dir + "A2.txt"s, ssd3 = source_dir + "A3.txt"s, ssd4 = source_dir + "B1.txt"s, ssd5 = source_dir + "B2.txt"s, ssd6 = source_dir + "B3.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str())); paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));

//  PCC measures
//    std::vector<double> edge_lengths_vector;
    string ssd7 = source_dir + "face_areas.txt"s;
    paths.push_back(const_cast<char*>(ssd7.c_str()));
    string ssd8 = source_dir + "polyhedron_volumes.txt"s;
    paths.push_back(const_cast<char*>(ssd8.c_str()));

//  PCC geometry (barycenter coordinates)
//    vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;
    string ssd9 = source_dir + "polyhedron_seeds.txt"s; // grain (polyhedron) seeds
    paths.push_back(const_cast<char*>(ssd9.c_str()));
    string ssd10 = source_dir + "vertex_coordinates.txt"s; // vertex coordinates
    paths.push_back(const_cast<char*>(ssd10.c_str()));
    string ssd11 = source_dir + "face_normals.txt"s; // face normal vectors
    paths.push_back(const_cast<char*>(ssd11.c_str()));
///** add here edge_coordinates and face_coordinates

/// Vector with rhe numbers of PCC k-cells for k\in{0,1,2,3} from file
    if (is_file_exists(source_dir + "voro_Ncells.txt"s)) {
        string ncells = source_dir + "voro_Ncells.txt"s;
        char* number_of_cells = const_cast<char *>(ncells.c_str());
// :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Polyhedrons number
        CellNumbs = VectorIReader(number_of_cells); // VectorReader is a function from the PCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file

// Special feature for the 2D case ("grains" -> 2-cells; "faces" -> 1-cells; "triple junctions (lines) -> 0-cells" in 2D):
/// WARNING: (Very important!) // It needs to make everything similar in the code for 2D and 3D cases with the same output of the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) software (the same set of A0, A1..., B2 matrices)
        if (dim == 2) {
            CellNumbs.push_back(CellNumbs.at(2)); CellNumbs.at(2) = CellNumbs.at(1); CellNumbs.at(1) = CellNumbs.at(0); CellNumbs.at(0) = NULL;
        }

/// CellNumbs output
cout << "=====================================================================================" << endl;  Out_logfile_stream << "==========================================================================================================================================================================" << endl;
        unsigned int t_length = 0, l_length = 0;
        for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;
        for (int j : CellNumbs) Out_logfile_stream << l_length++ << "-cells #\t" << j << endl;

    } // if voro_Ncells exists
    else cout << "ERROR: The file " << source_dir + "voro_Ncells.txt"s << " does not exists (!)" << endl;

/// Initial state::
State_p_vector.resize(CellNumbs.at(3 + (dim-3)),0); State_f_vector.resize(CellNumbs.at(2 + (dim-3)),0); State_e_vector.resize(CellNumbs.at(1 + (dim-3)),0); State_n_vector.resize(CellNumbs.at(0),0);
Configuration_State.push_back(State_n_vector); Configuration_State.push_back(State_e_vector); Configuration_State.push_back(State_f_vector); Configuration_State.push_back(State_p_vector); // the order here matters!
// Configuration_State in 3D: [0] - nodes, [1] - edges, [2] - faces, [3] - polyhedrons

/// Output paths.vector to console and logfile out
    int npath = 0, lpath=0;
    cout << "_____________________________________________________________________________________" << endl;
    Out_logfile_stream << "_____________________________________________________________________________________" << endl;
    for (auto m : paths) cout <<"[" << npath++ << "]" << " paths:\t" << m << endl;
    for (auto p : paths) Out_logfile_stream <<"[" << lpath++ << "]" << " paths:\t" << p << endl;

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
// In the LIST mode all the functions are calling one after another without loops
double S_time = 0.0, P_time = 0.0, C_time = 0.0, M_time = 0.0, K_time = 0.0, W_time = 0.0; // time interval variables for different modulus
CellsDesign new_cells_design; // an object (described in PCC_Objects.h) containing (1) all k-cell sequences and (2) all design_x_vectors for all k-cells

    if ( main_type == "LIST"s ) {
// I: PCC_Section.h module
        if (ConfigVector.at(1) == 1) { // if PCC_Section is SWITCH ON in the main.ini file
            cout << " START of the PCC Subcomplex module " << endl;
            Out_logfile_stream << " START of the PCC Subcomplex module " << endl;

     // Subcomplex();

            // ===== Elapsing time Subcomplex ================
            unsigned int Subcomplex_time = clock();
            S_time = (double) Subcomplex_time;
            cout << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Section time is equal to  " << S_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(SectionON)

/// II: PCC_Processing.h module
        if (ConfigVector.at(2) == 1) { // if PCC_Processing is SWITCH ON in the main.ini file
            cout << " START of the PCC Processing module " << endl;
            Out_logfile_stream << " START of the PCC Processing module " << endl;
            new_cells_design = PCC_Processing(Configuration_State);
/*
            // update global vectors
            if(new_cells_design.Get_p_sequence().size() > 0) {
                special_p_Sequence = new_cells_design.Get_p_sequence();
                special_cells_design.at(3) = special_p_Sequence;
            }
            if(new_cells_design.Get_f_sequence().size() > 0) {
                special_f_Sequence = new_cells_design.Get_f_sequence();
                special_cells_design.at(2) = special_f_Sequence;
            }
            if(new_cells_design.Get_e_sequence().size() > 0) {
                special_e_Sequence = new_cells_design.Get_e_sequence();
                special_cells_design.at(1) = special_e_Sequence;
            }
            if(new_cells_design.Get_n_sequence().size() > 0) {
                special_n_Sequence = new_cells_design.Get_n_sequence();
                special_cells_design.at(0) = special_n_Sequence;
            }
*/

// ===== Elapsing time Processing ================
            unsigned int Processing_time = clock();
            P_time = (double) Processing_time - S_time;
            cout << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Processing time is equal to  " << P_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(ProcessingON)

/// III: DCC_Characterisation module
        ProcessedComplex pcc_processed;
        if (ConfigVector.at(3) == 1) {
            cout << "START of the PCC Structure Characterisation module" << endl;
            Out_logfile_stream << "START of the PCC Structure Characterisation module" << endl;

            pcc_processed = PCC_StructureCharacterisation(new_cells_design);

// ===== Elapsing time DCC_Characterisation ================
            unsigned int Characterisation_time = clock();
            C_time = (double) Characterisation_time - S_time - P_time;
            cout << "Characterisation time is equal to  " << C_time/pow(10.0,6.0) <<  "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
        }// end if(CharacterisationON)

// IV: PCC_Multiphysics.h module
        if (ConfigVector.at(4) == 1) { // if PCC_Multiphysics is SWITCH ON in the main.ini file
            cout << "START of the PCC Multiphysics module" << endl;
            Out_logfile_stream << "START of the PCC Multiphysics module" << endl;

            // Multiphysics();

            // ===== Elapsing time Multiphysics ================
            unsigned int Multiphysics_time = clock();
            M_time = (double) Multiphysics_time - C_time - S_time - P_time;
            cout << "Section time is equal to  " << M_time/ pow(10.0,6.0) <<  "  seconds" << endl; cout << "-------------------------------------------------------" << endl;
            Out_logfile_stream << "Section time is equal to  " << M_time/ pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------" << endl;
        } // end if(MultiphysicsON)

// V: DCC_Kinetic module
        if (ConfigVector.at(5) == 1) { // if PCC_Kinetic is SWITCH ON in the main.ini file
            cout << "START of the PCC Kinetic module" << endl;
            Out_logfile_stream << "START of the PCC Kinetic module" << endl;

            //kface_sequence = DCC_Kinetic(special_faces_sequence, K_type);

// ===== Elapsing time Kinetic ================
            unsigned int Kinetic_time = clock();
            K_time = (double) Kinetic_time - C_time - S_time - P_time - M_time;
            cout << "Kinetic time is equal to  " << K_time/ pow(10.0,6.0) <<  "  seconds" << endl;
            cout << "-------------------------------------------------------" << endl;
        }// end if(KineticON)

/// VI: PCC_Writing module
        if (ConfigVector.at(6) == 1) { // simply output ON/OFF for the PCC_Writer module on the screen
            cout << "START of the DCC Writer module" << endl;
            PCC_Writer(new_cells_design, pcc_processed);

// ===== Elapsing time Writer ================
            unsigned int Writer_time = clock();
            double W_time = (double) Writer_time - C_time - S_time - P_time - M_time - K_time;
            cout << "Writer time is equal to  " << W_time/ pow(10.0,6.0) <<  "  seconds" << endl;
        } // end if(WriterON)

    } /// end SIMULATION MODE else (LIST)

/// ==========================================================================================================================================
/// ================================================= THE TASK MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
// #include sim_task (??? instead preprocessor related #include - call during the simulation time )

// # 1 # Simple sfaces processing
//#include "tasks/task_sFacesProcessing.cpp"

// # 2 # Processing design for Strips of Inclusions (with the parameters of average \mu and dispersion \sigma of the strip lengths distribution)
//#include "tasks/task_StripsProcessingDesign.cpp"

// # 3 # Macrocrack growth with multiple cracking simulations
//#include "tasks/task_macrocrack.cpp"

/// ==========================================================================================================================================

/// ================ Total Elapsing time ================ ///
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "-------------------------------------------------------------------------" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
    cout << dim << "D " << "runtime is equal to  " << fulltime/pow(10.0,6.0) <<  "  seconds" << endl;
    Out_logfile_stream << dim << "D " << "runtime is equal to  " << fulltime/pow(10.0,6.0) <<  "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;
    Out_logfile_stream << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the PCC Processing\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;

/// Closing log_file ofstream
    Out_logfile_stream.close();

    return 0;
} /// The END of Main function

/// ================================== FUNCTIONS DEFINED IN MAIN MODULE ==================================///

/// NEW HEAP
/*
/// Principal variables
///_____________________________________________________________________________________
std::vector<unsigned int> kface_sequence; // special Faces sequence for the Kinetic module
/// All grain coordinates
string GCpath_string = source_dir + "grain_seeds.txt"s;
char* GCpath = const_cast<char*>(GCpath_string.c_str());
vector<tuple<double, double, double>> grain_coordinates(CellNumbs.at(3));
grain_coordinates = TuplesReader(GCpath);
/// Set grain_vector_coordinates
for(auto grain_coord_tuple : grain_coordinates)
grain_coordinates_vector.push_back(grain_coord_tuple);

//REPAIR    cout << "Size of grain_seeds: " << gs_size << " CellNumbs.at(3): " << CellNumbs.at(3) << endl;

/// All face coordinates /// LONG CALCULATION HERE
string FSpath_string = source_dir + "face_seeds.txt"s;
char* FSpath = const_cast<char*>(FSpath_string.c_str());
face_coordinates_vector.clear();
face_coordinates_vector.resize(CellNumbs.at(2), make_tuple(0,0,0));
//REPAIR    cout << "face vector size: " << face_coordinates_vector.size() << endl;

    if (is_file_exists(FSpath_string)) {
        face_coordinates_vector = TuplesReader(FSpath);
    }
    else {
        cout << "Calculation of face coordinates started now by find_aGBseed function" << endl;
        #pragma omp parallel for // parallel execution by OpenMP
        for (unsigned int fnumber = 0; fnumber < CellNumbs.at(2)-1; ++fnumber) {
            face_coordinates_vector.at(fnumber) = find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates);
//REPAIR
     if(fnumber % 10 == 0)  cout << "Total # of faces " << CellNumbs.at(2)<< " face # " << fnumber << " "<< get<0>(face_coordinates_vector.at(fnumber)) << " "
                 << get<1>(face_coordinates_vector.at(fnumber)) << " " << get<2>(face_coordinates_vector.at(fnumber)) << endl;
        }
        cout << "Calculations was done!" << endl;
        ofstream OutFaceCoord; OutFaceCoord.open(FSpath, ios::trunc);
        if(OutFaceCoord.is_open())
            #pragma omp parallel for // parallel execution by OpenMP
            for (auto fcoord : face_coordinates_vector)
                OutFaceCoord << get<0>(fcoord) << " " << get<1>(fcoord) << " " << get<2>(fcoord) << endl;
    }


/// All vertex coordinates
string VCpath_string = source_dir + "vertex_seeds.txt"s;
char* VCpath = const_cast<char*>(VCpath_string.c_str());
vertex_coordinates_vector.resize(CellNumbs.at(0), make_tuple(0,0,0));
vertex_coordinates_vector = TuplesReader(VCpath);
//REPAIR    cout << "vertex_coordinates_vector size: " << vertex_coordinates_vector.size() << endl; //exit(0);

/// All grain volumes
string GVpath_string = source_dir + "grain_volumes.txt"s;
char* GVpath = const_cast<char*>(GVpath_string.c_str());
polyhedron_volumes_vector.resize(CellNumbs.at(3));
polyhedron_volumes_vector = dVectorReader(GVpath);

/// All face areas
string GBApath_string = source_dir + "face_areas.txt"s;
char* GBApath = const_cast<char*>(GBApath_string.c_str());
face_areas_vector.resize(CellNumbs.at(3));
face_areas_vector = dVectorReader(GBApath);

//---------------------------------------------------------------------------
*/

/// OLD HEAP ///
/*
string line;
ifstream inConf(sourcepath);
/// I /// General paths and initial parameters
if (inConf) { // If the file was successfully open, then
while (getline(inConf, line, '\n')) {
// Number of types
for (auto it: line) {
// integers
if (it == '^') res.push_back(line.at(0) - '0'); // dimension of the problem; res[0]
// strings
if (it == '!') {
stringstream line1_stream(line);
line1_stream >> main_type;
} // LIST or TASK simulation type
else if (it == '&' && main_type == "TASK") {
stringstream line2_stream(line);
line2_stream >> sim_task;
} // sumulation task path (if TASK main type only!)
else if (it == '@') {
stringstream line3_stream(line);
line3_stream >> source_dir;
} // input folder path input_dir = const_cast<char*>(input.c_str());
else if (it == '$') {
stringstream line4_stream(line);
line4_stream >> output_dir;
} // output folder path input_dir = const_cast<char*>(input.c_str());

} // end for (auto it: line)
} // end while (getline..)
} else cout << "The file " << sourcepath << " cannot be read" << endl; // If something goes wrong
// end if (inConf)

/// II /// Different modulus On/Off states///
bool isSectionON = 0, isProcessingON = 0, isCharacterisationON = 0, isKineticON = 0, isMultiphysicsON = 0, isWriterON = 0;

if (inConf) { //If the file stream was successfully open, then
while (getline(inConf, line, '\n')) {   // REPAIR            cout << line << endl;
// (1)
if (!line.compare("PCC_Section SWITCHED ON"s))
isSectionON = 1;
// (2)
else if (!line.compare("PCC_Processing SWITCHED ON"s))
isProcessingON = 1;
// (3)
else if (!line.compare("PCC_Characterisation SWITCHED ON"s))
isCharacterisationON = 1;
// (4)
else if (!line.compare("PCC_Multiphysics SWITCHED ON"s))
isMultiphysicsON = 1;
// (5)
else if (!line.compare("PCC_Kinetic SWITCHED OFF"s))
isKineticON = 1;
// (6)
else if (!line.compare("PCC_Writer SWITCHED OFF"s))
isWriterON = 1;

} // end while (getline..)
} // end if (inConf)
else
cout << "ProcessingON() error: The file " << sourcepath << " cannot be read" << endl; // If something goes wrong
*/

/*
// Reading and Output of the configuration file
std::vector<double> config_reader_main(char* config, string &Subcomplex_type, string &Processing_type, string &Kinetic_type, string &source_dir, string &output_dir) {
    vector<double> res;

    string line;
    double p_max = 0, f_max = 0;

    ifstream inConf(config);
    if (inConf) { // If the file was successfully open, then
        while (getline(inConf, line, '\n')) {
            // Number of types
            for (auto it: line) {
                if (it == '}') Subcomplex_type = line.at(0); // simulation Section type
                else if (it == '&') Processing_type = line.at(0); // simulation Processing type
                else if (it == '`') Kinetic_type = line.at(0); // simulation Kinetic type

                else if (it == '@') res.push_back(line.at(0) - '0'); // dimension of the problem; res[1]
                else if (it == '!') res.push_back(line.at(0) - '0'); // number of i(X) designs; res[0]
                    ///?? does it properly working with 10, 100 etc
                    //else if (it == '!') res.push_back(line.at(0) - '0'); // number of special Face types; res[0]

                else if (it == '?') {
                    stringstream line1_stream(line);
                    line1_stream >> p_max;
                    res.push_back(p_max);
                } // MAX fraction of Faces (calculation limit); res[2]
                else if (it == '^') {
                    stringstream line2_stream(line);
                    line2_stream >> f_max;
                    res.push_back(f_max);
                } // MAX fraction of Cracks (calculation limit); res[3]
                else if (it == '~') {
                    stringstream line3_stream(line);
                    line3_stream >> source_dir;
                } // input folder path input_dir = const_cast<char*>(input.c_str());
                else if (it == '$') {
                    stringstream line4_stream(line);
                    line4_stream >> output_dir;
                } // output folder path input_dir = const_cast<char*>(input.c_str());

                else if (it == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
                else if (it == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
            }
        }

    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

    res.size();
    cout << "The problem dimension that is the maximum dimension k_max of k-Cells\t | "s << res.at(0) << endl;
    cout << "The number of i(X) designs (0 - random, 1 - random + smax,.. )\t\t\t | " << res.at(1) << endl;
    if (SubcomplexON(config, 0)) cout << "Subcomplex type ('P' or 'H')\t\t\t\t\t\t\t\t\t | "s << Subcomplex_type << endl;
    if (ProcessingON(config, 0)) cout << "Processing type ('R', 'S', 'C' or 'X')\t\t\t\t\t\t\t\t\t | "s << Processing_type << endl;
    //cout << "Calculation type ('R', 'RW', 'S', 'F', 'I' or 'X'):\t\t\t\t\t\t | "s << Processing_type << endl;
    //cout << "The number of special Face (2-cells) types:\t\t\t\t\t\t\t\t | " << res.at(1) << endl;
    cout << "MAX fraction of Faces (calculation limit) \t\t\t\t\t\t\t\t | " << res.at(2) << endl;
    if (KineticON(config, 0)) cout << "Kinetic type ('W' or 'I')\t\t\t\t\t\t\t\t\t\t\t\t | "s << Kinetic_type << endl;
    if (KineticON(config, 0)) cout << "MAX fraction of Cracks (calculation limit)\t\t\t\t\t\t\t\t | " << res.at(3) << endl;
    cout << endl;
    cout << "Input folder:\t" << source_dir << endl;
    cout << "Output folder:\t" << output_dir << endl;
    cout << endl;
    if (CharacterisationON(config, 0)) {
        cout << "Nodes types statistics, indices and configuration entropy:     "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Edges types statistics, indices and configuration entropy:     "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Faces types statistics and structural indices:                 "; if (res.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Grain types statistics, indices and configuration entropy:     "; if (res.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Node Laplacian with its spectrum for the Nodes network:       "; if (res.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Edge Laplacian with its spectrum for the Edges network:       "; if (res.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Face  Laplacian with its spectrum for the Faces network:       "; if (res.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Grain Laplacian with its spectrum for the Grains network:      "; if (res.at(11) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
        cout << "Tutte polynomial for the special network:                      "; if (res.at(12) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    }

    return res;
}
*/

// void eraseSubStr(std::string & mainStr, const std::string & toErase); // support (technical) function: erase the first occurrence of given substring from main string
// Still not working in Windows:(
// using std::__fs::filesystem::current_path; // to obtain current working directory
// string MainPath = current_path(); // MainPath variables
// eraseSubStr(MainPath, "cmake-build-debug"s); // to delete .../cmake-build-debug/ from the path
