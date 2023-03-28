///******************************************************************************************************************************///
///*                                                                                                                           *///
///****************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   ***************************************///
///**************************************************************************************************************************///
///************************   Polyhedral Cell Complex (PCC) Processing Design :: (CPD code) (c)   **************************///
///************************************************************************************************************************///

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

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2022)
// #include <Eigen/Core>
// #include <Eigen/Dense>
#include <Eigen/SparseCore>

/// Spectra source: https://spectralib.org/ (2022)
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>

/// Open MP library https://www.openmp.org/resources/openmp-compilers-tools/
#include "/usr/local/Cellar/libomp/15.0.7/include/omp.h"

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

/// * ------------------------------------------------------------------------------- -------------------------- *
/// * ============================================ GLOBAL VARIABLES ============================================ * ///
/// * ----------- Declaration of GLOBAL variables (can be used in all the project modules and libraries)-------- *
int dim; // problem dimension (2D or 3D)
std::vector<unsigned int> CellNumbs; // the vector named CellNumbs with the number of k-cells of different types; will be read from file
string source_dir, output_dir; // input and output directories read from the file config_main.txt file
vector<char*> paths; // the vector with paths to all the PCC's matrices, sizes and other elements related to the input and output directories

std::vector<vector<unsigned int>> special_cells_design; // is a list of special_<cells>_sequences as an output of the Processing module. T
// First line here are numbers of special vertices;
// Second line -- numbers of special edges,
// Third -- numbers of special faces;
// Forth -- numbers of special polyhedrons.

/// PCC Complex Geometry
// Global vectors of Descartes Coordinates for: (1) vertices, (2) barycenters of edges, (3) barycenters of faces and (4) barycenters of polyhedrons
vector<tuple<double, double, double>> vertex_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, grain_coordinates_vector;

/// Measures
std::vector<double> grain_volumes_vector, face_areas_vector, edge_lengths_vector;

/// PCC special structures
std::vector<unsigned int> face_strip_distribution; // a vector containing length distribution of special faces (strips)

/// Global .log file output
ofstream Out_logfile_stream; //log.txt file output for the whole computation process

/// * MODULES and LIBRARIES * ///
// #include all the main project's modules and libraries
// IMPORTANT! Each module here is just a library (*_lib.h) or simple function(s) definition without any additional code around it

/* Various useful functions (it must be here - first in this list! ) */
#include "lib/PCC_SupportFunctions.h"

/* Objects library contains classes of various objects related to PCC substructures and elements (it must be here - second in this list! ) */
#include "lib/PCC_Objects.h"

/* Processing module assigned special IDs for the various elements (Vertices, Edges, Faces, Polyhedrons) of the space tessellation */
/* Output: module generates an s_sequences as the lists with the sequences of k-cells possessing "special" IDs including
 * (1) ASSIGNED: k-Cells, k= {0,1,2,3}, corresponding to different generation principles (random, maximum entropy, etc.),
 * (2) IMPOSED m-Cells (where m < k) by HIGH-ORDER k-Cells, and
 * (3) INDUCED i-Cells by some KINETIC process DEPENDENT on the assigned special k-Cells and imposed special m-Cells structures */
/* Based on each sequences in the list, the State_<*>vector's of special cells can be calculated, including State_nvector, State_evector, State_fvector and State_pvector */
/* The whole State_vectors_list can be described by the list of all the State vectors */
#include "lib/PCC_Processing/PCC_Processing.h"

/* Kinetic module assign some new values for the scalar or vector variables defined on the (Vertices, Edges, Faces, Polyhedrons) of the space tessellation (or new types of identifiers (IDs) of PCC's k-cells) */
/* As its output module currently generates one or several sequences of "fractured" cells' ID's; the new "fractured" State_cvector of faces (2-cells) can be generated */
#include "lib/PCC_Kinetic/PCC_Kinetic.h"

/* Section module calculates reduced PCC subcomplexes (including plain cuts) inheriting reduced sequences of special cells and State vectors of the original PCC */
#include "lib/PCC_Section/PCC_Subcomplex.h"

/* Multiphysics module assign various physical quantities (energies, temperature, electrical conductivity) to the k-Cells of the PCC's subcomplexes */
#include "lib/PCC_Multiphysics/PCC_Multiphysics.h"

/* Characterisation module calculates various structural characteristics of special substructures defined on the PCC elements */
#include "lib/PCC_Characterisation/StructureCharacterisation.h"

/* Writer module perform formatted output of various data structures generated by other modules */
#include "lib/PCC_Writer/PCC_Writer.h"

/// * PCC Special Elements * ///
// #1# 0-cells Precipitates
// #2# 1-cells Disclinations
// #3# 2-cells Macro-cracks
// std::vector <macrocrack> large_cracks_vector;
// #4# 3-cells Phases

///* ........................................................................................    Main    ................................................................ *///
int main() {
    cout << "-------------------------------------------------------------------------------------" << endl;

/// File path to the configuration profile (MainPath is the path to the code's main directory)
    string MainPath = "../";
    string config = MainPath + "config/config_main.txt"s; char* confpath = const_cast<char*>(config.c_str()); // file path with config_main.txt

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    // The source directory and simulation type from file config.txt
    string M_type; /// SIMULATION MODE(LIST) or 'SIMULATION MODE(TASK) :: This define the global simulation mode: "LIST" for the "list" of modules implementing one by one (if ON) and "TASK" for the user-defined task scripts with modules and functions from the project's libraries

/// ConfigVector contains all the control variables of the program readed from the config.txt
    std::vector<double> config_reader_main(char* config, string &Subcomplex_type, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder); // Read and output the initial configuration from the config.txt file
    vector<double> ConfigVector = config_reader_main(confpath, S_type, P_type, K_type, source_dir, output_dir);
    char* indir = const_cast<char*>(source_dir.c_str()); // const_cast for input directory //    char* odir = const_cast<char*>(output_dir.c_str()); // const_cast for output directory
    dim = ConfigVector.at(0); // space dimension of the problem (dim = 2 or 3);
    std::vector<unsigned int> kface_sequence; // special Faces sequence for the Kinetic module
/// Below the file names with the sparse PCC matrices which must already exit in the input_dir and have the same names (!)
/// They can be obtain by the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) based on the Neper output (*.tess file) (https://neper.info/) or PCC Structure Generator tool (https://github.com/PRISBteam/PCC_Structure_Generator) for plasticity problems
    string ssd0 = source_dir + "A0.txt"s, ssd1 = source_dir + "A1.txt"s, ssd2 = source_dir + "A2.txt"s, ssd3 = source_dir + "A3.txt"s, ssd4 = source_dir + "B1.txt"s, ssd5 = source_dir + "B2.txt"s, ssd6 = source_dir + "B3.txt"s,
            seeds = source_dir + "seeds.txt"s; //, NewSeeds = input_dir + "NewSeeds/NewSeeds.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str()));
    paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));
    paths.push_back(const_cast<char*>(seeds.c_str()));  //paths.push_back(const_cast<char*>(NewSeeds.c_str()));

    // Checking file names
    char* number_of_cells; string  ncells;
    // File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    if (is_file_exists(source_dir + "voro_Ncells.txt"s)) { ncells = source_dir + "voro_Ncells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }
    else if (is_file_exists(source_dir + "delau_Ncells.txt"s)) { ncells = source_dir + "delau_Ncells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }
    else { ncells = source_dir + "number_of_cells.txt"s; number_of_cells = const_cast<char*>(ncells.c_str()); }

    /// Useful ofstreams cleaning for data output
    // 1 LogFile ofstream
    Out_logfile_stream.open(output_dir + "log_file.txt"s, ios::trunc); // will close at the end of the main (!)
    Out_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    // Writer module global ofstreams //    OutJFile.open(output_dir + "TJp_random_theory.txt"s, ios::trunc); OutJFile.close(); //    OutJ2File.open(output_dir + "TJp_crystalline_theory.txt"s, ios::trunc); OutJ2File.close(); //    OutSFile.open(output_dir + "Sp_random_theory.txt"s, ios::trunc); OutSFile.close(); //    OutJ2File.open(output_dir + "Sp_crystalline_theory.txt"s, ios::trunc); OutS2File.close(); //    OutAgglStatistics.open(output_dir + "Stat_Agglomerations.txt"s, ios::trunc); OutAgglStatistics.close(); // Agglomerations statistics output

/// Output paths.vector to console and logfile out
    int npath = 0;
    cout << "_____________________________________________________________________________________" << endl;
    for (auto m : paths) cout <<"[" << npath++ << "]" << " paths:\t" << m << endl; for (auto p : paths) Out_logfile_stream <<"[" << npath++ << "]" << " paths:\t" << p << endl;

/// Principal variables
///_____________________________________________________________________________________

/// Read all the number of cells from file number_of_cells
// :: CellNumbs vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
    CellNumbs = VectorReader(number_of_cells); // VectorReader is a function from the PCC_SupportFunctions.h library; "number_of_cells" here if the PATH to file

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
/*
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
    */

    /// All vertex coordinates
    string VCpath_string = source_dir + "vertex_seeds.txt"s;
    char* VCpath = const_cast<char*>(VCpath_string.c_str());
    vertex_coordinates_vector.resize(CellNumbs.at(0), make_tuple(0,0,0));
    vertex_coordinates_vector = TuplesReader(VCpath);
//REPAIR    cout << "vertex_coordinates_vector size: " << vertex_coordinates_vector.size() << endl; //exit(0); exit(0);

    /// All grain volumes
    string GVpath_string = source_dir + "grain_volumes.txt"s;
    char* GVpath = const_cast<char*>(GVpath_string.c_str());
    grain_volumes_vector.resize(CellNumbs.at(3));
    grain_volumes_vector = dVectorReader(GVpath);

    /// All face areas
    string GBApath_string = source_dir + "face_areas.txt"s;
    char* GBApath = const_cast<char*>(GBApath_string.c_str());
    face_areas_vector.resize(CellNumbs.at(3));
    face_areas_vector = dVectorReader(GBApath);

/// Special feature for the 2D case ("grains" -> 2-cells; "faces" -> 1-cells; "triple junctions (lines) -> 0-cells" in 2D):
// Very important!
// It needs to make everything similar in the code for 2D and 3D cases with the same output of the PCC Generator tool (https://github.com/PRISBteam/Voronoi_PCC_Analyser) software (the same set of A0, A1..., B2 matrices)
    if (dim == 2) { CellNumbs.push_back(CellNumbs.at(2)); CellNumbs.at(2) = CellNumbs.at(1); CellNumbs.at(1) = CellNumbs.at(0); CellNumbs.at(0) = NULL; }

/// CellNumbs output
    cout << "=====================================================================================" << endl;  Out_logfile_stream << "==========================================================================================================================================================================" << endl;
    unsigned int t_length = 0;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl; for (int j : CellNumbs) Out_logfile_stream << t_length++ << "-cells #\t" << j << endl;
    //---------------------------------------------------------------------------
/// The PRIMAL DATA STRUCTURES and concepts - special (SFS) and ordinary (OFS) face sequences, defect state vectors (DSVs) and defect configurations (DCs)
    // Order of newly generated special faces | Process time
    std::vector <unsigned int> special_faces_sequence, ordinary_faces_sequence, current_sfaces_sequence, current_ofaces_sequence, crack_faces_sequence, current_cracks_sequence; // Variable sequences (in order of their generation) of special and ordinary Faces and Cracks
    vector<unsigned int> sub_sfaces_sequence, frac_sfaces_sequence;
    // State vector | Special faces IDs
    vector <unsigned int> State_sVector(CellNumbs.at(2), 0), current_State_sVector(CellNumbs.at(2), 0), State_cVector(CellNumbs.at(2), 0), current_State_cVector(CellNumbs.at(2), 0);
    std::vector <double> face_elastic_energies;
    cout << "=========================================================================" << endl << "\t\t\t\t\t[\tStart of the PCC Processing Tool\t]\t\t\t\t\t" << endl << "-------------------------------------------------------------------------" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl << "\t\t\t\t\t[\tStart of the PCC Processing Tool\t]\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

/// ==========================================================================================================================================
/// ================================================= THE LIST MODE STARTS HERE ==============================================================
/// ==========================================================================================================================================
/// In the LIST mode all the functions are calling one after another without loops
    double P_time = 0, K_time = 0, C_time = 0, W_time = 0;
    double mu_f = 0, sigm_f = 0; // only for L type simulations (!)
    vector<vector<int>> RW_series_vector; // only for L type simulations (!)

/// Clearance of files
//REPAIR    OfStreams_trancator();

/// THE CURRENT TASKS ARE PLACED HERE ///
/// ==========================================================================================================================================

// # 1 # Simple sfaces processing
#include "tasks/task_sFacesProcessing.cpp"

// # 2 # Processing design for Strips of Inclusions (with the parameters of average \mu and dispersion \sigma of the strip lengths distribution)
//#include "tasks/task_StripsProcessingDesign.cpp"

// # 3 # Macrocrack growth with multiple cracking simulations
//#include "tasks/task_macrocrack.cpp"

/// ==========================================================================================================================================

// For experiments

/// ==========================================================================================================================================


/// ===== Elapsing time ================>
    unsigned int end_time = clock();
    double fulltime = (double) end_time;
    cout << "-------------------------------------------------------------------------" << endl;     Out_logfile_stream << "-------------------------------------------------------------------------" << endl;
    cout << dim << "D " << "runtime is equal to  " << fulltime/pow(10.0,6.0) <<  "  seconds" << endl; Out_logfile_stream << dim << "D " << "runtime is equal to  " << fulltime/pow(10.0,6.0) <<  "  seconds" << endl;
    cout << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the VoroCAnalyser program\t]\t\t\t\t\t" << endl << "=========================================================================" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl << "\t\t\t\t\t[\tThe end of the VoroCAnalyser program\t]\t\t\t\t\t" << endl << "=========================================================================" << endl;

    // closing logfile ofstream
    Out_logfile_stream.close();

    return 0;
} /// The END of Main function

/// ================================== DEFINED IN MAIN FUNCTIONS ==================================///

/// ================== Initial configuration - reading and output =================>

/// Reading and Output of the configuration file
std::vector<double> config_reader_main(char* config, string &Subcomplex_type, string &Processing_type, string &Kinetic_type, string &input_folder, string &output_folder) {
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
                    line3_stream >> input_folder;
                } // input folder path input_dir = const_cast<char*>(input.c_str());
                else if (it == '$') {
                    stringstream line4_stream(line);
                    line4_stream >> output_folder;
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
    cout << "Input folder:\t" << input_folder << endl;
    cout << "Output folder:\t" << output_folder << endl;
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


/// OLD HEAP ///
// void eraseSubStr(std::string & mainStr, const std::string & toErase); // support (technical) function: erase the first occurrence of given substring from main string
// Still not working in Windows:(
// using std::__fs::filesystem::current_path; // to obtain current working directory
// string MainPath = current_path(); // MainPath variables
// eraseSubStr(MainPath, "cmake-build-debug"s); // to delete .../cmake-build-debug/ from the path
