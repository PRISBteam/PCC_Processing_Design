///******************************************************************************************************************************///
///*                                                                                                                            *///
///******************************************************************************************************************************///
///********************************** Dr Elijah Borodin, Manchester, UK, Spring 2022   ***********************************************///
///******************************************************************************************************************************///
///****************************************   DCCAnalyser(c) utility   ********************************************************///
///******************************************************************************************************************************///

/// Standard (STL) C++ libraries:
///------------------------------
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>
#include<numeric>
#include <string>
#include <map>
#include <algorithm>
///-------------------------------
/// Attached user defined C++ libraries:
///------------------------------------------
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/SymEigsSolver.h>
///------------------------------------------
using namespace std; //Standard namespace
using namespace Eigen;
using namespace Spectra;

// Triplets in the form T = T(i,j,value), where i and j element's indices in the corresponding dense matrix
typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type name - Tr
typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles name - SpMat

#include "src/DCC_SupportFunctions.h" // Various useful functions
#include "src/DCCProcessing/DCCProcessing.h" // Change element types of the DCC itself
#include "src/DCCKinetic/DCCKinetic.h" // Generate a process on the elements of the DCC without any changes in the complex
#include "src/DCCCharacterisation/StructureCharacterisation.h" // Characterisation of special substructures in the DCC

/// Declaration of FUNCTIONS, see the function bodies at the end of file.
std::vector<double> confCount(char* config, char* type, char* Kinetic_type, string& input_folder, string& output_folder); // Read and output the initial configuration
bool ProcessingON(char* config); // Check the Processing Module Status
bool CharacterisationON(char* config); // Check the Structure Characterisation Module Status
bool KineticON(char* config); // Check the Kinetic Module Status
bool SIMULATION_MODE(char* config); // Check the SIMULATION MODE
std::string get_working_path(); // Get the current working directory (containing Main.cpp)
void eraseSubStr(std::string & mainStr, const std::string & toErase); //Erase First Occurrence of given  substring from main string

///*.....................................................................    Main    .............................................*///
int main() {
/** Two functions (dependent on the problem's spatial dimension - 2D or 3D) are launching here with the arguments of the paths for
 * adjacency AN, AE, AF, (AG) and incidence (boundary operators) BEN, BFE, (BGF) matrices. **/

    cout << "=============================================" << endl << "Start of the DCC Processing" << endl << "-------------------------------------------------------------" << endl;
/// File path to the configuration profile
    using std::__fs::filesystem::current_path; //to obtain current working directory
    string MainPath = current_path();
    eraseSubStr(MainPath, "cmake-build-debug"s);

    string  config = MainPath + "config.txt"s; char* confpath = const_cast<char*>(config.c_str());

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
// The source directory and simulation type from file config.txt
    char simulation_type; // 'R', 'S', 'W', 'F' or 'I', 'E' :: This char define the process type: 'R' for Random, 'S' for maximum configuration Entropy production, 'F' for the 3D one-layer film, 'I' for the Ising-like model, 'E' for experimental data obtained by EBSD
    char Kinetic_type; // 'W' or 'I' :: This char define the Kinetic process type: 'W' for Wear, 'I' for Ising
    string input_folder, output_folder; // input and output folders from file config
    vector<bool> SChar_config; // Characterisation module clonfiguration
    std:vector<double> ConfigVector = confCount(confpath, &simulation_type, &Kinetic_type, input_folder, output_folder);
    int dim = ConfigVector.at(0); // Space dimension of the problem (dim);

/// Below the file names with the sparse DCC matrices must be defined
    string ssd0 = input_folder + "A0.txt"s, ssd1 = input_folder + "A1.txt"s, ssd2 = input_folder + "A2.txt"s, ssd3 = input_folder + "A3.txt"s, ssd4 = input_folder + "B1.txt"s, ssd5 = input_folder + "B2.txt"s, ssd6 = input_folder + "B3.txt"s,
            seeds = input_folder + "seeds.txt"s, NewSeeds = input_folder + "NewSeeds/NewSeeds.txt"s;
//The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    vector<char*> paths;
    char* odir = const_cast<char*>(output_folder.c_str());
    char* indir = const_cast<char*>(input_folder.c_str());
    paths.push_back(const_cast<char*>(ssd0.c_str())); paths.push_back(const_cast<char*>(ssd1.c_str())); paths.push_back(const_cast<char*>(ssd2.c_str())); paths.push_back(const_cast<char*>(ssd3.c_str()));
    paths.push_back(const_cast<char*>(ssd4.c_str())); paths.push_back(const_cast<char*>(ssd5.c_str())); paths.push_back(const_cast<char*>(ssd6.c_str()));
    paths.push_back(const_cast<char*>(seeds.c_str())); paths.push_back(const_cast<char*>(NewSeeds.c_str()));
// File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    string  ncells = input_folder + "number_of_cells.txt"s; char* number_of_cells = const_cast<char*>(ncells.c_str());

/// Principal variables
///_____________________________________________________________________________________
    // Read all the number of cells from file number_of_cells
    // :: vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
    std::vector<unsigned int> CellNumbs = VectorReader(number_of_cells);
    // State vector | Special faces IDs
    std::vector <unsigned int> State_Vector(CellNumbs.at(2), 0), current_State_Vector(CellNumbs.at(2), 0);
    // Order of newly generated special faces | Process time
    std::vector <unsigned int> special_faces_sequence, current_sfaces_sequence; // Sequence (in order of their generation time) of special Faces
    std::vector <unsigned int> crack_faces_sequence, crack_current_sequence;
    // The number of special Face (2-cells) types
    int number_of_types = ConfigVector.at(1);
    // MAX fraction of Faces | Calculation limit
    double max_sFaces_fraction = ConfigVector.at(3);

/// ===========================================================================================================================================
/// ========================================== SCIENTIFIC MODE STARTS HERE ===================================================================
/// ===========================================================================================================================================
    if (SIMULATION_MODE(confpath)) {
        ///  ====== Analytical solutions ========>
       // TJsAnalytics(1000, odir); // A function from TJsLab.h
        /// ====== Processing ========>
        cout << "START of DCC Processing Module" << endl;
        if (ProcessingON(confpath))
            DCC_Processing3D(State_Vector, special_faces_sequence, simulation_type, max_sFaces_fraction, number_of_types, CellNumbs, paths);
        //for (auto sfe : special_faces_sequence) cout << sfe << "\t"; cout << endl;
//        special_faces_sequence = VectorReader("/Users/user/Dropbox/OFFICE/NEPER/resultsJune2022/1kCells/I/SpecialGrainBoundaries.txt"); //all Faces

        //string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s;
        //string output_TJs_dir = output_folder + TJs_output_filename, output_Entropy_dir = output_folder + Entropy_output_filename;
        string output_TJs_dir = output_folder + "TJsLab_TJsTypes.txt"s, output_Entropy_dir = output_folder + "TJsLab_ConTJsEntropy.txt"s;
        char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

        /// Creation of the two files with TJs fraction and Configuration entropies output
        ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::trunc); OutTJsFile.close();
        ofstream OutSFile; OutSFile.open(cEntropy_dir, ios::trunc); OutSFile.close();

        /// Loop over special_faces_sequence with the step = Face_fraction/ number_of_steps
    double Face_fraction = 0.05, Crack_fraction = 0.005; //initial fraction
    long number_of_steps = 10, number_of_csteps = 10;
    double dp = 0, df = 0, max_cFaces_fraction = 0.3; //increments (dp - special Faces fraction, df - cracked faces fraction) and max fraction of fractured (cracked) faces

        /// Output to file Special Face Laplacian eigenvalues
        // Special Face Laplacian
        ofstream OutElCondfile; OutElCondfile.open(odir + "Face_Conductivity.txt"s, ios::trunc);
        OutElCondfile << "\tFraction of special Faces\t" << " " << "\tFraction of cracks\t" << " " << "\tSpecial Face conductivity\t" << " " << "\tNumber of Special Face components\t"
                      << " " << "\tSpecial Face median entropy\t" << " " << "\tSpecial Face skrew entropy\t" << " " << "\tSpecial Face informativeness\t" << " " << "\tCracked Face conductivity\t" << " " << "\tNumber of Cracked Face components\t" << " " << "\tCracked Face median entropy\t" << " " << "\tSpecial Face skrew entropy\t" << " " << "\tCracked Face informativeness\t" << endl;
        OutElCondfile.close();
        ofstream OutFLfile; OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
        OutFLfile << "Laplacian Matrix of all Special Faces" << endl;
        OutFLfile.close();

        OutElCondfile.open(odir + "Face_Conductivity.txt"s, ios::app);
        OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::app);

        /// The first loop over all Face_fraction with the step dp
        for(long i = 0; i < number_of_steps; ++i) { // loop over all Face_fraction
         Face_fraction += dp;
         unsigned int special_Faces_numb = Face_fraction * CellNumbs.at(2);

         unsigned int itu = 0;
         current_sfaces_sequence.clear();
             for (unsigned int sfe : special_faces_sequence) {
                     if (itu <= special_Faces_numb) current_sfaces_sequence.push_back(sfe);
                     ++itu;
                 } // for (auto sfe : current_sfaces_sequence) cout << sfe << "\t"; cout << endl;

             unsigned int itr = 0;
             fill(current_State_Vector.begin(), current_State_Vector.end(), 0);
             for(auto sv : State_Vector) {
                 if (sv != 0 && find(current_sfaces_sequence.begin(),current_sfaces_sequence.end(),itr) != current_sfaces_sequence.end())
                 current_State_Vector.at(itr) = sv;
                 itr++;
             } //  for (auto sfe : current_State_Vector) cout << sfe ; cout << endl;

            /// ====== Fracture Kinetic ========>
            //cout << "START of DCC Kinetic Module" << endl;
            if (KineticON(confpath))
                crack_faces_sequence = DCC_Kinetic(Kinetic_type, current_sfaces_sequence, paths, indir, odir);
//REPAIR        for (auto cn : crack_sequence)  cout << cn << "\t"; cout << endl; //        exit(19);

            /// The second loop over all crack_Faces_numb
            for(long y = 0; y < number_of_csteps; ++y) { // loop over all Face_crcks
                Crack_fraction += df;

                unsigned int crack_Faces_numb = Crack_fraction * CellNumbs.at(2);
                unsigned int itc = 0;
                crack_current_sequence.clear();
                for (unsigned int sfe : crack_faces_sequence) {
                    if (itc <= crack_Faces_numb) crack_current_sequence.push_back(sfe);
                    ++itc;
                }
// REPAIR       for (auto sfe : crack_current_sequence) cout << sfe << "\t"; cout << endl;

                /// ====== Structure Characterisation module ========>
//                if (Face_fraction == 0) cout << "START of DCC Structure Characterisation Module" << endl;
                if (CharacterisationON(confpath)) {
                    DCC_StructureCharacterisation(current_State_Vector, current_sfaces_sequence, crack_current_sequence,
                                                  ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile);
                } else break;
            /// Increment of crack fraction
                df = max_cFaces_fraction / (double) number_of_csteps; /// Crack fraction increment - 0.95 is an arbitrary limit close to 1 !!!
            } // for (number_of_csteps)
            Crack_fraction = 0; // Start each dp iteration with Zero crack fraction

            /// Increment of special faces fraction
            dp = max_sFaces_fraction / (double) number_of_steps; // Face fraction increment
        } // for (number_of_steps)

        OutFLfile.close();
        OutElCondfile.close();

    }// SIMULATION MODE if

/// ==========================================================================================================================================
/// ================================================= LIST MODE STARTS HERE ===================================================================
/// ===========================================================================================================================================
    else {
/// I: DCC_Processing module
        if (ProcessingON(confpath)) { // if DCC_Processing is SWITCH ON in the config.txt file
            DCC_Processing3D(State_Vector, special_faces_sequence, simulation_type, max_sFaces_fraction, number_of_types, CellNumbs, paths);
            //cout << "new:" <<endl; for (auto itd : special_faces_sequence) cout << itd << endl;
        } // if(ProcessingON(confpath))

/// III: DCC_Kinetic module
        if (KineticON(confpath)) { // if DCC_Kinetic is SWITCH ON in the config.txt file
            DCC_Kinetic(Kinetic_type, special_faces_sequence, paths, indir, odir);
            }// if(KineticON(confpath))
/// II: DCC_Characterisation module
        if (CharacterisationON(confpath)) {
            cout << "START of DCC Structure Characterisation Module" << endl;
            ofstream OutFLfile; ofstream OutElCondfile;
            DCC_StructureCharacterisation(State_Vector, special_faces_sequence, crack_faces_sequence, ConfigVector, CellNumbs, paths, odir, OutFLfile, OutElCondfile);
        }// if(CharacterisationON(confpath))
    }// SIMULATION MODE else

/// ===== Elapsing time ================>
    unsigned int end_time = clock();
    double fulltime = (double) end_time/ 60000.0;
    cout << "HAGBsProbability " << dim << "D " << "runtime is equal to  " << fulltime <<  "  seconds" << endl;

    cout << "-------------------------------------------------------------" << endl << "The end of the VoroCAnalyser program" << endl << "=============================================" << endl ;

    return 0;
} /// The end of Main function


/// ================================== Related functions ==================================///

///Initial configuration - reading and output
bool ProcessingON(char* config) {
    std::string line;
    ifstream inConf(config);
    bool isProcessingON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
//            cout << line << endl;
            if (!line.compare("DCC_Processing SWITCHED ON"s)) {
                isProcessingON = 1;
                cout << "DCC_Processing SWITCHED ON"s << endl;
                return isProcessingON;
            }
    } else cout << "ProcessingON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "DCC_Processing SWITCHED OFF"s << endl;
    return isProcessingON;
}

bool KineticON(char* config) {
    std::string line;
    ifstream inConf(config);
    bool isKineticON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
//            cout << line << endl;
            if (!line.compare("DCC_Kinetic SWITCHED ON"s)) {
                isKineticON = 1;
                cout << "DCC_Kinetic SWITCHED ON"s << endl;
                return isKineticON;
            }
    } else cout << "KineticON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "DCC_Kinetic SWITCHED OFF"s << endl;
    return isKineticON;
}

bool CharacterisationON(char* config) {
    std::string line;
    ifstream inConf(config);
    bool isCharacterisationON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
//            cout << line << endl;
            if (!line.compare("DCC_Characterisation SWITCHED ON"s)) {
                isCharacterisationON = 1;
      //          cout << "DCC_StructureCharacterisation SWITCHED ON"s << endl;
                return isCharacterisationON;
            }
    } else cout << "isCharacterisationON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "DCC_StructureCharacterisation SWITCHED OFF"s << endl;
    return isCharacterisationON;
}

bool SIMULATION_MODE(char* config) {
    std::string line;
    ifstream inConf(config);
    bool isCharacterisationON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
//            cout << line << endl;
            if (!line.compare("SIMULATION MODE: #LOOP"s)) {
                isCharacterisationON = 1;
                cout << "LOOP SIMULATION MODE"s << endl;
                return isCharacterisationON;
            }
    } else cout << "SIMULATION_MODE() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "LIST SIMULATION MODE"s << endl;
    return isCharacterisationON;
}


/// Reading and Output of the configuration file
std::vector<double> confCount(char* config, char* type, char* Kinetic_type, string &input_folder, string &output_folder) {
    std::string line, input;
    double p_max = 0, output_step = 0;
    vector<double> res;

    ifstream inConf(config);
    if (inConf) { //If the file was successfully open, then
        while (getline(inConf, line, '\n')) {
            // Number of types
            for (auto it: line)
                if (it == '@')  res.push_back(line.at(0)-'0'); // dimension
                else if (it == '&') *type = line.at(0); // simulation type
                else if (it =='`') *Kinetic_type = line.at(0); // simulation Kinetic type
                else if (it == '!')  res.push_back(line.at(0)-'0'); // number of special Face types
                else if (it == '?')  { stringstream line1_stream(line); line1_stream >> p_max; res.push_back(p_max);} // MAX fraction of Faces (calculation limit)
                else if (it == '^')  { stringstream line2_stream(line); line2_stream >> output_step; res.push_back(output_step);} // the output step (number of new special steps) between the outputs
                else if (it == '~')  { stringstream line3_stream(line); line3_stream >> input_folder; } // input folder path input_folder = const_cast<char*>(input.c_str());
                else if (it == '$')  { stringstream line4_stream(line); line4_stream >> output_folder; } // output folder path input_folder = const_cast<char*>(input.c_str());

                    //        if(line.at(1) == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
                    //        if(line.at(1) == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
                else if (it == '#')  res.push_back(1); // 1 and # means accept - the parameter will be calculated
                else if (it == '%')  res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
        }
    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong
    res.size();
    cout << "The problem dimension that is the maximum dimension k_max of k-Cells:\t | "s << res.at(0) << endl;
    cout << "Calculation type ('R', 'W', 'S', 'F', 'I' or 'E'):\t\t\t\t\t\t | "s << *type << endl;
    cout << "Kinetic type ('W' or 'I'):\t\t\t\t\t\t\t\t\t\t\t\t | "s << *Kinetic_type << endl;
    cout << "The number of special Face (2-cells) types:\t\t\t\t\t\t\t\t | " << res.at(1) << endl;
    cout << "MAX fraction of Faces (calculation limit): \t\t\t\t\t\t\t\t | " << res.at(2) << endl;
    cout << "A number of new special Faces converted between the outputs: \t\t\t | " << res.at(3) << endl;
    cout << endl;
    cout << "Input folder:\t" << input_folder << endl;
    cout << "Output folder:\t" << output_folder << endl;
    cout << endl;
    cout << "Nodes types statistics, indices and configuration entropy:     "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges types statistics, indices and configuration entropy:     "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Faces types statistics and structural indices:                 "; if (res.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain types statistics, indices and configuration entropy:     "; if (res.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Node Laplacian with its spectrum for the Nodes network:       "; if (res.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edge Laplacian with its spectrum for the Edges network:       "; if (res.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Face  Laplacian with its spectrum for the Faces network:       "; if (res.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain Laplacian with its spectrum for the Grains network:      "; if (res.at(11) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Tutte polynomial for the special network:                      "; if (res.at(12) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;

    return res;
}

//Erase First Occurrence of given  substring from main string
void eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    // Search for the substring in string
    size_t pos = mainStr.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

/// Archive
/*void Manual_input () {
 do { ///Manual user input of the space dimension value (dim)
    cout << " Please, input the dimension of the problem: 3 for 3D and 2 for 2D "s << endl;
    cin >> dim;
    cout << "The problem dimension is\t" << dim << endl;
    if (dim != 2 && dim != 3) cout << "Input Error: Please retype 2 or 3" << endl;
    }while (dim != 2 && dim != 3); */
/* do { ///Manual user input of the simulation type
    cout << " Please, input the symbol of particular simulation type of the problem: E (experiment), R (rndom), S (entropy maximisation) or I (Ising model):"s << endl;
    cin >> simulation_type;
    cout <<"The simulation type is\t"<< simulation_type << endl;
    if (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I') cout << "Input Error: Please retype 'E', 'R', 'S' or 'I' for the specific simulation type" << endl;
    }while (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I');  */
/* /// Manual user input of the DCC files folder path
    cout << " Please, input the name of the folder where the DCC source files are (like 1k3cells/ ):"s << endl;
    cin >> problem_folder_path; */
/* /// Manual user input of the simulation results output folder path
    cout << " Please, input the simulation results output folder path (like test/ ):"s << endl;
    cin >> results_output_folder_path;
return;
} */


/*
///=============================================================================================================================================////
///=================================== Characterisation module =============================>>

if( OCellAmount % output_step == 0 && ordinary_faces_fraction > 0.05) { // Periodicity of characterisation output
unsigned int SAM_size = Structure_Characterisation(numerator, CellNumbs, SFaces_Triplet_list, SpecialCellMap, special_faces_sequence, FES, odir, special_faces_fraction, configuration);

cout << "Fraction of special faces:" << 1.0 - ordinary_faces_fraction << ";\t Size of S-Face Adjacency matrix \t" << SAM_size << endl;
} // End of analysis and output iterator ( IF: iterator % X == 0 )
*/