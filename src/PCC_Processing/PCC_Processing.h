///================================ PCC Processing module ===================================================================================///
///=========================================================================================================================================///
///* The interface using functions from Processing_Functions.h C++ library to generate quasi-random or non-random processes of changes    *///
///* in the elements of the pre-constructed polyhedral cell complex (PCC) 'M' with the whole set of incidence and adjacency martices     *///
///* ---------------------------------------------------------------------------------------------------------------------------------- *///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design tool (PDT code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1б 2022-2023 years   *///
///===================================================================================================================================///

/// Attached user-defined C++ libraries:
// Here a Processing_Functions.h library of C++ functions for advanced random and non-random generations of special chains of PCC elements
///-------------------------------------
#include "Processing_Functions.h"
///-------------------------------------

string SIMULATION_MODE(char* config); // Check the SIMULATION MODE type in the config.txt file

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
// The source directory and simulation type from file config_processing.txt
string S_type; // 'P' or 'H' :: This char define the subsection type: 'P' for the whole Plane cut, 'H' for the half-plane cut like a crack


/// User-defined function definitions
bool ProcessingON(char* config, bool time_step_one); // Check the Processing module status (On/Off) in the config.txt file

///* ========================================================= PCC PROCESSING FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///

// ARCHIVE sface_design :: A list of special_faces_sequence as an output of the module. // It can be any sequence specified in the "processing_task" file, but normally the first line here is the random case (zero-design), then following S_max (1-sequence) and Sd_min (Smin), and then all the designs in between S_max and S_min
// material :: quadruple points, triple junctions, grain boundaries, and grains
// tessellation :: points, lines, polygones, and polytopes
// PCC :: nodes, edges, faces, and volumes
// IDEA special_cells_design - vector of vectors containing :: (1) snodes_sequence, (2) sedges_sequence, (3) sfaces_sequence, and (4) svolumes_sequence
/*!
 * @param special_faces_sequence // base object - a sequence of generated 2-cell numbers
 * @param State_sVector
 *  @param P_type
 * @return special_cells_design
 */

/// Off-streams :
ofstream OutFLfile, OutJFile, OutJ2File, OutSFile, OutS2File, OutPowersADistributions, OutAvLengthsADistributions, OutAgglStatistics;

/// ConfigVector contains all the control variables of the program readed from the config.txt
vector<double> ConfigVector = confCount(confpath, S_type, P_type, input_dir, output_dir);

// MAX fraction of special Faces | Calculation limit
double max_sFaces_fraction = ConfigVector.at(2); // maximum fractions of special (defect) faces (2-cells)
 
std::vector<vector<unsigned int>> DCC_Processing(std::vector<unsigned int> &special_faces_sequence, vector<unsigned int> &State_sVector, const string P_type) {
// CellNumbs :: vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
// Maximal fraction (max_sFaces_fraction) for simulation loop max_sFaces_fraction = [0,1]
// State_Vector in the form : [Element index] - > [Type] = kind of a code related to the microstructure of PCC

// Type of the Processing tool from config.txt (const_cast<char*> for processing type)
    char* stype = const_cast<char*>(P_type.c_str());
    char simulation_type = *stype;

/// Cases for Processing types
    if (simulation_type == 'R') { ///  Random generation case
        Processing_Random(State_sVector, special_faces_sequence, max_sFaces_fraction);
        special_cells_design.push_back(special_faces_sequence);
    } ///End of 'R' type simulations

    else if (simulation_type == 'S') { // Maximum entropy production
        Processing_maxEntropy(State_sVector, special_faces_sequence);
        special_cells_design.push_back(special_faces_sequence);
        cout << "Processing_maxEntropy - successfully finished" << endl;

    } ///End of 'S' type simulations

/*
    else if (simulation_type == 'S') { // Maximum entropy production
        Processing_Random(State_sVector, special_faces_sequence, max_sFaces_fraction);
        special_cells_design.push_back(special_faces_sequence);
        cout << "Processing_Random - successfully finished" << endl;

        std::fill(State_sVector.begin(), State_sVector.end(), 0);
        special_faces_sequence.clear();
        Processing_maxEntropy(State_sVector, special_faces_sequence);
        special_cells_design.push_back(special_faces_sequence);
        cout << "Processing_maxEntropy - successfully finished" << endl;

    } ///End of 'R+S' type simulations

    else if (simulation_type == 'S') { // Maximum entropy production
        Processing_Random(State_sVector, special_faces_sequence, max_sFaces_fraction);
        special_cells_design.push_back(special_faces_sequence);
        cout << "Processing_Random - successfully finished" << endl;
        Processing_maxEntropy(State_sVector, special_faces_sequence);
        cout << "Processing_maxEntropy - successfully finished" << endl;
        special_cells_design.push_back(special_faces_sequence);
//        Processing_minEntropy(State_sVector, special_faces_sequence);
//        cout << "Processing_minEntropy - successfully finished" << endl;
//        sface_design.push_back(special_faces_sequence);

    } ///End of 'S' type simulations

    else if (simulation_type == 'S') { // Maximum entropy production
        Processing_Random(State_sVector, special_faces_sequence, max_sFaces_fraction);
        cout << "Processing_Random - successfully finished" << endl;
        special_cells_design.push_back(special_faces_sequence);
        Processing_maxEntropy(State_sVector, special_faces_sequence);
        cout << "Processing_maxEntropy - successfully finished" << endl;
        special_cells_design.push_back(special_faces_sequence);
//        Processing_minEntropy(State_sVector, special_faces_sequence);
//        cout << "Processing_minEntropy - successfully finished" << endl;
//        sface_design.push_back(special_faces_sequence);

        for (int it = 3; it < design_number; ++it) {
            //ip_index = Smax - Smin / (design_number - 1.0);
            //int Processing_ipIndex(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, int index_type, double ip_index) {
        }
    } ///End of 'S' type simulations
    */

    else cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl;

    return special_cells_design;
} /// The end of HAGBsProbability3D()

std::vector<vector<unsigned int>> DCC_Processing(std::vector<unsigned int>  &special_faces_sequence, vector<unsigned int> &State_sVector, const string P_type, double &mu_f, double &sigm_f, vector<vector<int>> &RW_series_vector) {
/// sface_design :: A list of special_faces_sequence as an output of the module. The first line here is always the random case (zero-design), then following Smax (1-sequence) and Smin (2-sequence), and then all the designs in between Smax and Smin
// CellNumbs :: vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
// Maximal fraction (max_sFaces_fraction) for simulation loop max_sFaces_fraction = [0,1]
// State_Vector in the form : [Element index] - > [Type]
///vector<unsigned int> State_sVector(CellNumbs.at(2),0); // State vector filling with zeros

/// Used functions
    std::vector<double>  Log_normal_distribution(double &mu_f, double &sigm_f, int baskets);

/// Type of the Processing tool from config.txt (const_cast for processing type; const_cast for output directory)
    char* stype = const_cast<char*>(P_type.c_str());
    char simulation_type = *stype;

    //// Cases for Processing types ////

    if (simulation_type == 'R') { //  Random generation case
        Processing_Random(State_sVector, special_faces_sequence, max_sFaces_fraction);
        special_cells_design.push_back(special_faces_sequence);

    } ///End of 'R' type simulations

    else if (simulation_type == 'L') { //  Random lengthy strips generation case
        // Probability density function
//REPAIR        cout << "ATTENTION: mu and sigm" << endl; cout << mu_f << "  " << sigm_f << endl;
        std::vector<double> strip_lenght_distribution = Log_normal_distribution(mu_f, sigm_f,30);

        // From fractions to number of Faces:
        face_strip_distribution.clear(); // clearing file
        for (auto  itr = strip_lenght_distribution.begin(); itr != strip_lenght_distribution.end(); ++itr)
//        face_strip_distribution.push_back(*itr); //           face_strip_distribution.push_back(max_sFaces_fraction * CellNumbs.at(2) * (*itr));
            face_strip_distribution.push_back(max_sFaces_fraction * CellNumbs.at(2) * (*itr)/ (distance(strip_lenght_distribution.begin(), itr) + 1.0));

        /// Generation of special Faces structure
        RW_series_vector = RStrips_Distribution(face_strip_distribution, State_sVector, special_faces_sequence);
        special_cells_design.push_back(special_faces_sequence);

    } ///End of 'L' type simulations

    else cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl;

    return special_cells_design;
} /// The end of HAGBsProbability3D()

/// * Related functions * ///
// 4
bool ProcessingON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isProcessingON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("DCC_Processing SWITCHED ON"s)) {
                isProcessingON = 1;
                if (time_step_one == 1) cout << "ON    | DCC_Processing"s << endl;
                return isProcessingON;
            }
    } else cout << "ProcessingON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF    | DCC_Processing"s << endl;
    return isProcessingON;
}

/// Heap
// bool time_step_one = 1; // (technical support variable) special ID for the first step of iteration MAX fraction of Faces

/*     /// Ordinary face sequence (!!works!!)
   ordinary_faces_sequence.clear();
   for (auto  itr = State_sVector.begin(); itr != State_sVector.end(); ++itr)
       if (*itr == 0) ordinary_faces_sequence.push_back(distance(State_sVector.begin(), itr));

*
*
* else if (*stype == 'D') { // DDRX recrystalisation process
       Processing_DDRX(State_sVector, special_faces_sequence);
   } /// End of 'D' type simulations

   else if (*stype == 'X') { // structural index-based calculation
       /// Loop over all the index values
       enum index_types { Smin = 0, ip_ind, p_star, labda_ind};
       /// Actual value of the index !!! (move after to config.txt!)
       int index_type = Smin;

       int ip_stepsnumb = 1;
       double ip_index = 1.0, dip = 1.0;
       /// Informativeness
       if (index_type == ip_ind) {
           int ip_stepsnumb = 100;
           double ip_index = -1, dip = 2.0 / ip_stepsnumb; // ip index can changes from -1 to 1
       }

       for (int index = 0; index < ip_stepsnumb; ++index) {
            Processing_ipIndex(State_sVector, special_faces_sequence, index_type, ip_index);

           ip_index += dip; // index requested value increasing
       } // end for ( index < ip_stepnumb)
   } /// End of 'X' (ip index) type of simulations

  else if (*stype == 'E') { // Experimental data
       Processing_ExperimentalData(State_Vector, special_faces_sequence, max_sFaces_fraction, number_of_types, CellNumbs, paths);

       /// ====== Reading from file instead of DCC_Processing3D ( ) if it is off ============>
       string SFS_path = in + "Smax/SpecialGrainBoundaries.txt"s;
       cout << "-----------------------------------------------------------------------------------------------------------------------------"s << endl;
       cout << "Warning!!: special_faces_sequence successfully loaded from file:\t"s << SFS_path << " because the Processing is OFF"s << endl;
       cout << "------------------------------------------------------------------------------------------------------------------------------"s << endl;
       char* SFS_dir = const_cast<char*>(SFS_path.c_str());
       special_faces_sequence = VectorReader(SFS_dir); //all Faces
       for (auto itr: special_faces_sequence) State_Vector.at(itr) = 1;
   }
*/