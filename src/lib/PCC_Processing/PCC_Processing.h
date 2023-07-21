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
///-----------------------------------------------------

//#include "PCC_Processing.h"

 #include "functions/Processing_Imposed_functions.h"
 #include "functions/Processing_Assignment_functions.h"
 #include "functions/Processing_Induced_functions.h"
///-----------------------------------------------------

///* ========================================================= PCC PROCESSING FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
// special_cells_design - vector of vectors containing :: (1) snodes_sequence, (2) sedges_sequence, (3) sfaces_sequence, and (4) svolumes_sequence
/*! processing.ini
 * @param P_type
 * @param Configuration_State [not mandatory - only is there is some "non-zero" initial state]
 * @return cells_design
 */
CellsDesign PCC_Processing(std::vector<std::vector<int>> &Configuration_State,std::vector<std::vector<int>> &Configuration_cState) {
// CellNumbs :: vector components: [0] - Nodes number, [1] - Edges number, [2] - Faces number, [3] - Grains number
// Maximal fraction (max_sFaces_fraction) for simulation loop max_sFaces_fraction \in [0,1]
/// Main output of the module
CellsDesign CD;

/// The PRIMAL DATA STRUCTURES and concepts - special (SFS) and ordinary (OFS) face sequences, defect state vectors (DSVs) and defect configurations (DCs)
// Processing vector variables
//special_n_Sequence, special_e_Sequence, special_f_Sequence, special_p_Sequence & special_cells_design
// State_p_vector, State_f_vector, State_e_vector, State_n_vector & Configuration_State
///
std::vector <unsigned int> special_x_sequence, current_sx_sequence; // variable sequences (in order of their generation) of special and ordinary k-cells
std::vector <unsigned int> special_c_sequence; // variable sequences (in order of their generation) of "fractured" k-cells as the result of a kinetic process
// State vector | Special faces IDs
//////vector <unsigned int> State_sVector(CellNumbs.at(2), 0), current_State_sVector(CellNumbs.at(2), 0), State_cVector(CellNumbs.at(2), 0), current_State_cVector(CellNumbs.at(2), 0);
cout << "=========================================================================" << endl << "\t\t\t\t\t[\tStart of the PCC Processing Tool\t]\t\t\t\t\t" << endl << "-------------------------------------------------------------------------" << endl;
Out_logfile_stream << "==============================================================================================================================================================" << endl << "\t\t\t\t\t[\tStart of the PCC Processing Tool\t]\t\t\t\t\t" << endl << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

/// Read configuration from processing.ini file :: the number of special cell types and calculation parameters.
double mu = 1.0, sigma = 0.0; // mean and dispersion for a distribution
std::vector<vector<double>> max_fractions_vectors(4), max_cfractions_vectors(4); // [0][..] - nodes, [1][..] - edges, [2][..] - faces, [3][..] - polyhedrons
std::vector<string> ptype_vector(4), ctype_vector(4); // vector of strings corresponding to the processing types from file processing.ini :
std::vector<double> pindex_vector(4); // vector of indeces of the Processing module: S(p_index = 1) - Smax, S(p_index = 0) - Smin, I(p_index = x.x) - Index mode
// in both vectors: [0] - nodes, [1] - edges, [2] - faces, [3] - polyhedrons
std::vector<string> sequence_source_paths(4); // k-sequence pathes for the reading them from file(s) instead of the new Processing routes
config_reader_processing(source_path, sequence_source_paths, max_fractions_vectors, max_cfractions_vectors, mu, sigma, ptype_vector, ctype_vector, pindex_vector, Out_logfile_stream); // void

/// Cases for Processing types
// cell type k :: 0 - nodes, 1 - edges, 2 - faces, 3 -polyhedrons - must coincide with the indexing of the CellNumbs.at(cell_type) vector (!)
/// here dim !
for (int cell_type = (3 - (dim - 3)); cell_type >= 0; --cell_type) { /// loop over all types of k-Cells in the complex
// both for 2D and 3D cases
    special_x_sequence = {0}; special_c_sequence = {0};

    if (ptype_vector.at(cell_type + (3 - 3)) == "R" && max_fractions_vectors.at(cell_type + (3 - 3)).size() > 0) { //  Random generation case
        cout << "Random processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "Random processing in operation: cell_type : "s << cell_type  + (3 - 3) << endl;
        special_x_sequence = Processing_Random(cell_type, Configuration_State, max_fractions_vectors);
    } // End of 'R' type simulations

    else if (ptype_vector.at(cell_type + (3 - 3)) == "F" && max_fractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
//        if (pindex_vector.at(cell_type) == 0) {
// processing index :: 0 - direct special faces assignment;  1 - crystallographic ; 2 - configurational TJs-based entropy (deviatoric);
            cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
            Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
            special_x_sequence = Processing_maxFunctional(2, Configuration_State, max_fractions_vectors, pindex_vector.at(2));
            // cell type = 2 -> faces
//        } else if (pindex_vector.at(cell_type) == 1) {
//            cout << "MaxFunctional = Configuration Entropy processing in operation: cell_type : "s << cell_type << endl;
//            special_x_sequence =  Processing_maxFunctional(cell_type, Configuration_State, max_fractions_vectors, pindex_vector.at(cell_type));
//        }
    } // End of 'F' type simulations (elseif)
    else if (ptype_vector.at(cell_type + (3 - 3)) == "D" && max_fractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
        cout << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
    if (max_fractions_vectors.at(cell_type).size() > 0)
        special_x_sequence = Processing_minConfEntropy(2, Configuration_State, max_fractions_vectors, pindex_vector.at(2));

    } // End of 'D' [S min] type simulations (elseif)
    else if (ptype_vector.at(cell_type + (3 - 3)) == "Cm" && max_fractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
        cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
        if (max_fractions_vectors.at(cell_type).size() > 0)
            special_x_sequence = Processing_maxF_crystallographic(2, Configuration_State, max_fractions_vectors, pindex_vector.at(2));
        // cell type = 2 -> faces
    } // End of 'Cm' type simulations (elseif)

    else if (ptype_vector.at(cell_type + (3 - 3)) == "Cr" && max_fractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
        cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
        if (max_fractions_vectors.at(cell_type).size() > 0)
            special_x_sequence = Processing_maxP_crystallographic(2, Configuration_State, max_fractions_vectors, pindex_vector.at(2));
///        special_x_sequence = Processing_Random_crystallographic(2, Configuration_State, max_fractions_vectors, pindex_vector.at(2));
    } // End of 'Cr' type simulations (elseif)

    else if (ptype_vector.at(cell_type + (3 - 3)) == "S") {
        char* kseq_sourcepath = const_cast<char*>(sequence_source_paths.at(cell_type + (3 - 3)).c_str());
            special_x_sequence = VectorIReader(kseq_sourcepath);

        /// (!) Output +1 like in Neper, so he numbers should be modified back as -1
        for (auto it = special_x_sequence.begin(); it != special_x_sequence.end(); ++it)
            special_x_sequence.at(distance(special_x_sequence.begin(), it)) = *it - 1;
/// Update of the corresponding Configuration State vector
        Configuration_State[cell_type + (3 - 3)].clear();
        std::vector<int> State_vector(CellNumbs.at(cell_type + (3 - 3)), 0);
        for (unsigned int k_cell : special_x_sequence) // fill state vector from cells_sequence
            State_vector.at(k_cell) = 1;

        for (int var : State_vector)
         Configuration_State[cell_type + (3 - 3)].push_back(var);
    } // End of 'S' [reading from file] type simulations (elseif)
//    else if (ctype_vector.at(cell_type + (3 - 3)) == "Km" && max_cfractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
    else cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'F' or 'I'..!" << endl;

    if (ctype_vector.at(cell_type + (3 - 3)) == "Km" && max_cfractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
        cout << "Induced processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
        if (max_cfractions_vectors.at(cell_type).size() > 0)
            special_c_sequence = PCC_Kinematic_cracking(cell_type, special_x_sequence, Configuration_cState, max_cfractions_vectors);
        // cell type = 2 -> faces
    } // End of 'Km' type simulations (elseif)

    else if (ptype_vector.at(cell_type + (3 - 3)) == "Kn" && max_cfractions_vectors[cell_type + (3 - 3)].size() > 0) { // Maximum <functional> production
        cout << "Induced processing in operation: cell_type : "s << cell_type << endl;
        Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type + (3 - 3) << endl;
//        if (max_cfractions_vectors.at(cell_type).size() > 0)
        //        special_c_sequence = PCC_Kinetic_cracking(Configuration_State, face_elastic_energies, large_crack);
    } // End of 'Km' type simulations (elseif)
    else cout << "ERROR [HAGBsProbability3D] : unknown induced simulation type - please replace with 'Km' or 'Kn'..!" << endl;

//REPAIR    cout << "ctype_vector " << ctype_vector.at(cell_type + (3 - 3)) << "  " << max_cfractions_vectors[cell_type + (3 - 3)].size() << endl;

/// Assigned sequences:
    CD.Set_sequence(special_x_sequence, cell_type + (3 - 3)); // (sequence, id)
    CD.Set_design(Configuration_State.at(cell_type + (3 - 3)), cell_type + (3 - 3)); // (design, id) - design vector with types
/// Induced sequences:
    CD.Set_induced_sequence(special_c_sequence, cell_type + (3 - 3)); // (sequence, ctype)

    } // for (int cell_type = 3; cell_type >= 0; --cell_type)

    cout << endl; Out_logfile_stream << endl;
    cout << "p-sequence size: " << CD.Get_p_sequence().size() << endl; Out_logfile_stream << "p-sequence size: " << CD.Get_p_sequence().size() << endl;
    cout << "f-sequence size: " << CD.Get_f_sequence().size() << endl; Out_logfile_stream << "f-sequence size: " << CD.Get_f_sequence().size() << endl;
    cout << "e-sequence size: " << CD.Get_e_sequence().size() << endl; Out_logfile_stream << "e-sequence size: " << CD.Get_e_sequence().size() << endl;
    cout << "n-sequence size: " << CD.Get_n_sequence().size() << endl; Out_logfile_stream << "n-sequence size: " << CD.Get_n_sequence().size() << endl;
    cout << endl; Out_logfile_stream << endl;

    cout << "p-induced-sequence size: " << CD.Get_p_induced_sequence().size() << endl; Out_logfile_stream << "p-sequence size: " << CD.Get_p_induced_sequence().size() << endl;
    cout << "f-induced-sequence size: " << CD.Get_f_induced_sequence().size() << endl; Out_logfile_stream << "f-sequence size: " << CD.Get_f_induced_sequence().size() << endl;
    cout << "e-induced-sequence size: " << CD.Get_e_induced_sequence().size() << endl; Out_logfile_stream << "e-sequence size: " << CD.Get_e_induced_sequence().size() << endl;
    cout << "n-induced-sequence size: " << CD.Get_n_induced_sequence().size() << endl; Out_logfile_stream << "n-sequence size: " << CD.Get_n_induced_sequence().size() << endl;
    cout << endl; Out_logfile_stream << endl;

    cout << "p-design vector size: " << CD.Get_p_design().size() << endl; Out_logfile_stream << "p-design vector size: " << CD.Get_p_design().size() << endl;
    cout << "f-design vector size: " << CD.Get_f_design().size() << endl; Out_logfile_stream << "f-design vector size: " << CD.Get_f_design().size() << endl;
    cout << "e-design vector size: " << CD.Get_e_design().size() << endl;  Out_logfile_stream << "e-design vector size: " << CD.Get_e_design().size() << endl;
    cout << "n-design vector size: " << CD.Get_n_design().size() << endl; Out_logfile_stream << "n-design vector size: " << CD.Get_n_design().size() << endl;
    cout << endl; Out_logfile_stream << endl;

    return CD;
} /// The end of Processing()

/*
 * /// Off-streams :
ofstream OutFLfile, OutJFile, OutJ2File, OutSFile, OutS2File, OutPowersADistributions, OutAvLengthsADistributions, OutAgglStatistics;

 *
 * std::vector<vector<unsigned int>> DCC_Processing(std::vector<unsigned int>  &special_faces_sequence, vector<unsigned int> &State_sVector, const string P_type, double &mu_f, double &sigm_f, vector<vector<int>> &RW_series_vector) {
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
*/


/// * Heap * ///
// Related functions

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


/*
/// For L mode only !
// PCC special structures
    std::vector<unsigned int> face_strip_distribution; // a vector containing length distribution of special faces chains (strips)
//    double mu_f = 0, sigm_f = 0; // only for L type simulations (!)
    vector<vector<int>> RW_series_vector; // only for L type simulations (!)

 * bool ProcessingON(char* config, bool time_step_one) {
    std::string line;
    ifstream inConf(config);
    bool isProcessingON = 0;

    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (!line.compare("PCC_Processing SWITCHED ON"s)) {
                isProcessingON = 1;
                if (time_step_one == 1) cout << "ON    | PCC_Processing"s << endl;
                return isProcessingON;
            }
    } else cout << "ProcessingON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF    | PCC_Processing"s << endl;
    return isProcessingON;
}
*/

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