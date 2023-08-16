///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifiers for k-Cells, k={0,1,2,3}     **/
/**  in the PCC Processing module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                        **/
///==============================================================================================================================///

//#include "Processing_Assignment_functions.h"

///------------------------------
/// Two basic function first (choosing of a/the particular face(s))
///------------------------------
/// Standard C++ libraries (STL):
#include <random> // Require C++ 11 and above
// #include <execution> // Require C++ 17 and above

/// (1) Quasi-random choice of the element with # New2CellNumb from the list of numbers {0 to OCellsNumb}
/*!
 * @param OCellsNumb // the number of ordinary faces
 * @return New2CellNumb
 */
unsigned int NewCellNumb_R(unsigned int OCellsNumb){ // Random generation machine for a new 2-Cell number
    unsigned int New2CellNumb;
    uniform_int_distribution<size_t> uni_rand (0, OCellsNumb - 1); // uniformly distributed from 0 to OCellsNumb-1 inclusive
    std::random_device rd; // generates unsigned random integers
    std::mt19937 mt(rd());

    New2CellNumb = uni_rand(mt); // Random generation of the boundary number in the range from 0 to OrdinaryCellNumbs.size()-1
    return New2CellNumb;
// Fast standard generator:  return New2CellNumb = ::rand() % OCellsNumb;
}

/// (2) Quasi-random choice of the element with # New2CellNumb from the list of numbers {0 to OCellsNumb}
/*!
 * @param <int> IniFaceNumber // initial number of face
 * @param <int> strip_length // the number of faces in one strip
 * @param <int> Leap_probability = 0 and <double> Leap_dist = 1 // for the future possibility of leaps: probability (hence frequency) and distance ((?)fraction of the complex size in grains)
 * @return vector<int> NewFacesStrip_RW // vector of faces joining in one strip
 */
// strip_sface_distribution - is a vector containin the numbers of strips of each length starting with 1 expressed in # of faces
// Example: vector<int> strip_sface_distribution = {2 4 5 27 8 6 3 1} means 2 strips of length 1 faces each, 4 strips of length 2 faces each,... , 1 strip of length 8 faces each
// int strip_length - is the number of faces in the current strip or RW path length in # of faces
/// This RW choose ANY faces, not necessary only ordinary ones (!)
vector<int> NewFacesStrip_RW( int iniFaceNumber, int strip_length, int Leap_friquency = 1, double Leap_dist = 1) { // Random generation machine for a strips of new 2-Cells
    unsigned int New2CellNumb;  vector<int> NewStripVector_RW;
    vector<double> neigh_Faces; // vector for all neighbours of each face
    /// Sparse Face Adjacency matrix - reading from the file of the considered PCC
    SpMat AFS = SMatrixReader(paths.at(2 + (dim0 - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose()));     //  Full symmetric AFS matrix instead of triagonal
    // Find ordinary-ONLY faces: Calculation OrdinaryCellNumbs vector based on a given S_Vector std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(2), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++) OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces // S_Vector with its non-zero elements set any pre-define structure of special element feeding to the function Processing_Random for( unsigned int itr : S_Vector) if(itr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + itr); // !!! Delete its element from the vector decreasing its size BUT

    /// (1) One random face first as the first face of the strip
    NewStripVector_RW.push_back(iniFaceNumber); // add new element to the NewFacesStrip_RW face vector

    /// (2) Loop over strip_length (in the current basket of the strip lengths distribution)
//    #pragma omp parallel for // parallel execution by OpenMP
    for (int strip_length_counter = 0; strip_length_counter < strip_length; strip_length_counter++) {
        // Looking for all the face neighbours
 //       #pragma omp parallel for // parallel execution by OpenMP
        for (int k = 0; k < CellNumbs.at(2); ++k) // Loop over all the Faces in the PCC
            if (AFS.coeff(New2CellNumb, k) == 1) neigh_Faces.push_back(k); // set of all the face neighbours
        // new random choice between all the face neighbours
        New2CellNumb = neigh_Faces.at(NewCellNumb_R(neigh_Faces.size()));
        NewStripVector_RW.push_back(New2CellNumb); // add new element to the NewFacesStrip_RW face vector
        neigh_Faces.clear(); // clear the vector for the next face neighbours
    } // end of  for (int strip_length_counter = 0; strip_length_counter < f_length; strip_length_counter++) {

    return NewStripVector_RW;
}

///--------------------------------------------------------------------------------
/// Several more complex generation function (generation of face sequences)
///--------------------------------------------------------------------------------
/// (3) The Random generation process function
/*!
 * @param S_Vector
 * @param s_faces_sequence
 * @param max_sFaces_fraction
 */
/// S_Vector with its non-zero elements set any pre-define structure of special element feeding to the function Processing_Random
std::vector<unsigned int> Processing_Random(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> max_fractions_vectors) {

///================================================================= 'R' =======================================================================////
/// ====================================================>  Random generation process  <========================================================////
///===========================================================================================================================================////
std::vector<unsigned int> special_cells_sequence; // output of the function

///=============================================================================================================================================////
/// =====> Initial initialisation with the previous calculation step (if any) based on the "special_faces_sequence" file

// OrdinaryCellNumbs is just a tricky way to acceleration the random process
    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type + (dim0 - 3)), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    vector<int> S_Vector(CellNumbs.at(cell_type + (dim0 - 3)), 0); // S_Vector - State Vector for a given type of k-cells
    if (Configuration_State.at(cell_type + (dim0 - 3)).size() > 0)
        S_Vector = Configuration_State.at(cell_type + (dim0 - 3)); // initial predefined system, if exists
    /// S_Vector with its non-zero elements set any pre-define structure of special element feeding to the function Processing_Random

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr)
        if(*istr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + distance(S_Vector.begin(),istr)); // !!! Delete its element from the vector decreasing its size BUT

    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++j)
        if(max_fractions_vectors[cell_type + (dim0 - 3)][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type + (dim0 - 3)][j];

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3))); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
        cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
        Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

// number of cell types
    int number_of_types = std::count_if(max_fractions_vectors[cell_type + (dim0 - 3)].begin(), max_fractions_vectors[cell_type + (dim0 - 3)].end(), [](int c){return c > 0;});

/// ================= Loop over all ordinary Faces before sFaces_fraction = @parameter max_sFaces_fraction fractions of special cells  =======================>
//    #pragma omp parallel do // parallel execution by OpenMP
    do { // do{ ... }while(output_step) loop starting point
        int NewFaceType = 1;
        unsigned int NewCellNumb = 0;
        NewCellNumb = NewCellNumb_R(OrdinaryCellNumbs.size()); // advanced random generator of  pecial faces
// REPAIR   cout << "in Random function! "s << OrdinaryCellNumbs.size() <<" ncn: " << OrdinaryCellNumbs.at(NewCellNumb) << endl;

// special_cells_fraction
        /// Random generation of types !!! with IDs < number_of_types
        if (number_of_types > 1) // int NewFaceType = 1;
            NewFaceType = rand() % number_of_types; // Random chose for the chose of cell's type to be assigned over all special types

        /// Changes in State vectors
        S_Vector.at(OrdinaryCellNumbs.at(NewCellNumb)) = NewFaceType; // change element of the State Vector according to the NewFaceType
        special_cells_sequence.push_back(OrdinaryCellNumbs.at(NewCellNumb)); // add new element to the s_faces_sequence

        // It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
        OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + NewCellNumb); /// !!! Delete its element from the vector decreasing its size

// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3)); // type (i+1) of special x_cells

 //REPAIR cout << "special_faces_fraction: \t" << special_faces_fraction << "\t\t" << endl;
        /// output calculation progress
        if ((int) (special_cells_fraction*CellNumbs.at(cell_type + (dim0 - 3))) % (int) 0.1 * CellNumbs.at(cell_type + (dim0 - 3)) == 0) cout << "special " << cell_type + (dim0 - 3) << "-cells fraction:      " << special_cells_fraction << endl;
        if ((int) (special_cells_fraction*CellNumbs.at(cell_type + (dim0 - 3))) % (int) 0.1 * CellNumbs.at(cell_type + (dim0 - 3)) == 0) Out_logfile_stream << "special " << cell_type + (dim0 - 3) << "-cells fraction:      " << special_cells_fraction << endl;

        /// test output (!)
        std::vector<double> j_edge_fractions(dim0 + 1, 0), d_edge_fractions(dim0, 0);
//        vector<int> e_state = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_edge_fractions, d_edge_fractions);
//    for (auto fr : e_state)
//        std::cout << fr << " ";
//        cout << endl;

 //       std::cout << "j sum : " << j_edge_fractions.at(1) + j_edge_fractions.at(2) + j_edge_fractions.at(3) << endl;
 //       std::cout << "j1 : " << j_edge_fractions.at(1) << "  j2 : " << j_edge_fractions.at(2) << "  j3 : " << j_edge_fractions.at(3) << endl;
 //       std::cout << "p : " << special_cells_fraction << "  ConFEnt : " << std::get<0>(Configuration_Entropy_tuple(j_edge_fractions)) + std::get<1>(Configuration_Entropy_tuple(j_edge_fractions)) << endl;

    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type].push_back(var);

    return special_cells_sequence;
} // end of Random generation mode

/*
/// (2) Lengthy special strips with the distribution of lengths taken from file
vector<vector<int>> RStrips_Distribution( std::vector<unsigned int> const &face_strip_distribution, std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence) {
///================================================================= 'L' =======================================================================////
/// ==============================================>  Random lengthy strips generation process  <===============================================////
///===========================================================================================================================================////
    int NewFaceType = 1; // Random generation of types with IDs < number_of_types
    vector<vector<int>> RW_series;

    /// Random Walker (RW) start
    // OrdinaryCellNumbs is just a tricky way to acceleration the random process
    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(2), 1); // Vector of the size equal to the total number of faces in DCC initialised with '0's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for(unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++) OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    /// S_Vector with its non-zero elements set any pre-define structure of special element feeding to the function Processing_Random
    for(unsigned int itr : S_Vector)
        if(itr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + itr); // !!! Delete its element from the vector decreasing its size BUT

    // Initial face fractions
    double ordinary_faces_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(2);
    double special_faces_fraction = 1.0 - ordinary_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) return RW_series; // if, after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied

    /// (1) Loop over the vector of the strips distribution (several "baskets")
    for (auto  itr = face_strip_distribution.begin(); itr != face_strip_distribution.end(); ++itr) {

        /// (2) Inside each basket
        int strip_length = (distance(face_strip_distribution.begin(), itr) + 1.0); // strip lengths, starting with 1
        // Example: vector<int> strip_sface_distribution = {2 4 5 27 8 6 3 1} means 2 strips of length 1 faces each, 4 strips of length 2 faces each,... , 1 strip of length 8 faces each
        for (int number_of_lstrips = 0; number_of_lstrips < (*itr); number_of_lstrips++) { // Number of strips of size *itr > 0

            /// Random choice from ALL 2-Cells the iniFaceNumber - initial face for Random Walker start
            int iniFaceNumber = NewCellNumb_R(CellNumbs.at(2)); // random choice function (!) can choose already special face

            /// Ranfdom Walker giving the sequence of faces vector<int> NewStripVector_RW of length strip_length as a result
            vector<int> NewStripVector_RW = NewFacesStrip_RW(iniFaceNumber, strip_length);

            /// Add a new strip to the vector of sfaces-strips
            RW_series.push_back(NewStripVector_RW);

            /// Changes in vectors from Main function - first element of the strip
            for ( int val : NewStripVector_RW) {
                S_Vector.at(val) += 1; // change element of the State Vector
                s_faces_sequence.push_back(val); // add new element to the s_faces_sequence
            }

        } // end of for (int number_of_lstrips = 0; number_of_lstrips < *itr; number_of_lstrips++) { // Number of strips of size *itr > 0

        /// OrdinaryCellNumbs update
        for(unsigned int itr : S_Vector)
            if(itr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + itr); // !!! Delete its element from the vector decreasing its size BUT

        // Special and Ordinary Faces fraction calculation
        ordinary_faces_fraction = OrdinaryCellNumbs.size() / (double) CellNumbs.at(2);
        special_faces_fraction = 1.0 - ordinary_faces_fraction;
    } // end of    for (auto  itr = face_strip_distribution.begin(); itr != face_strip_distribution.end(); ++itr) {
//REPAIR    for (auto a_vector: S_Vector) cout << a_vector << endl;

    return RW_series;
} // end  of Random lengthy inclusions
*/
/// (3) Maximum Functional based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_maxFunctional(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index) {
///=============================================================================================================================================////
///==================================================================== 'F' ===================================================================////
/// ===============================================>  Maximum functional production process   <================================================////
///============================================================================================================================================////
std::vector<unsigned int> special_cells_sequence; // output of the function
    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type].size(); ++j)
        if(max_fractions_vectors[cell_type][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type][j];
//REPAIR    for (auto k : max_fractions_vectors[cell_type]) cout << k << endl;  for (auto k : CellNumbs) cout << k << endl;

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    vector<int> S_Vector(CellNumbs.at(cell_type), 0); // S_Vector - State Vector for a given type of k-cells

    /// ================> Initial Face seeds - initial state for the MAX Functional production algorithm (!)
    if (std::count(Configuration_State.at(cell_type).begin(), Configuration_State.at(cell_type).end(), 1) / (double) CellNumbs.at(cell_type) > 0.05)
        S_Vector = Configuration_State.at(cell_type); // initial predefined system, if exists
    else {
        ///***function std::vector<unsigned int> Processing_Random(cell_type, &Configuration_State, max_fractions_vectors) from PCC_SupportFunctions.h
        std::vector<vector<double>> seed_fractions_vector = {{0.05}, {0.05}, {0.05}, {0.05}}; // max fractions for the initial seed
        special_cells_sequence = Processing_Random(cell_type, Configuration_State, seed_fractions_vector);
    }
// REPAIR cout << "s_faces_sequence.size(): " << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl;
    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's

    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    for (unsigned int itr : special_cells_sequence)
        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr) {
        if (*istr != 0) OrdinaryCellNumbs.erase(std::find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(),distance(S_Vector.begin(), istr)));
// REPAIR        cout << distance(S_Vector.begin(), istr) << "  " << *istr << "  " << OrdinaryCellNumbs.size() << " HERE!!!" << S_Vector.size() << endl;
    } //    exit(0);

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type);
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type)); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
            cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

// number of cell types
int number_of_types = std::count_if(max_fractions_vectors[cell_type].begin(), max_fractions_vectors[cell_type].end(), [](int c){return c > 0;});

/// Vectors for Edges types and Edges-related configuration entropy
int sub_cell_type = 0;
if (cell_type != 0) sub_cell_type = cell_type - 1;
vector<int> EdgeTypes(CellNumbs.at(sub_cell_type), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values

/// WARNING: 4 and 3 here must be calculated in the future !!!
std::vector<double> j_fractions(4, 0), d_fractions(3, 0); // fractions of (1) edges of different types and (2) edges of different degrees
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);
//REPAIR for(auto kl : EdgesTypes) cout << " " <<kl ; cout << endl; exit(10);

    // Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))); // Edges-Faces (3D) or Nodes-Edges (2D) sparse incidence matrix

/// ================ Loop over all Faces ===================
int index_new = 0;
    do { // do{ ... }while(output_step) loop starting point
        // exit(0);

/// An Entropy Increase List calculation for all the Faces in the given DCC
std::vector<vector<int>> cases_list;
std::map<unsigned int, std::vector<unsigned int>> cases_to_sfaces;

    /// new EdgeTypes
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);

    cases_list.clear(); cases_to_sfaces.clear();
    cases_list = Get_cases_list(S_Vector, EdgeTypes, FES,cases_to_sfaces ,p_index);
    /// WARNING: Only one possible Face type (binary model) (!)

std::vector<double> EntropyIncreaseList; // vector with values of configuration entropy increases at conversion of each Face
j_fractions = j_fractions_vector(EdgeTypes);

// cout << "Current Conf Edges Entropy:   " << get<0>(Configuration_Entropy_tuple(j_fractions)) + get<1>(Configuration_Entropy_tuple(j_fractions)) << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        cout << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
    EntropyIncreaseList.clear();
//        cout << "place 2" << endl;
    for (auto NewEdgeTypes : cases_list) {
        EntropyIncreaseList.push_back(Configuration_Entropy_change(j_fractions, NewEdgeTypes));
//        cout << "EntropyIncreaseList: "  << EntropyIncreaseList.back() << endl;
    }
// REPAIR for (auto EIE :   EntropyIncreaseList)  cout << "EIList " << EIE << endl;  exit(0);
//        cout << "place 3" << endl;
/// Number of the cell giving the maximum increase in configuration entropy at its conversion
double case_chosen = 0;
    // the very first special element with S_max
    case_chosen = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR cout << "s_faces_sequence.size(): " << EntropyIncreaseList.size() << " New2CellNumb:  " << New2CellNumb << endl;
//        cout << "place 4" << endl;
/// Form the max_set of the faces with the same value of EntropyIncreaseList
    std::vector<unsigned int> max_set; // possible set of element with equal values of the entropy increase
    max_set.clear();
    max_set.push_back(case_chosen);
    for (auto  itr = EntropyIncreaseList.begin(); itr != EntropyIncreaseList.end(); ++itr) {
        if (*itr == EntropyIncreaseList.at(case_chosen) && distance(EntropyIncreaseList.begin(), itr) != case_chosen)
            max_set.push_back(distance(EntropyIncreaseList.begin(), itr));
    }
    /// Random choice the number of the element in the max_set with the equal Entropy Increase NewCellNumb_R(max_set.size())
/// and then choose the final element New2CellNumb as max_set[NewCellNumb_R(max_set.size())]
    if(max_set.size() > 1)
        case_chosen = max_set.at(NewCellNumb_R(max_set.size()));
//        cout << "place 5" << endl;

    for (unsigned int itr : cases_to_sfaces [case_chosen]) {
//        cout << " cases_to_sfaces [case_chosen]:  " << itr << endl;
        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector
//        cout << "place 5.1 " << endl;
        OrdinaryCellNumbs.erase(std::remove(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr), OrdinaryCellNumbs.end());

        //     auto o_iterator = find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr);
    //    cout << "place 5.2 " << endl;
     //   if (o_iterator <= OrdinaryCellNumbs.end()) /// ?
     //{
       //     cout << " OrdinaryCellNumbs.at(o_iterator):  " << OrdinaryCellNumbs.at(o_iterator) << endl;
        //    OrdinaryCellNumbs.erase(o_iterator); // !!! Delete its element from the vector decreasing its size
    //    }
//    cout << "place 6" << endl;
        /// Add the new element to s_faces_sequence if it is still not here
        if (std::find(special_cells_sequence.begin(), special_cells_sequence.end(), itr) == special_cells_sequence.end())
            special_cells_sequence.push_back(itr);

        EntropyIncreaseList.at(case_chosen) = 0.0; // zero entropy increase for this element

    } // end for(unsigned int c = 0; c < cases_list.size(); ++c)

    //        cout << "place 7" << endl;
/// Special and Ordinary Faces fraction calculation
// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type);
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type); // type (i+1) of special x_cells
 //       cout << "place 8" << endl;
//REPAIRS
        /// ca
// REPAIR        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "  SV: " << S_Vector.size() << "  ctf: " << cases_to_sfaces.size() << "  ms :  " <<  max_set.size() << "  eel: " << EntropyIncreaseList.size() << "  sss :" << special_cells_sequence.size() << " OCN:  " << OrdinaryCellNumbs.size() << endl;
        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "special " << cell_type << "-cells fraction:  " <<  special_cells_fraction << endl;
//        cout << "place 9" << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
 //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type].push_back(var);

    return special_cells_sequence;
} // end of S_max

/// (3) Maximum Functional based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_minConfEntropy(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index) {
///=============================================================================================================================================////
///==================================================================== 'F' ===================================================================////
/// ===============================================>  Maximum functional production process   <================================================////
///============================================================================================================================================////
    std::vector<unsigned int> special_cells_sequence; // output of the function

    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type].size(); ++j)
        if(max_fractions_vectors[cell_type][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type][j];

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    vector<int> S_Vector(CellNumbs.at(cell_type), 0); // S_Vector - State Vector for a given type of k-cells

    /// ================> Initial Face seeds - initial state for the MAX Functional production algorithm (!)
    if (std::count(Configuration_State.at(cell_type).begin(), Configuration_State.at(cell_type).end(), 1) / (double) CellNumbs.at(cell_type) > 0.05)
        S_Vector = Configuration_State.at(cell_type); // initial predefined system, if exists
    else {
        ///***function std::vector<unsigned int> Processing_Random(cell_type, &Configuration_State, max_fractions_vectors) from PCC_SupportFunctions.h
        std::vector<vector<double>> seed_fractions_vector = {{0.05}, {0.05}, {0.05}, {0.05}}; // max fractions for the initial seed
        special_cells_sequence = Processing_Random(cell_type, Configuration_State, seed_fractions_vector);
    }
// REPAIR cout << "s_faces_sequence.size(): " << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl;

    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    for (unsigned int itr : special_cells_sequence)
        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr) {
        if (*istr != 0) OrdinaryCellNumbs.erase(std::find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(),distance(S_Vector.begin(), istr)));
// REPAIR        cout << distance(S_Vector.begin(), istr) << "  " << *istr << "  " << OrdinaryCellNumbs.size() << " HERE!!!" << S_Vector.size() << endl;
    } //    exit(0);
// Old_wrong_version    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr) if(*istr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + distance(S_Vector.begin(),istr)); // !!! Delete its element from the vector decreasing its size BUT

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type);
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type)); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
            cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

// number of cell types
    int number_of_types = std::count_if(max_fractions_vectors[cell_type].begin(), max_fractions_vectors[cell_type].end(), [](int c){return c > 0;});

/// Vectors for Edges types and Edges-related configuration entropy
    int sub_cell_type = 0;
    if (cell_type != 0) sub_cell_type = cell_type - 1;
    vector<int> EdgeTypes(CellNumbs.at(sub_cell_type), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values

    std::vector<double> j_fractions(4, 0), d_fractions(3, 0); // fractions of (1) edges of different types and (2) edges of different degrees
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);
//REPAIR for(auto kl : EdgesTypes) cout << " " <<kl ; cout << endl; exit(10);

    // Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))); // Edges-Faces sparse incidence matrix

/// ================ Loop over all Faces ===================
    do { // do{ ... }while(output_step) loop starting point

/// An Entropy Increase List calculation for all the Faces in the given DCC
        std::vector<vector<int>> cases_list;
        std::map<unsigned int, std::vector<unsigned int>> cases_to_sfaces;

        /// new EdgeTypes
        EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);

        cases_list.clear(); cases_to_sfaces.clear();
        cases_list = Get_cases_list(S_Vector, EdgeTypes, FES,cases_to_sfaces ,p_index);
        /// WARNING: Only one possible Face type (binary model) (!)

//    cout << "place 1" << endl;

        std::vector<double> DevEntropyIncreaseList; // vector with values of configuration entropy increases at conversion of each Face
        j_fractions = j_fractions_vector(EdgeTypes);

// cout << "Current Conf Edges Entropy:   " << get<0>(Configuration_Entropy_tuple(j_fractions)) + get<1>(Configuration_Entropy_tuple(j_fractions)) << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        cout << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        DevEntropyIncreaseList.clear();
//        cout << "place 2" << endl;
        for (auto NewEdgeTypes : cases_list) {
            DevEntropyIncreaseList.push_back(Configuration_DevEntropy_change(j_fractions, NewEdgeTypes));
//        cout << "EntropyIncreaseList: "  << EntropyIncreaseList.back() << endl;
        }
// REPAIR for (auto EIE :   EntropyIncreaseList)  cout << "EIList " << EIE << endl;  exit(0);
//        cout << "place 3" << endl;
/// Number of the cell giving the maximum increase in configuration entropy at its conversion
        double case_chosen = 0;
        // the very first special element with S_max
        case_chosen = std::max_element(std::begin(DevEntropyIncreaseList), std::end(DevEntropyIncreaseList)) - std::begin(DevEntropyIncreaseList); // gives index of the max element
//REPAIR cout << "s_faces_sequence.size(): " << EntropyIncreaseList.size() << " New2CellNumb:  " << New2CellNumb << endl;
//        cout << "place 4" << endl;
/// Form the max_set of the faces with the same value of EntropyIncreaseList
        std::vector<unsigned int> max_set; // possible set of element with equal values of the entropy increase
        max_set.clear();
        max_set.push_back(case_chosen);
        for (auto  itr = DevEntropyIncreaseList.begin(); itr != DevEntropyIncreaseList.end(); ++itr) {
            if (*itr == DevEntropyIncreaseList.at(case_chosen) && distance(DevEntropyIncreaseList.begin(), itr) != case_chosen)
                max_set.push_back(distance(DevEntropyIncreaseList.begin(), itr));
        }
        /// Random choice the number of the element in the max_set with the equal Entropy Increase NewCellNumb_R(max_set.size())
/// and then choose the final element New2CellNumb as max_set[NewCellNumb_R(max_set.size())]
        if(max_set.size() > 1)
            case_chosen = max_set.at(NewCellNumb_R(max_set.size()));
//        cout << "place 5" << endl;
        /////////////////////////////////////////////////////////
        for (unsigned int itr : cases_to_sfaces [case_chosen]) {
//        cout << " cases_to_sfaces [case_chosen]:  " << itr << endl;
            S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector
//        cout << "place 5.1 " << endl;
            OrdinaryCellNumbs.erase(std::remove(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr), OrdinaryCellNumbs.end());

            //     auto o_iterator = find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr);
            //    cout << "place 5.2 " << endl;
            //   if (o_iterator <= OrdinaryCellNumbs.end()) /// ?
            //{
            //     cout << " OrdinaryCellNumbs.at(o_iterator):  " << OrdinaryCellNumbs.at(o_iterator) << endl;
            //    OrdinaryCellNumbs.erase(o_iterator); // !!! Delete its element from the vector decreasing its size
            //    }
//    cout << "place 6" << endl;
            /// Add the new element to s_faces_sequence if it is still not here
            if (std::find(special_cells_sequence.begin(), special_cells_sequence.end(), itr) == special_cells_sequence.end())
                special_cells_sequence.push_back(itr);

            DevEntropyIncreaseList.at(case_chosen) = 0.0; // zero entropy increase for this element

        } // end for(unsigned int c = 0; c < cases_list.size(); ++c)

        //        cout << "place 7" << endl;
/// Special and Ordinary Faces fraction calculation
// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type);
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type); // type (i+1) of special x_cells
        //       cout << "place 8" << endl;
//REPAIRS
        /// ca
// REPAIR        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "  SV: " << S_Vector.size() << "  ctf: " << cases_to_sfaces.size() << "  ms :  " <<  max_set.size() << "  eel: " << EntropyIncreaseList.size() << "  sss :" << special_cells_sequence.size() << " OCN:  " << OrdinaryCellNumbs.size() << endl;
        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "special " << cell_type << "-cells fraction:  " <<  special_cells_fraction << endl;
//        cout << "place 9" << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type].push_back(var);

    return special_cells_sequence;
} // end of S_min

/// (q) Maximum p crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */

std::vector<unsigned int> Processing_maxP_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index) {
///=============================================================================================================================================////
///==================================================================== 'F' ===================================================================////
/// ===============================================>  Maximum functional production process   <================================================////
///============================================================================================================================================////
    std::vector<unsigned int> special_cells_sequence; // output of the function

    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++j)
        if(max_fractions_vectors[cell_type + (dim0 - 3)][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type + (dim0 - 3)][j];

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    vector<int> S_Vector(CellNumbs.at(cell_type + (dim0 - 3)), 0); // S_Vector - State Vector for a given type of k-cells

    /// ================> Initial Face seeds - initial state for the MAX Functional production algorithm (!)
    if (std::count(Configuration_State.at(cell_type + (dim0 - 3)).begin(), Configuration_State.at(cell_type + (dim0 - 3)).end(), 1) / (double) CellNumbs.at(cell_type + (dim0 - 3)) > 0.05)
        S_Vector = Configuration_State.at(cell_type + (dim0 - 3)); // initial predefined system, if exists
    else {
        ///***function std::vector<unsigned int> Processing_Random(cell_type, &Configuration_State, max_fractions_vectors) from PCC_SupportFunctions.h
        std::vector<vector<double>> seed_fractions_vector = {{0.05}, {0.05}, {0.05}, {0.05}}; // max fractions for the initial seed
        special_cells_sequence = Processing_Random(cell_type, Configuration_State, seed_fractions_vector);
    }
// REPAIR cout << "s_faces_sequence.size(): " << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl;

// number of cell types
    int number_of_types = std::count_if(max_fractions_vectors[cell_type + (dim0 - 3)].begin(), max_fractions_vectors[cell_type + (dim0 - 3)].end(), [](int c){return c > 0;});

/// Vectors for Edges types and Edges-related configuration entropy
    vector<int> EdgeTypes(CellNumbs.at(1 + (dim0 - 3)), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    std::vector<double> j_fractions(dim0 + 1, 0), d_fractions(dim0, 0); // fractions of (1) edges of different types and (2) edges of different degrees
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);
//REPAIR for(auto kl : EdgesTypes) cout << " " <<kl ; cout << endl; exit(10);

    // Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim0 - 3)), CellNumbs.at(1 + (dim0 - 3)), CellNumbs.at(2 + (dim0 - 3))); // Edges-Faces sparse incidence matrix
    SpMat GFS = SMatrixReader(paths.at(6 + (dim0 - 3)), CellNumbs.at(2 + (dim0 - 3)), CellNumbs.at(3 + (dim0 - 3))); // Faces-Grains sparse incidence matrix

    SpMat AGS = SMatrixReader(paths.at(3 + (dim0 - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Faces
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

    std::vector<vector<double>> grain_quaternions(CellNumbs.at(3 + (dim0 - 3)), std::vector<double>(4)),
            grain_triangle_quaternions(CellNumbs.at(3 + (dim0 - 3)), std::vector<double>(3)); // grain orientations 2-quaternions and 3-quternions
    double HAGBs_threshold = 15.0; // treshold for the definition of the "high-angle" disorientations

/// Initial TRIANGLE LATTICE
// q1 \in [0, 1], q2 \in [0, 1], q3 \in [0, 1]
// n_grains - number of points, dq - "lattice parameter"
    int n_lattice_points = 10; // is an arbitrary user-defined parameter here (!)
    double  dq = 1.0 / (double) n_lattice_points;

    std::vector<vector<double>> q_coord_vector; // grain rotation lattice
    for(int qi = 0; qi <= n_lattice_points; ++qi) // 1
        for(int qk = 0; qk <= (n_lattice_points - qi); ++qk) // 2
            q_coord_vector.push_back({qi*dq, qk*dq, 1.0 - qi*dq - qk*dq});
// REPAIR    for(auto gtq : q_coord_vector)
//        cout << "q_coord_vector " << gtq.at(0) << " " << gtq.at(1) << " " << gtq.at(2) << " "<< gtq.at(0) + gtq.at(1) + gtq.at(2) << endl;

/// Assigning of the INITIAL grain orientations for each grain in the PCC
    double q0_coord = 0;
    unsigned int new_grain_triangle_coords = 0;
    grain_triangle_quaternions.clear();
    grain_quaternions.clear();

    for (int i = 0; i < CellNumbs.at(3 + (dim0 - 3)); ++i) {

// random choice of a point in a triangle coordinate space
        new_grain_triangle_coords = NewCellNumb_R(q_coord_vector.size());

        std::vector<double> grain_q_quaternion(3), grain_full_quaternion(4);
/// triangle Q3-quaternions
        grain_q_quaternion.at(0) = q_coord_vector[new_grain_triangle_coords][0];
        grain_q_quaternion.at(1) = q_coord_vector[new_grain_triangle_coords][1];
        grain_q_quaternion.at(2) = q_coord_vector[new_grain_triangle_coords][2];

        grain_triangle_quaternions.push_back({grain_q_quaternion.at(0), grain_q_quaternion.at(1), grain_q_quaternion.at(2)});
// REPAIR        cout << "grain_triangle_quaternions " << grain_triangle_quaternions.back().at(0) << " " << grain_triangle_quaternions.back().at(1) << " " << grain_triangle_quaternions.back().at(2) << " "<< grain_triangle_quaternions.back().at(0) + grain_triangle_quaternions.back().at(1) + grain_triangle_quaternions.back().at(2) << endl;

/// full grain quaternions
        /// angle - an additional random generation
        q0_coord = 0.2*rand() / (RAND_MAX + 1.0); //        q0_coord = std::sqrt(1.0 - pow(grain_full_quaternion.at(0),2) - pow(grain_full_quaternion.at(1),2) - pow(grain_full_quaternion.at(2),2));
        // axis
        grain_full_quaternion.at(0) = std::sqrt(grain_q_quaternion.at(0)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(1) = std::sqrt(grain_q_quaternion.at(1)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(2) = std::sqrt(grain_q_quaternion.at(2)*(1.0 - pow(q0_coord,2)));

// full grain quaternions vector
/// (1) RANDOM ASSIGHMENT (Already Most GBs are HAGBs) (!)
//        grain_quaternions.push_back({q0_coord, grain_full_quaternion.at(0), grain_full_quaternion.at(1), grain_full_quaternion.at(2)});
/// (2) Low HAGBs initial ASSIGHMENT (!!)
//        q0_coord = 0.3;
//        grain_quaternions.push_back({q0_coord, std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.4*(1.0 - pow(q0_coord,2)))});
        grain_quaternions.push_back({q0_coord, std::sqrt(0.1*(1.0 - pow(q0_coord,2))), std::sqrt(0.6*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2)))});
// REPAIR cout << "full grain quaternions " << grain_quaternions.back().at(0) << " " << grain_quaternions.back().at(1) << " " << grain_quaternions.back().at(2) << " " << grain_quaternions.back().at(3) << "  " << pow(grain_quaternions.back().at(0),2) + pow(grain_quaternions.back().at(1),2) + pow(grain_quaternions.back().at(2),2) + pow(grain_quaternions.back().at(3),2)<< endl;
    } // end for (int i = 0; i < CellNumbs.at(3 + (dim - 3)); ++i)

    for(auto gtq : grain_quaternions)
        grain_triangle_quaternions.push_back({pow(gtq.at(1), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(2), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(3), 2) / (1.0 - pow(gtq.at(0), 2))});

// REPAIR    for(auto gtq : grain_triangle_quaternions)
//        cout << gtq.at(0) << "  " << gtq.at(1) << "  " << gtq.at(2) << endl;
    std::vector<unsigned int> gb_set; // g_set, gb_special_set;
    std::map<unsigned int, std::vector<unsigned int>> g_gbs_map; // map of grain boundaries for each grain

    g_gbs_map.clear();
    for (unsigned int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // loop over all Grains in a PCC
        gb_set.clear();
        for(unsigned int f = 0; f < CellNumbs.at(2 + (dim0 - 3)); ++f)// loop over all Edges
            if (GFS.coeff(f, g) != 0)
                gb_set.push_back(f); // for each grain in a PCC
/// new map g_gb element
        g_gbs_map.insert(std::pair<unsigned int, std::vector<unsigned int>> (g, gb_set));
    } // end of for (unsigned int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g)

/// Initial Faces State_Vector
    fill(S_Vector.begin(),S_Vector.end(),0); /// Zeroing S-vector (!)
    for (int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // grains
        for(unsigned int ngb : g_gbs_map[g]) { // all the special GBs related with the rotated grain
            for (int g2 = 0; g2 < CellNumbs.at(3 + (dim0 - 3)); ++g2) { // grains
                if (GFS.coeff(ngb, g2) != 0 && g2 != g) {
                    if (isHAGB(grain_quaternions.at(g), grain_quaternions.at(g2), HAGBs_threshold))
                        S_Vector.at(ngb) = 1;
                } // end if()
            } // end if()
        } // end of for(unsigned int ngb : g_gbs_map[g]) {
    } // end of for (int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g) { // grains

    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type + (dim0 - 3)), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

///    for (unsigned int itr : special_cells_sequence)
///        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr)
        if(*istr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + distance(S_Vector.begin(),istr)); // !!! Delete its element from the vector decreasing its size BUT

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3))); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
            cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

/// ================ Loop over all Faces ===================
    do { // do{ ... }while(output_step) loop starting point

/// An Entropy Increase List calculation for all the Faces in the given DCC
        std::vector<vector<int>> cases_list;
        std::map<unsigned int, std::vector<unsigned int>> cases_to_sfaces; // map from cases to special faces (set of special faces for each case from the list)
        std::map<unsigned int, unsigned int> cases_to_grains; // map from cases to grains
        std::map<unsigned int, std::vector<double>> cases_to_new_quaternions; // map from cases to their new quaternions
        /// new EdgeTypes
        EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);

        cases_list.clear(); cases_to_sfaces.clear();
        cases_list = Get_crystallographic_cases_random_list(grain_quaternions, g_gbs_map, dq, HAGBs_threshold, S_Vector, EdgeTypes, GFS, FES, cases_to_grains, cases_to_sfaces, cases_to_new_quaternions, p_index);
        /// WARNING: Only one possible Face type (binary model) (!)

//    cout << "place 1" << endl;
        std::vector<double> EntropyIncreaseList; // vector with values of configuration entropy increases at conversion of each Face
        j_fractions = j_fractions_vector(EdgeTypes);

// cout << "Current Conf Edges Entropy:   " << get<0>(Configuration_Entropy_tuple(j_fractions)) + get<1>(Configuration_Entropy_tuple(j_fractions)) << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
//        cout << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        EntropyIncreaseList.clear();
//        cout << "place 2" << endl;
/// Configuration_Entropy_change function here (!)
//        for (auto NewEdgeTypes : cases_list) {
        for (auto cltr = cases_list.begin(); cltr != cases_list.end(); ++cltr) {
            EntropyIncreaseList.push_back(std::count(cltr->begin(), cltr->end(), 1)/ (double) CellNumbs.at(2 - (dim0 - 3)));
// REPAIR cout << "EElist  " <<  EntropyIncreaseList.back() << endl;
//            EntropyIncreaseList.push_back(cases_to_sfaces[distance(cases_list.begin(),cltr)].size());
///        cout << "EntropyIncreaseList: "  << EntropyIncreaseList.back() << endl;
        }
// REPAIR for (auto EIE :   EntropyIncreaseList)  cout << "EIList " << EIE << endl;  exit(0);
//        cout << "place 3" << endl;
/// Number of the cell giving the maximum increase in configuration entropy at its conversion
        double case_chosen = 0;
        // the very first special element with S_max
        case_chosen = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR
cout << "s_faces.size(): " << cases_to_sfaces[case_chosen].size() << " case_chosen:  " << case_chosen << endl;

//        cout << "place 4" << endl;
/// Form the max_set of the faces with the same value of EntropyIncreaseList
        std::vector<unsigned int> max_set; // possible set of element with equal values of the entropy increase
        max_set.clear();
        max_set.push_back(case_chosen);
        for (auto  itr = EntropyIncreaseList.begin(); itr != EntropyIncreaseList.end(); ++itr) {
            if (*itr == EntropyIncreaseList.at(case_chosen) && distance(EntropyIncreaseList.begin(), itr) != case_chosen)
                max_set.push_back(distance(EntropyIncreaseList.begin(), itr));
        }
        /// Random choice the number of the element in the max_set with the equal Entropy Increase NewCellNumb_R(max_set.size())
/// and then choose the final element New2CellNumb as max_set[NewCellNumb_R(max_set.size())]
        if(max_set.size() > 1)
            case_chosen = max_set.at(NewCellNumb_R(max_set.size()));
//        cout << "place 5" << endl;

/// New MAPS-related changes (!): grain_quaternions update
//REPAIR cout << " check_0 " << case_chosen << "  " <<  cases_to_sfaces[case_chosen].size() << endl;
//cout << " check_1 " << grain_quaternions.at(cases_to_grains[case_chosen]).at(0) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(1) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(2) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(3) << endl;
//cout << " check_2 " << cases_to_new_quaternions[case_chosen].at(0) << "  " << cases_to_new_quaternions[case_chosen].at(1) <<  "  " << cases_to_new_quaternions[case_chosen].at(2) <<  "  " << cases_to_new_quaternions[case_chosen].at(3) << endl;
        grain_quaternions.at(cases_to_grains[case_chosen]) = cases_to_new_quaternions[case_chosen];
//cout << " check_3 " << grain_quaternions.at(cases_to_grains[case_chosen]).at(0) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(1) <<"  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(2) <<"  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(3) << endl;
//exit(0);

/// new special faces
        S_Vector = cases_list[case_chosen];
        for (unsigned int sface : cases_to_sfaces [case_chosen]) {
//        cout << " cases_to_sfaces [case_chosen]:  " << itr << endl;
//            S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector
//        cout << "place 5.1 " << endl;
            OrdinaryCellNumbs.erase(std::remove(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), sface), OrdinaryCellNumbs.end());

            //     auto o_iterator = find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr);
            //    cout << "place 5.2 " << endl;
            //   if (o_iterator <= OrdinaryCellNumbs.end()) /// ?
            //{
            //     cout << " OrdinaryCellNumbs.at(o_iterator):  " << OrdinaryCellNumbs.at(o_iterator) << endl;
            //    OrdinaryCellNumbs.erase(o_iterator); // !!! Delete its element from the vector decreasing its size
            //    }
//    cout << "place 6" << endl;
            /// Add the new element to s_faces_sequence if it is still not here
            if (std::find(special_cells_sequence.begin(), special_cells_sequence.end(), sface) == special_cells_sequence.end())
                special_cells_sequence.push_back(sface);

            EntropyIncreaseList.at(case_chosen) = 0.0; // zero entropy increase for this element

        } // end for(unsigned int c = 0; c < cases_list.size(); ++c)

//                cout << "place 7" << endl;
/// Special and Ordinary Faces fraction calculation
// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3)); // type (i+1) of special x_cells
//               cout << "place 8" << endl;
//REPAIRS
        /// ca
// REPAIR        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "  SV: " << S_Vector.size() << "  ctf: " << cases_to_sfaces.size() << "  ms :  " <<  max_set.size() << "  eel: " << EntropyIncreaseList.size() << "  sss :" << special_cells_sequence.size() << " OCN:  " << OrdinaryCellNumbs.size() << endl;
        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "special " << cell_type << "-cells fraction:  " <<  special_cells_fraction << endl;

        cout << special_cells_fraction << endl;

        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type + (dim0 - 3)].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type + (dim0 - 3)].push_back(var);

    return special_cells_sequence;
} // END of max_p_crystallographic

/// (q) Maximum p crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_Random_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index) {
///=============================================================================================================================================////
///==================================================================== 'F' ===================================================================////
/// ===============================================>  Maximum functional production process   <================================================////
///============================================================================================================================================////
    std::vector<unsigned int> special_cells_sequence; // output of the function

    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type].size(); ++j)
        if(max_fractions_vectors[cell_type][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type][j];

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    vector<int> S_Vector(CellNumbs.at(cell_type), 0); // S_Vector - State Vector for a given type of k-cells

    /// ================> Initial Face seeds - initial state for the MAX Functional production algorithm (!)
    if (std::count(Configuration_State.at(cell_type).begin(), Configuration_State.at(cell_type + (dim0 - 3)).end(), 1) / (double) CellNumbs.at(cell_type + (dim0 - 3)) > 0.05)
        S_Vector = Configuration_State.at(cell_type); // initial predefined system, if exists
    else {
        ///***function std::vector<unsigned int> Processing_Random(cell_type, &Configuration_State, max_fractions_vectors) from PCC_SupportFunctions.h
        std::vector<vector<double>> seed_fractions_vector = {{0.05}, {0.05}, {0.05}, {0.05}}; // max fractions for the initial seed
        special_cells_sequence = Processing_Random(cell_type, Configuration_State, seed_fractions_vector);
    }
// REPAIR cout << "s_faces_sequence.size(): " << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl;

// number of cell types
    int number_of_types = std::count_if(max_fractions_vectors[cell_type].begin(), max_fractions_vectors[cell_type].end(), [](int c){return c > 0;});

/// Vectors for Edges types and Edges-related configuration entropy
    vector<int> EdgeTypes(CellNumbs.at(1 ), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    std::vector<double> j_fractions(4, 0), d_fractions(3, 0); // fractions of (1) edges of different types and (2) edges of different degrees
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);
//REPAIR for(auto kl : EdgesTypes) cout << " " <<kl ; cout << endl; exit(10);

    // Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))); // Edges-Faces sparse incidence matrix
    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), CellNumbs.at(2 + (dim - 3)), CellNumbs.at(3 + (dim - 3))); // Faces-Grains sparse incidence matrix

    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3 + (dim - 3))), (CellNumbs.at(3 + (dim - 3)))); //all Faces
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

    std::vector<vector<double>> grain_quaternions(CellNumbs.at(3 + (dim - 3)), std::vector<double>(4)),
            grain_triangle_quaternions(CellNumbs.at(3 + (dim - 3)), std::vector<double>(3)); // grain orientations 2-quaternions and 3-quternions
    double HAGBs_threshold = 15.0; // treshold for the definition of the "high-angle" disorientations

/// Initial TRIANGLE LATTICE
// q1 \in [0, 1], q2 \in [0, 1], q3 \in [0, 1]
// Number of grains which will be taken from the "EntropyIncreaseList" at each calculation step (!)
    unsigned int grain_set_size = std::floor(0.05*CellNumbs.at(3 + (dim - 3)));

// n_grains - number of points, dq - "lattice parameter"
    int n_lattice_points = 30; // is an arbitrary user-defined parameter here (!)
    double  dq = 1.0 / (double) n_lattice_points;

    std::vector<vector<double>> q_coord_vector; // grain rotation lattice
    for(int qi = 0; qi <= n_lattice_points; ++qi) // 1
        for(int qk = 0; qk <= (n_lattice_points - qi); ++qk) // 2
            q_coord_vector.push_back({qi*dq, qk*dq, 1.0 - qi*dq - qk*dq});
// REPAIR    for(auto gtq : q_coord_vector)
//        cout << "q_coord_vector " << gtq.at(0) << " " << gtq.at(1) << " " << gtq.at(2) << " "<< gtq.at(0) + gtq.at(1) + gtq.at(2) << endl;

/// Assigning of the INITIAL grain orientations for each grain in the PCC
    double q0_coord = 0;
    unsigned int new_grain_triangle_coords = 0;
    grain_triangle_quaternions.clear();
    grain_quaternions.clear();

    for (int i = 0; i < CellNumbs.at(3 + (dim0 - 3)); ++i) {

// random choice of a point in a triangle coordinate space
        new_grain_triangle_coords = NewCellNumb_R(q_coord_vector.size());

        std::vector<double> grain_q_quaternion(3), grain_full_quaternion(4);
/// triangle Q3-quaternions
        grain_q_quaternion.at(0) = q_coord_vector[new_grain_triangle_coords][0];
        grain_q_quaternion.at(1) = q_coord_vector[new_grain_triangle_coords][1];
        grain_q_quaternion.at(2) = q_coord_vector[new_grain_triangle_coords][2];

        grain_triangle_quaternions.push_back({grain_q_quaternion.at(0), grain_q_quaternion.at(1), grain_q_quaternion.at(2)});
// REPAIR        cout << "grain_triangle_quaternions " << grain_triangle_quaternions.back().at(0) << " " << grain_triangle_quaternions.back().at(1) << " " << grain_triangle_quaternions.back().at(2) << " "<< grain_triangle_quaternions.back().at(0) + grain_triangle_quaternions.back().at(1) + grain_triangle_quaternions.back().at(2) << endl;

/// full grain quaternions
        /// angle - an additional random generation
        q0_coord = 0.2*rand() / (RAND_MAX + 1.0); //        q0_coord = std::sqrt(1.0 - pow(grain_full_quaternion.at(0),2) - pow(grain_full_quaternion.at(1),2) - pow(grain_full_quaternion.at(2),2));
        // axis
        grain_full_quaternion.at(0) = std::sqrt(grain_q_quaternion.at(0)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(1) = std::sqrt(grain_q_quaternion.at(1)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(2) = std::sqrt(grain_q_quaternion.at(2)*(1.0 - pow(q0_coord,2)));

// full grain quaternions vector
/// (1) RANDOM ASSIGHMENT (Already Most GBs are HAGBs) (!)
//        grain_quaternions.push_back({q0_coord, grain_full_quaternion.at(0), grain_full_quaternion.at(1), grain_full_quaternion.at(2)});
/// (2) Low HAGBs initial ASSIGHMENT (!!)
//        q0_coord = 0.3;
//        grain_quaternions.push_back({q0_coord, std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.4*(1.0 - pow(q0_coord,2)))});
        grain_quaternions.push_back({q0_coord, std::sqrt(0.1*(1.0 - pow(q0_coord,2))), std::sqrt(0.6*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2)))});
// REPAIR cout << "full grain quaternions " << grain_quaternions.back().at(0) << " " << grain_quaternions.back().at(1) << " " << grain_quaternions.back().at(2) << " " << grain_quaternions.back().at(3) << "  " << pow(grain_quaternions.back().at(0),2) + pow(grain_quaternions.back().at(1),2) + pow(grain_quaternions.back().at(2),2) + pow(grain_quaternions.back().at(3),2)<< endl;
    } // end for (int i = 0; i < CellNumbs.at(3 + (dim - 3)); ++i)

    for(auto gtq : grain_quaternions)
        grain_triangle_quaternions.push_back({pow(gtq.at(1), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(2), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(3), 2) / (1.0 - pow(gtq.at(0), 2))});

// REPAIR for(auto gtq : grain_triangle_quaternions)
//        cout << gtq.at(0) << "  " << gtq.at(1) << "  " << gtq.at(2) << endl; exit(0);
    std::vector<unsigned int> gb_set; // g_set, gb_special_set;
    std::map<unsigned int, std::vector<unsigned int>> g_gbs_map; // map of grain boundaries for each grain

    g_gbs_map.clear();
    for (unsigned int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // loop over all Grains in a PCC
        gb_set.clear();
        for(unsigned int f = 0; f < CellNumbs.at(2 + (dim0 - 3)); ++f)// loop over all Edges
            if (GFS.coeff(f, g) != 0)
                gb_set.push_back(f); // for each grain in a PCC
/// new map g_gb element
        g_gbs_map.insert(std::pair<unsigned int, std::vector<unsigned int>> (g, gb_set));
    } // end of for (unsigned int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g)

/// Initial Faces State_Vector
    fill(S_Vector.begin(),S_Vector.end(),0); /// Zeroing S-vector (!)
    for (int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // grains
        for(unsigned int ngb : g_gbs_map[g]) { // all the special GBs related with the rotated grain
            for (int g2 = 0; g2 < CellNumbs.at(3 + (dim0 - 3)); ++g2) { // grains
                if (GFS.coeff(ngb, g2) != 0 && g2 != g) {
                    if (isHAGB(grain_quaternions.at(g), grain_quaternions.at(g2), HAGBs_threshold))
                        S_Vector.at(ngb) = 1;
                } // end if()
            } // end if()
        } // end of for(unsigned int ngb : g_gbs_map[g]) {
    } // end of for (int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g) { // grains

    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type + (dim0 - 3)), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

///    for (unsigned int itr : special_cells_sequence)
///        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr)
        if (*istr != 0) OrdinaryCellNumbs.erase(std::find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(),distance(S_Vector.begin(), istr)));

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type);
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3))); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
            cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

//    cout << "HerE!" << endl; exit(11);

/// ================ Loop over all Faces ===================
    do { // do{ ... }while(output_step) loop starting point

/// An Entropy Increase List calculation for all the Faces in the given DCC
        std::vector<vector<int>> cases_list;
        std::map<unsigned int, std::vector<unsigned int>> cases_to_sfaces; // map from cases to special faces (set of special faces for each case from the list)
        std::map<unsigned int, unsigned int> cases_to_grains; // map from cases to grains
        std::map<unsigned int, std::vector<double>> cases_to_new_quaternions; // map from cases to their new quaternions
        /// new EdgeTypes
        EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);

        cases_list.clear(); cases_to_sfaces.clear();
        cases_list = Get_crystallographic_cases_random_list(grain_quaternions, g_gbs_map, dq, HAGBs_threshold, S_Vector, EdgeTypes, GFS, FES, cases_to_grains, cases_to_sfaces, cases_to_new_quaternions, p_index);
        /// WARNING: Only one possible Face type (binary model) (!)

//    cout << "place 1  " << grain_quaternions.size() << endl;

        std::vector<double> EntropyIncreaseList; // vector with values of configuration entropy increases at conversion of each Face
        j_fractions = j_fractions_vector(EdgeTypes);

// cout << "Current Conf Edges Entropy:   " << get<0>(Configuration_Entropy_tuple(j_fractions)) + get<1>(Configuration_Entropy_tuple(j_fractions)) << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
//        cout << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        EntropyIncreaseList.clear();
//        cout << "place 2" << endl;
/// Configuration_Entropy_change function here (!)
//        for (auto NewEdgeTypes : cases_list) {
        for (auto cltr = cases_list.begin(); cltr != cases_list.end(); ++cltr) {
            EntropyIncreaseList.push_back(std::count(cltr->begin(), cltr->end(), 1)/ (double) CellNumbs.at(cell_type));
// REPAIR cout << "EElist  " <<  EntropyIncreaseList.back() << endl;
//            EntropyIncreaseList.push_back(cases_to_sfaces[distance(cases_list.begin(),cltr)].size());
///        cout << "EntropyIncreaseList: "  << EntropyIncreaseList.back() << endl;
        }
// REPAIR for (auto EIE :   EntropyIncreaseList)  cout << "EIList " << EIE << endl;  exit(0);
//        cout << "place 3" << endl;
/// Number of the cell giving the maximum increase in configuration entropy at its conversion
//        double case_chosen = 0;
        // the very first special element with S_max
//        case_chosen = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR
        //      cout << "s_faces.size(): " << cases_to_sfaces[case_chosen].size() << " case_chosen:  " << case_chosen << endl;

//       cout << "place 4  " << cases_list.size() << endl;
/// sorting the "EntropyIncreaseList"
        std::sort(EntropyIncreaseList.begin(),EntropyIncreaseList.end());

/// Form the max_set of the faces with the max values of EntropyIncreaseList
        std::vector<double> EntropyIncreaseList_temp = EntropyIncreaseList;
        std::vector<unsigned int> case_chosen_set; // possible set of element with equal values of the entropy increase
        case_chosen_set.clear();
//        cout << "place 4.1  " << EntropyIncreaseList_temp.size() << endl;

//        for(auto g_iter = EntropyIncreaseList.begin(); g_iter < grain_set_size; ++g_iter) {
        for(unsigned int g_it = 0; g_it < grain_set_size; ++g_it) {
            case_chosen_set.push_back(NewCellNumb_R(EntropyIncreaseList_temp.size()));
            EntropyIncreaseList_temp.erase(EntropyIncreaseList_temp.begin() + case_chosen_set.back());
        }
//        cout << "place 4.2" << endl;

        for (auto ccs : case_chosen_set) {
/// New MAPS-related changes (!): grain_quaternions update
//REPAIR cout << " check_0 " << case_chosen << "  " <<  cases_to_sfaces[case_chosen].size() << endl;
//cout << " check_1 " << grain_quaternions.at(cases_to_grains[case_chosen]).at(0) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(1) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(2) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(3) << endl;
//cout << " check_2 " << cases_to_new_quaternions[case_chosen].at(0) << "  " << cases_to_new_quaternions[case_chosen].at(1) <<  "  " << cases_to_new_quaternions[case_chosen].at(2) <<  "  " << cases_to_new_quaternions[case_chosen].at(3) << endl;
            grain_quaternions.at(cases_to_grains[ccs]) = cases_to_new_quaternions[ccs];
//cout << " check_3 " << grain_quaternions.at(cases_to_grains[case_chosen]).at(0) << "  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(1) <<"  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(2) <<"  " << grain_quaternions.at(cases_to_grains[case_chosen]).at(3) << endl;
//exit(0);

/// new special faces
            S_Vector = cases_list[ccs];
            for (unsigned int sface : cases_to_sfaces [ccs]) {
//        cout << " cases_to_sfaces [case_chosen]:  " << itr << endl;
//            S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector
//        cout << "place 5.1 " << endl;
                OrdinaryCellNumbs.erase(std::remove(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), sface),
                                        OrdinaryCellNumbs.end());

                //     auto o_iterator = find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr);
                //    cout << "place 5.2 " << endl;
                //   if (o_iterator <= OrdinaryCellNumbs.end()) /// ?
                //{
                //     cout << " OrdinaryCellNumbs.at(o_iterator):  " << OrdinaryCellNumbs.at(o_iterator) << endl;
                //    OrdinaryCellNumbs.erase(o_iterator); // !!! Delete its element from the vector decreasing its size
                //    }
//    cout << "place 6" << endl;
                /// Add the new element to s_faces_sequence if it is still not here
                if (std::find(special_cells_sequence.begin(), special_cells_sequence.end(), sface) ==
                    special_cells_sequence.end())
                    special_cells_sequence.push_back(sface);

                EntropyIncreaseList.at(ccs) = 0.0; // zero entropy increase for this element
            } // end for(unsigned int c = 0; c < cases_list.size(); ++c)

        } // end for (auto ccs : case_chosen_set) {

//                cout << "place 7" << endl;
/// Special and Ordinary Faces fraction calculation
// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3)); // type (i+1) of special x_cells
//               cout << "place 8" << endl;
//REPAIRS
        /// ca
// REPAIR        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "  SV: " << S_Vector.size() << "  ctf: " << cases_to_sfaces.size() << "  ms :  " <<  max_set.size() << "  eel: " << EntropyIncreaseList.size() << "  sss :" << special_cells_sequence.size() << " OCN:  " << OrdinaryCellNumbs.size() << endl;
        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "special " << cell_type << "-cells fraction:  " <<  special_cells_fraction << endl;
        cout << special_cells_fraction << endl;

        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type].push_back(var);

    return special_cells_sequence;
} // END of Random_crystallographic

/// (y) Maximum Functional crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_maxF_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index) {
///=============================================================================================================================================////
///==================================================================== 'F' ===================================================================////
/// ===============================================>  Maximum functional production process   <================================================////
///============================================================================================================================================////
    std::vector<unsigned int> special_cells_sequence; // output of the function

    // calculation of the total max special cell fraction
    double total_max_sCell_fraction = 0;
    for (int j = 0; j < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++j)
        if(max_fractions_vectors[cell_type + (dim0 - 3)][j] > 0)
            total_max_sCell_fraction += max_fractions_vectors[cell_type + (dim0 - 3)][j];

    if (total_max_sCell_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type <<" total_max_sCell_fraction of " << cell_type << "-cells in the processing.ini file = " << total_max_sCell_fraction << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    }
    else if ( total_max_sCell_fraction == 0.0) return special_cells_sequence;

    vector<int> S_Vector(CellNumbs.at(cell_type + (dim0 - 3)), 0); // S_Vector - State Vector for a given type of k-cells

    /// ================> Initial Face seeds - initial state for the MAX Functional production algorithm (!)
    if (std::count(Configuration_State.at(cell_type + (dim0 - 3)).begin(), Configuration_State.at(cell_type + (dim0 - 3)).end(), 1) / (double) CellNumbs.at(cell_type + (dim0 - 3)) > 0.05)
        S_Vector = Configuration_State.at(cell_type + (dim0 - 3)); // initial predefined system, if exists
    else {
        ///***function std::vector<unsigned int> Processing_Random(cell_type, &Configuration_State, max_fractions_vectors) from PCC_SupportFunctions.h
        std::vector<vector<double>> seed_fractions_vector = {{0.05}, {0.05}, {0.05}, {0.05}}; // max fractions for the initial seed
        special_cells_sequence = Processing_Random(cell_type, Configuration_State, seed_fractions_vector);
    }
// REPAIR cout << "s_faces_sequence.size(): " << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl;

// number of cell types
    int number_of_types = std::count_if(max_fractions_vectors[cell_type + (dim0 - 3)].begin(), max_fractions_vectors[cell_type + (dim0 - 3)].end(), [](int c){return c > 0;});

/// Vectors for Edges types and Edges-related configuration entropy
    vector<int> EdgeTypes(CellNumbs.at(1 + (dim0 - 3)), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    std::vector<double> j_fractions(dim0 + 1, 0), d_fractions(dim0, 0); // fractions of (1) edges of different types and (2) edges of different degrees
    EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);
//REPAIR for(auto kl : EdgesTypes) cout << " " <<kl ; cout << endl; exit(10);

    // Obtaining Faces (coloumns) - Edges (rows) Incidence matrix B2 using the file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim0 - 3)), CellNumbs.at(1 + (dim0 - 3)), CellNumbs.at(2 + (dim0 - 3))); // Edges-Faces sparse incidence matrix
    SpMat GFS = SMatrixReader(paths.at(6 + (dim0 - 3)), CellNumbs.at(2 + (dim0 - 3)), CellNumbs.at(3 + (dim0 - 3))); // Faces-Grains sparse incidence matrix

    SpMat AGS = SMatrixReader(paths.at(3 + (dim0 - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Faces
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

    std::vector<vector<double>> grain_quaternions(CellNumbs.at(3 + (dim0 - 3)), std::vector<double>(4)),
                                    grain_triangle_quaternions(CellNumbs.at(3 + (dim0 - 3)), std::vector<double>(3)); // grain orientations 2-quaternions and 3-quternions
    double HAGBs_threshold = 15.0; // treshold for the definition of the "high-angle" disorientations

/// Initial TRIANGLE LATTICE
// q1 \in [0, 1], q2 \in [0, 1], q3 \in [0, 1]
// n_grains - number of points, dq - "lattice parameter"
    int n_lattice_points = 20; // is an arbitrary user-defined parameter here (!)
    double  dq = 1.0 / (double) n_lattice_points;

    std::vector<vector<double>> q_coord_vector; // grain rotation lattice
    for(int qi = 0; qi <= n_lattice_points; ++qi) // 1
        for(int qk = 0; qk <= (n_lattice_points - qi); ++qk) // 2
                q_coord_vector.push_back({qi*dq, qk*dq, 1.0 - qi*dq - qk*dq});
// REPAIR    for(auto gtq : q_coord_vector)
//        cout << "q_coord_vector " << gtq.at(0) << " " << gtq.at(1) << " " << gtq.at(2) << " "<< gtq.at(0) + gtq.at(1) + gtq.at(2) << endl;

/// Assigning of the INITIAL grain orientations for each grain in the PCC
    double q0_coord = 0;
    unsigned int new_grain_triangle_coords = 0;
    grain_triangle_quaternions.clear();
    grain_quaternions.clear();

    for (int i = 0; i < CellNumbs.at(3 + (dim0 - 3)); ++i) {

// random choice of a point in a triangle coordinate space
        new_grain_triangle_coords = NewCellNumb_R(q_coord_vector.size());

        std::vector<double> grain_q_quaternion(3), grain_full_quaternion(4);
/// triangle Q3-quaternions
        grain_q_quaternion.at(0) = q_coord_vector[new_grain_triangle_coords][0];
        grain_q_quaternion.at(1) = q_coord_vector[new_grain_triangle_coords][1];
        grain_q_quaternion.at(2) = q_coord_vector[new_grain_triangle_coords][2];

        grain_triangle_quaternions.push_back({grain_q_quaternion.at(0), grain_q_quaternion.at(1), grain_q_quaternion.at(2)});
// REPAIR        cout << "grain_triangle_quaternions " << grain_triangle_quaternions.back().at(0) << " " << grain_triangle_quaternions.back().at(1) << " " << grain_triangle_quaternions.back().at(2) << " "<< grain_triangle_quaternions.back().at(0) + grain_triangle_quaternions.back().at(1) + grain_triangle_quaternions.back().at(2) << endl;

/// full grain quaternions
    /// angle - an additional random generation
        q0_coord = 0.2*rand() / (RAND_MAX + 1.0); //        q0_coord = std::sqrt(1.0 - pow(grain_full_quaternion.at(0),2) - pow(grain_full_quaternion.at(1),2) - pow(grain_full_quaternion.at(2),2));
        // axis
        grain_full_quaternion.at(0) = std::sqrt(grain_q_quaternion.at(0)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(1) = std::sqrt(grain_q_quaternion.at(1)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(2) = std::sqrt(grain_q_quaternion.at(2)*(1.0 - pow(q0_coord,2)));

// full grain quaternions vector
/// (1) RANDOM ASSIGHMENT (Already Most GBs are HAGBs) (!)
//        grain_quaternions.push_back({q0_coord, grain_full_quaternion.at(0), grain_full_quaternion.at(1), grain_full_quaternion.at(2)});
/// (2) Low HAGBs initial ASSIGHMENT (!!)
//        q0_coord = 0.3;
//        grain_quaternions.push_back({q0_coord, std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2))), std::sqrt(0.4*(1.0 - pow(q0_coord,2)))});
        grain_quaternions.push_back({q0_coord, std::sqrt(0.1*(1.0 - pow(q0_coord,2))), std::sqrt(0.6*(1.0 - pow(q0_coord,2))), std::sqrt(0.3*(1.0 - pow(q0_coord,2)))});
// REPAIR cout << "full grain quaternions " << grain_quaternions.back().at(0) << " " << grain_quaternions.back().at(1) << " " << grain_quaternions.back().at(2) << " " << grain_quaternions.back().at(3) << "  " << pow(grain_quaternions.back().at(0),2) + pow(grain_quaternions.back().at(1),2) + pow(grain_quaternions.back().at(2),2) + pow(grain_quaternions.back().at(3),2)<< endl;
    } // end for (int i = 0; i < CellNumbs.at(3 + (dim - 3)); ++i)

    for(auto gtq : grain_quaternions)
        grain_triangle_quaternions.push_back({pow(gtq.at(1), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(2), 2) / (1.0 - pow(gtq.at(0), 2)), pow(gtq.at(3), 2) / (1.0 - pow(gtq.at(0), 2))});

// REPAIR    for(auto gtq : grain_triangle_quaternions)
//        cout << gtq.at(0) << "  " << gtq.at(1) << "  " << gtq.at(2) << endl;


    /* REPAIR
    /// Mackenzie check
    std::vector<double> disorientations_vector;
    for(unsigned int gr = 0; gr < CellNumbs.at(3 - (dim -3)); ++gr)
        for(unsigned int grn = 0; grn < gr; ++grn)
            if(AGS.coeff(gr,grn) != 0) {
                //        cout << "coeff " << AGS.coeff(gr,grn) << endl;
                disorientations_vector.push_back(Get_2grains_FCCdisorientation(grain_quaternions.at(gr), grain_quaternions.at(grn)));
                cout << disorientations_vector.back()*57.3 << endl;

            }
*/
    std::vector<unsigned int> gb_set; // g_set, gb_special_set;
    std::map<unsigned int, std::vector<unsigned int>> g_gbs_map; // map of grain boundaries for each grain

    g_gbs_map.clear();
    for (unsigned int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // loop over all Grains in a PCC
        gb_set.clear();
        for(unsigned int f = 0; f < CellNumbs.at(2 + (dim0 - 3)); ++f)// loop over all Edges
            if (GFS.coeff(f, g) != 0)
                gb_set.push_back(f); // for each grain in a PCC
/// new map g_gb element
        g_gbs_map.insert(std::pair<unsigned int, std::vector<unsigned int>> (g, gb_set));
    } // end of for (unsigned int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g)

    /// Initial Faces State_Vector
    fill(S_Vector.begin(),S_Vector.end(),0); /// Zeroing S-vector (!)
    for (int g = 0; g < CellNumbs.at(3 + (dim0 - 3)); ++g) { // grains
        for(unsigned int ngb : g_gbs_map[g]) { // all the special GBs related with the rotated grain
            for (int g2 = 0; g2 < CellNumbs.at(3 + (dim0 - 3)); ++g2) { // grains
                if (GFS.coeff(ngb, g2) != 0 && g2 != g) {
                    if (isHAGB(grain_quaternions.at(g), grain_quaternions.at(g2), HAGBs_threshold))
                        S_Vector.at(ngb) = 1;
                } // end if()
            } // end if()
        } // end of for(unsigned int ngb : g_gbs_map[g]) {
    } // end of for (int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g) { // grains

    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(cell_type + (dim0 - 3)), 1); // Vector of the size equal to the total number of faces in PCC initialised with '1's
    // (!) all the cell Numbers start with 0, not 1 like in Neper, Matlab, Fortran and many other software
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++)
        OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

///    for (unsigned int itr : special_cells_sequence)
///        S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    for (auto istr = S_Vector.begin(); istr != S_Vector.end(); ++istr)
        if(*istr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + distance(S_Vector.begin(),istr)); // !!! Delete its element from the vector decreasing its size BUT

    // initial fractions of special cells
    double ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
    double special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector

    vector<double> scell_fractions_vector; // vector for all different types of special cells
    for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i) {
        scell_fractions_vector.push_back(std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3))); // type (i+1) of special x_cells
        if (special_cells_fraction >= total_max_sCell_fraction) {
            cout << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            Out_logfile_stream << "WARNING [Processing module]:" << "initial special cells fraction is already GREATER than the total max special cell fraction from processing.ini file!"s << endl;
            return special_cells_sequence;
        }     // (!) If after the initial set of special faces by their definition in S_Vector their fraction appeared to be larger than max_sFaces_fraction, so the condition for finishing the Processing module are satisfied
    } // end for (int i = 0; i < max_fractions_vectors.size(); ++i)

/// ================ Loop over all Faces ===================
    do { // do{ ... }while(output_step) loop starting point

/// An Entropy Increase List calculation for all the Faces in the given DCC
        std::vector<vector<int>> cases_list;
        std::map<unsigned int, std::vector<unsigned int>> cases_to_sfaces; // map from cases to special faces (set of special faces for each case from the list)
        std::map<unsigned int, unsigned int> cases_to_grains; // map from cases to grains
        std::map<unsigned int, std::vector<double>> cases_to_new_quaternions; // map from cases to their new quaternions
        /// new EdgeTypes
        EdgeTypes = Edge_types_byFaces(CellNumbs, special_cells_sequence, j_fractions, d_fractions);

        cases_list.clear(); cases_to_sfaces.clear();
        cases_list = Get_crystallographic_cases_list(grain_quaternions, g_gbs_map, dq, HAGBs_threshold, S_Vector, EdgeTypes, GFS, FES, cases_to_grains, cases_to_sfaces, cases_to_new_quaternions, p_index);
        /// WARNING: Only one possible Face type (binary model) (!)

//    cout << "place 1" << endl;
        std::vector<double> EntropyIncreaseList; // vector with values of configuration entropy increases at conversion of each Face
        j_fractions = j_fractions_vector(EdgeTypes);

// cout << "Current Conf Edges Entropy:   " << get<0>(Configuration_Entropy_tuple(j_fractions)) + get<1>(Configuration_Entropy_tuple(j_fractions)) << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        cout << "  Configuration_Entropy: " << Configuration_Entropy(EdgeTypes) <<endl;
        EntropyIncreaseList.clear();
//        cout << "place 2" << endl;
/// Configuration_Entropy_change function here (!)
        for (auto NewEdgeTypes : cases_list) {
            EntropyIncreaseList.push_back(std::abs(Configuration_Entropy_change(j_fractions, NewEdgeTypes)));
//        cout << "EntropyIncreaseList: "  << EntropyIncreaseList.back() << endl;
        }
// REPAIR for (auto EIE :   EntropyIncreaseList)  cout << "EIList " << EIE << endl;  exit(0);
//        cout << "place 3" << endl;
/// Number of the cell giving the maximum increase in configuration entropy at its conversion
        double case_chosen = 0;
        // the very first special element with S_max
        case_chosen = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR cout << "s_faces_sequence.size(): " << EntropyIncreaseList.size() << " New2CellNumb:  " << New2CellNumb << endl;
//        cout << "place 4" << endl;
/// Form the max_set of the faces with the same value of EntropyIncreaseList
        std::vector<unsigned int> max_set; // possible set of element with equal values of the entropy increase
        max_set.clear();
        max_set.push_back(case_chosen);
        for (auto  itr = EntropyIncreaseList.begin(); itr != EntropyIncreaseList.end(); ++itr) {
            if (*itr == EntropyIncreaseList.at(case_chosen) && distance(EntropyIncreaseList.begin(), itr) != case_chosen)
                max_set.push_back(distance(EntropyIncreaseList.begin(), itr));
        }
        /// Random choice the number of the element in the max_set with the equal Entropy Increase NewCellNumb_R(max_set.size())
/// and then choose the final element New2CellNumb as max_set[NewCellNumb_R(max_set.size())]
        if(max_set.size() > 1)
            case_chosen = max_set.at(NewCellNumb_R(max_set.size()));
//        cout << "place 5" << endl;

/// New MAPS-related changes (!): grain_quaternions update
        grain_quaternions.at(cases_to_grains[case_chosen]) = cases_to_new_quaternions[case_chosen];
// REPAIR
if (special_cells_fraction >  0.60 && special_cells_fraction <  0.61) {
    cout << "grain_triangle_quaternions" << " at p =  " << special_cells_fraction << endl;
    for (auto gtq: grain_quaternions)
        cout << pow(gtq.at(1), 2) / (1.0 - pow(gtq.at(0), 2)) << " " << pow(gtq.at(2), 2) / (1.0 - pow(gtq.at(0), 2))
             << " " << pow(gtq.at(3), 2) / (1.0 - pow(gtq.at(0), 2)) << endl;
}
// new special faces
        for (unsigned int itr : cases_to_sfaces [case_chosen]) {
//        cout << " cases_to_sfaces [case_chosen]:  " << itr << endl;
            S_Vector.at(itr) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector
//        cout << "place 5.1 " << endl;
            OrdinaryCellNumbs.erase(std::remove(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr), OrdinaryCellNumbs.end());

            //     auto o_iterator = find(OrdinaryCellNumbs.begin(), OrdinaryCellNumbs.end(), itr);
            //    cout << "place 5.2 " << endl;
            //   if (o_iterator <= OrdinaryCellNumbs.end()) /// ?
            //{
            //     cout << " OrdinaryCellNumbs.at(o_iterator):  " << OrdinaryCellNumbs.at(o_iterator) << endl;
            //    OrdinaryCellNumbs.erase(o_iterator); // !!! Delete its element from the vector decreasing its size
            //    }
//    cout << "place 6" << endl;
            /// Add the new element to s_faces_sequence if it is still not here
            if (std::find(special_cells_sequence.begin(), special_cells_sequence.end(), itr) == special_cells_sequence.end())
                special_cells_sequence.push_back(itr);

            EntropyIncreaseList.at(case_chosen) = 0.0; // zero entropy increase for this element

        } // end of for (unsigned int itr : cases_to_sfaces [case_chosen])

        //        cout << "place 7" << endl;
/// Special and Ordinary Faces fraction calculation
// Special and Ordinary Faces fraction recalculation
        ordinary_cells_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(cell_type + (dim0 - 3));
        special_cells_fraction = 1.0 - ordinary_cells_fraction; // special face vecror definition based on the ordinary face vector
        for (int i = 0; i < max_fractions_vectors[cell_type + (dim0 - 3)].size(); ++i)
            scell_fractions_vector.at(i) = std::count(S_Vector.begin(), S_Vector.end(), (i + 1)) / (double) CellNumbs.at(cell_type + (dim0 - 3)); // type (i+1) of special x_cells
        //       cout << "place 8" << endl;
//REPAIRS
        /// ca
// REPAIR        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "  SV: " << S_Vector.size() << "  ctf: " << cases_to_sfaces.size() << "  ms :  " <<  max_set.size() << "  eel: " << EntropyIncreaseList.size() << "  sss :" << special_cells_sequence.size() << " OCN:  " << OrdinaryCellNumbs.size() << endl;
        if ((int) (10.0*special_cells_fraction) % 40 == 0) cout << "special " << cell_type << "-cells fraction:  " <<  special_cells_fraction << endl;

        cout << special_cells_fraction << endl;

        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
        //       if ((int) (special_cells_fraction) % 100000*CellNumbs.at(3) == 0) Out_logfile_stream << "special " << cell_type << "-cells fraction :      " <<  special_cells_fraction << endl;
    } while(special_cells_fraction < total_max_sCell_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Update of the corresponding Configuration State vector
    Configuration_State[cell_type + (dim0 - 3)].clear();
    for (int var : S_Vector )
        Configuration_State[cell_type + (dim0 - 3)].push_back(var);

    return special_cells_sequence;
} // END of S_max_crystallographic

/*
/// (4) i(p) index govern processed
int Processing_ipIndex(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, int index_type, double ip_index) {
///=============================================================================================================================================////
///========================================================================= 'X' (ip INDEX) =================================================================////
/// ====================================================>  i(p) Index-based production process   <========================================================////
///=============================================================================================================================================////
    /// Local functions declarations
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Sstr0 = 0, Sstr00 = 0, SstrN = 0, Sstr0N = 0, Configuration_Face_Entropy = 0, Face_Entropy_Median = 0, Face_Entropy_Skrew = 0;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file paths.at(5 + (dim - 3))
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal

    double special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
    double ordinary_faces_fraction = 1.0 - special_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) {
        cout << "(!) Early exit: Fraction of special faces is already >= MAX values from config.txt" << endl;
        return 0;
    }
    else if (special_faces_fraction < 0.05) {
        /// ================> Initial Face seeds (initial state for the MAX entropy production algorithm)
        Processing_Random(S_Vector, s_faces_sequence, 0.05);
///***function
    }

    /// ========== Calculation configuration entropy at each calculation step ===============>
    /// Vectors for Edges types and Edges-related configuration entropy
    vector<double> TJsTypes(CellNumbs.at(1), 0); // vector int in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
    // Zeroing coefficients of TJsTypes vector
    std::fill(TJsTypes.begin(), TJsTypes.end(), 0);
    /// TJs type calculations
    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
///***function
//REPAIR for(auto kl : TJsTypes) cout << " " <<kl ; cout << endl; exit(10);
    J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    Configuration_Face_Entropy = 0;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions | (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0 / Jall;
    j1 = J1 / Jall;
    j2 = J2 / Jall;
    j3 = J3 / Jall;
    /// Configuration Entropy related with Faces
    Configuration_Face_Entropy = -(j0 * log2(j0) + j1 * log2(j1) + j2 * log2(j2) + j3 * log2(j3));

    /// ======== Loop over all Faces ===========>
    vector<double> EntropyIncreaseList(CellNumbs.at(2), 0); // vector with values of configuration entropy increases at conversion of each Face
    for (unsigned int k = 0; k < CellNumbs.at(2); k++) { // loop over all Faces in DCC
        double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_Entropy = 0, CFace_EntropyIncrease = 0;

        if (S_Vector.at(k) == 0) { // Loop over each still ORDINARY element neighbours
            J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;

            vector<double> j_types_neigh_fractions = GBIndex(k, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
///***function
//REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;
            // Values before conversion
            J00 = j_types_neigh_fractions.at(0);
            J10 = j_types_neigh_fractions.at(1);
            J20 = j_types_neigh_fractions.at(2);
//REPAIR    cout << " J00= " << J00<< " J10= " << J10 << " J20= " << J20 << endl;
            // Values after conversion
            J1N = J00;
            J2N = J10;
            J3N = J20;
            // The entropy increase calculation for a given Face
            CFace_EntropyIncrease = (J0 * log2(J0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                    + (J1 * log2(J1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) +
                                       J1N * log2(J1N + pow(10, -30)))
                                    + (J2 * log2(J2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) +
                                       J2N * log2(J2N + pow(10, -30)))
                                    + (J3 * log2(J3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
//REPAIR  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;

            // The result of one iteration
            EntropyIncreaseList.at(k) = CFace_EntropyIncrease;

        } // if OrdinaryCells (S_Vector.at(Face) == 0)
    } // for (..k < CellNumbs.at(2)..)

    double New2CellNumb = 0; // Only one possible Face type (binary model)
    /// Number of element giving the maximum increase in configuration entropy at its conversion
    New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
///min  New2CellNumb = std::min_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
//REPAIR        cout << s_faces_sequence.size() << "   " << New2CellNumb << endl;

    // Then all the corresponding maps chain
    S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

    /// Add the new element to s_faces_sequence if it is still not here
    if (find(s_faces_sequence.begin(), s_faces_sequence.end(), New2CellNumb) == s_faces_sequence.end())
        s_faces_sequence.push_back(New2CellNumb);

    EntropyIncreaseList.at(New2CellNumb) = 0.0; ///working_max

///=============================================================================================================================================////
/// ================= S_MAX LOOP over the fractions of special cells [0.05,1] =======================>
    vector<double> new_neigh_TJs, new_neigh_Faces;
    do { // do{ ... }while(output_step) loop starting point

        /// Neighbours of the converted Face
        new_neigh_TJs.clear();
        new_neigh_Faces.clear();
        for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the edges
            if (FES.coeff(k, New2CellNumb) == 1) new_neigh_TJs.push_back(k);
        if (new_neigh_TJs.size() > 0) for (auto itr: new_neigh_TJs) TJsTypes.at(itr)++;

        for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
            if (AFS.coeff(New2CellNumb, m) == 1 && S_Vector.at(m) == 0) new_neigh_Faces.push_back(m);

        if (new_neigh_Faces.size() > 0) {
            for (auto itr: new_neigh_Faces) {
                vector<double> j_types_neigh_fractions = GBIndex(itr, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                /// New EntropyIncreaseList changes
                double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
                // Values before conversion
                J00 = j_types_neigh_fractions.at(0);
                J10 = j_types_neigh_fractions.at(1);
                J20 = j_types_neigh_fractions.at(2);
                // Values after conversion
                J1N = J00;
                J2N = J10;
                J3N = J20;
                // The entropy increase calculation for a given Face
                CFace_EntropyIncrease = (J0 * log2(J0 + pow(10, -30)) - J00 * log2(J00 + pow(10, -30)))
                                        + (J1 * log2(J1 + pow(10, -30)) - J10 * log2(J10 + pow(10, -30)) + J1N * log2(J1N + pow(10, -30)))
                                        + (J2 * log2(J2 + pow(10, -30)) - J20 * log2(J20 + pow(10, -30)) + J2N * log2(J2N + pow(10, -30)))
                                        + (J3 * log2(J3 + pow(10, -30)) + J3N * log2(J3N + pow(10,-30)));
                /// The result of the one iteration
                EntropyIncreaseList.at(itr) = CFace_EntropyIncrease; // for(auto kl :EntropyIncreaseList) cout << kl << endl;

            } // for (auto itr : new_neigh_Faces )
        } //if

        double New2CellNumb = 0; // Only one possible Face type (binary model)
        double min_val = pow(10, 30);
        /// Number of element giving the maximum increase in configuration entropy at its conversion
        unsigned int ipl = 0;
        for ( auto pl : EntropyIncreaseList ) {
            if (pl < min_val) {
                min_val = pl;
                New2CellNumb = (double) ipl;
            }
            ++ipl;
        }
//REPAIR     cout  << "\t\t" <<  New2CellNumb << "\t\t" << EntropyIncreaseList.at(New2CellNumb) << endl;

        // Then all the corresponding maps chain
        S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 1 (special) instead of 0 (ordinary) in the State Faces vector

        /// Add the new element to s_faces_sequence if it is still not here
        if(find(s_faces_sequence.begin(),s_faces_sequence.end(),New2CellNumb) == s_faces_sequence.end())
            s_faces_sequence.push_back(New2CellNumb);

        EntropyIncreaseList.at(New2CellNumb) = 0.0; ///working_max
///min        EntropyIncreaseList.at(New2CellNumb) = pow(10,10);
//REPAIRS       cout << New2CellNumb << "  " << EntropyIncreaseList.at(New2CellNumb) << endl; //for(auto kl :s_faces_sequence) cout << " " <<kl ;
//INFO cout << "\tFace fraction = \t" << s_faces_sequence.size() / (double) CellNumbs.at(2) << endl; //exit(10);

        /// Special and Ordinary Faces fraction calculation
        special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
        ordinary_faces_fraction = 1.0 - special_faces_fraction;

    } while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process
//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;

/// Closing and deleting
    TJsTypes.clear();
    TJsTypes.shrink_to_fit();
    EntropyIncreaseList.clear();
    EntropyIncreaseList.shrink_to_fit();
    new_neigh_TJs.clear();
    new_neigh_TJs.shrink_to_fit();
    new_neigh_Faces.clear();
    new_neigh_Faces.shrink_to_fit();
    FES.makeCompressed();
    AFS.makeCompressed();

    return 0;
} /// end of Processing_ipIndex()
*/

/// (5) DDRX process
/*!
 *
 * @param State_Vector
 * @param special_faces_sequence
 * @return
 */
//int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence, double max_sFaces_fraction, std::vector<char*> const paths, int number_of_types, std::vector<unsigned int> &CellNumbs) {
int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence) {
///=============================================================================================================================================////
///========================================================================= 'D' =================================================================////
/// ==================================================================>  DDRX process   <========================================================////
///=============================================================================================================================================////
    /// Function
    tuple<double, double, double> find_aGBseed(unsigned int Facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates);
    /// Streams
    ofstream NewSeedsStream;
    NewSeedsStream.open(paths.at(8), ios::trunc);
    NewSeedsStream << "New generated seeds\t" << endl;
    NewSeedsStream.close();

    /// Vector
//    vector<Triplet<double>> AllSeeds_coordinates, NewSeeds_coordinates;
    vector<tuple<double, double, double>> AllSeeds_coordinates;
    tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);

    /// Face Seeds reading from file Seeds.txt
    AllSeeds_coordinates = TuplesReader(paths.at(7));

/// ============== Constants and model parameters ===================
    double kb = 1.23* pow(10,-23);
    double GSav = 200.0*pow(10,-9), Temp = 300, Stress = 4.0*pow(10.0,6), Rcr = 0.5*GSav/10.0,
            Q_energy = 213.0*pow(10,1), R_const = 8.31;
/// Conversion probability
    double prob_seed = exp(-Q_energy/(R_const*Temp));
    //    double prob_seed = exp(- Stress*(13.0*pow(Rcr,3)/3.0)/(kb*Temp));
    cout << "DDRX module:: prob_seed =\t" << prob_seed << endl;

    srand((unsigned) time(NULL)); // seed for random
    // Output stream opening
    NewSeedsStream.open(paths.at(8), ios::app);

    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber) {
        double rv = (rand() / (RAND_MAX + 1.0));
        //   cout << rv << "\t" << prob_seed << endl;
        if (rv <= prob_seed) {
            NewSeed_coordinates = find_aGBseed(fnumber, paths, CellNumbs, AllSeeds_coordinates);
            NewSeedsStream << fnumber << "\t" << get<0>(NewSeed_coordinates) << "\t" << get<1>(NewSeed_coordinates) << "\t" << get<2>(NewSeed_coordinates) << endl;
        }
    }
    NewSeedsStream.close();

    return 0;
}

std::vector <unsigned int> Smax_sequence_reader(char* SFS_dir) {
    std::vector <unsigned int> s_faces_sequence;
    cout << "-------------------------------------------------------f----------------------------------------------------------------------"s << endl;
    cout << "Warning!!: special_faces_sequence successfully loaded from file:\t"s << SFS_dir << " because the Processing is OFF"s << endl;
    cout << "------------------------------------------------------------------------------------------------------------------------------"s << endl;
    s_faces_sequence = VectorIReader(SFS_dir); //all Faces

    return s_faces_sequence;
// REPAIR   cout << "S_size =\t" << special_faces_sequence.size() << endl;

}

/// Log-normal distribution generator
std::vector<double> Log_normal_distribution (double &mu_f, double &sigm_f, int baskets) { // Log-normal distribution generator
    vector<double>  s_lenght_distribution;
    for (int x = 1; x < baskets; x++) {
        double X = 1.0 * x; // scale
        s_lenght_distribution.push_back((1.0 / (X * sigm_f * sqrt(2.0 * 3.14159))) * exp(-pow((log(X) - log(mu_f)), 2) / (2.0 * pow(sigm_f, 2))));
    }

    s_lenght_distribution.at(0) += 1.0 - std::accumulate(s_lenght_distribution.begin(), s_lenght_distribution.end(), decltype(s_lenght_distribution)::value_type(0));
    //cout << std::accumulate(s_lenght_distribution.begin(), s_lenght_distribution.end(), decltype(s_lenght_distribution)::value_type(0)) << endl;

    return s_lenght_distribution; // end of Stips_distribution function
}


/// Heap
/*
 * int Processing_ExperimentalData(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs, std::vector<char*> paths) {
    if ( ProcessingON(confpath, time_step_one) && *Processing_type == 'E') {
        string  pass1path = input_dir + "pass_1_misorientation.txt"s; char* p1path = const_cast<char*>(pass1path.c_str());
        string  pass1_B2path = input_dir + "pass_1_B2.txt"s; char* p1B2path = const_cast<char*>(pass1_B2path.c_str());
        string  pass2path = input_dir + "pass_2_misorientation.txt"s; char* p2path = const_cast<char*>(pass2path.c_str());
        string  pass2_B2path = input_dir + "pass_2_B2.txt"s; char* p2B2path = const_cast<char*>(pass2_B2path.c_str());
        string  pass4path = input_dir + "pass_4_misorientation.txt"s; char* p4path = const_cast<char*>(pass4path.c_str());
        string  pass4_B2path = input_dir + "pass_4_B2.txt"s; char* p4B2path = const_cast<char*>(pass4_B2path.c_str());
        string  pass8path = input_dir + "pass_8_misorientation.txt"s; char* p8path = const_cast<char*>(pass8path.c_str());
        string  pass8_B2path = input_dir + "pass_8_B2.txt"s; char* p8B2path = const_cast<char*>(pass8_B2path.c_str());

        vector<*char> e_paths = {p1path, p1B2path, p2path, p2B2path, p4path, p4B2path, p8path, p8B2path};

    }
    return 0;
}
*/