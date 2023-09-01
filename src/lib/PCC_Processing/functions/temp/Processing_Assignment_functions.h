///================================ A part of the PCC Processing module HEADER file =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifiers for k-Cells, k={0,1,2,3}     **/
/**  in the PCC Processing module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                        **/
///==============================================================================================================================///

/// (1) Quasi-random choice of the element with # New2CellNumb from the list of numbers {0 to OCellsNumb}
/*!
 * @param OCellsNumb // the number of ordinary faces
 * @return New2CellNumb
 */
unsigned int NewCellNumb_R(unsigned int OCellsNumb);

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
vector<int> NewFacesStrip_RW( int iniFaceNumber, int strip_length, int Leap_friquency = 1, double Leap_dist = 1); // Random generation machine for a strips of new 2-Cells

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
std::vector<unsigned int> Processing_Random(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> max_fractions_vectors);

/// (2) Lengthy special strips with the distribution of lengths taken from file
// vector<vector<int>> RStrips_Distribution( std::vector<unsigned int> const &face_strip_distribution, std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence) {

/// (3) Maximum Functional based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_maxFunctional(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

/// (3) Maximum Functional based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_minConfEntropy(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

/// (q) Maximum p crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_maxP_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

/// (q) Maximum p crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_Random_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

/// (y) Maximum Functional crystallographic based generation process
/*!
 * @param S_Vector
 * @param s_faces_sequence
 */
std::vector<unsigned int> Processing_maxF_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_State, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

/// (4) i(p) index govern processed
//int Processing_ipIndex(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, int index_type, double ip_index);

/// (5) DDRX process
/*!
 *
 * @param State_Vector
 * @param special_faces_sequence
 * @return
 */
int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence);

std::vector <unsigned int> Smax_sequence_reader(char* SFS_dir);

/// Log-normal distribution generator
std::vector<double> Log_normal_distribution (double &mu_f, double &sigm_f, int baskets); // Log-normal distribution generator
