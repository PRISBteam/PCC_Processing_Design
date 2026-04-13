#ifndef PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H
#define PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H

/*! ## 1 ##
 * @brief Mersenne Twister 19937 based algorithm implementing quasi-random generation function for a new k-cell number and many other purposes: quasi-random choice of an element from the list of numbers from 0 to OCellsNumb. OCellsNumb is the range of integer random numbers equal to the number of 'ordinary' k-cells.
 * @param OCellsNumb
 * @return random number
 */
unsigned int NewCellNumb_R(unsigned int OCellsNumb);

/*! ## 1.2 ##
 * @brief Basic quasi-random generation function for a new k-cell number and many other purposes: quasi-random choice of an element from the list of numbers from 0 to OCellsNumb. OCellsNumb is the range of integer random numbers equal to the number of 'ordinary' k-cells.
 * @param OCellsNumb
 * @return random number
 */
unsigned int NewCellNumb_R_fast(unsigned int OCellsNumb);

/*! ## 2 ##
 * @brief Mersenne Twister 19937 based algorithm implementing Random Walks generation process (with possible 'leaps') for a new k-cell chains of labels.
 * Leap_distance = 1 // for the future possibility of leaps: probability (hence frequency) and distance ((?)fraction of the complex size in grains).
 * @param cell_type
 * @param iniCellNumber
 * @param strip_length
 * @param Leap_frequency
 * @param Leap_distance
 * @return RW series of random numbers
 */
std::vector<unsigned int> NewCellsStrip_RW(int cell_type, unsigned int iniCellNumber, unsigned int strip_length, double Leap_frequency = 1, double Leap_distance = 1);

/*! ## 2.2 ##
 * @brief Basic Random Walks generation machine based on the 'rand()' function (with possible 'leaps') for a new k-cell chains of labels.
 * Leap_distance = 1 // for the future possibility of leaps: probability (hence frequency) and distance ((?)fraction of the complex size in grains).
 * @param cell_type
 * @param iniCellNumber
 * @param strip_length
 * @param Leap_frequency
 * @param Leap_distance
 * @return RW series of random numbers
 */
std::vector<unsigned int> NewCellsStrip_RW_fast(int cell_type, unsigned int iniCellNumber, unsigned int strip_length, double Leap_friquency = 1, double Leap_dist = 1);

/*! ## 3 ##
 * @brief Implement the quasi-random process of 'cell labelling', generating new structures on the PCC cells described by the corresponding State Vectors. Employ 'NewCellNumb_R' function.
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param multiplexity
 * @return random special_x_sequence
 */
std::vector<std::vector<unsigned int>> Processing_Random(int const cell_type, std::vector<std::vector<unsigned int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, bool multiplexity);

/*! ## 4 ##
 * @brief Implement the quasi-random process of 'cell labelling', generating new structures on the PCC cells described by the corresponding State Vectors. Unlike the 'Processing_Random' function it generates series special k-cells (series of strips or chains of k-cells) with the distribution of lengths taken from the 'face_strip_distribution' vector. Employs 'NewCellsStrip_RW' function.
 * @param cell_type
 * @param cell_strip_distribution
 * @param Configuration_State
 * @param max_fractions_vectors
 * @return special_x_sequence
 */
std::vector<std::vector<unsigned int>> Processing_Random_Strips(int cell_type, std::vector<unsigned int> &cell_strip_distribution, std::vector<std::vector<unsigned int>> &Configuration_State, std::vector<std::vector<double>> max_fractions_vectors);

/*! ## 5 ##
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return special_x_sequence
 */
///std::vector<unsigned int> Processing_maxFunctional(int cell_type, std::vector<std::vector<unsigned int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, bool multiplexity, double(*measure)(std::vector<double> const&));

/*! ## 6 ##
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return special_x_sequence
 */
///std::vector<unsigned int> Processing_maxF_crystallographic(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, double p_index);
// std::vector<unsigned int> Processing_maxP_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_sState, std::vector<vector<double>> const &max_fractions_vectors, double p_index);
// std::vector<unsigned int> Processing_Random_crystallographic(int cell_type, std::vector<vector<int>> &Configuration_sState, std::vector<vector<double>> const &max_fractions_vectors, double p_index);

#endif //PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNED_LABELLING_H

/*! ## 7 ##
 * @brief
 * @param cell_type
 * @param Configuration_State
 * @param max_fractions_vectors
 * @param p_index
 * @return special_x_sequence
 */
// std::vector<unsigned int> Processing_minConfEntropy(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, double p_index);

/// HEAP
///================================================================= 'R' =======================================================================////
/// ====================================================>  Random generation process  <========================================================////
///===========================================================================================================================================////

