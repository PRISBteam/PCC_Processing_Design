#ifndef PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNMENT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNMENT_FUNCTIONS_H

/// Random generation machine for a new 2-Cell number::
/// (1.1) Quasi-random choice of the element with # New2CellNumb from the list of numbers from 0 to OCellsNumb.
/*!
 * @param OCellsNumb :: the number of 'ordinary' faces
 * @return New2CellNumb :: a single k-cell number corresponding to the new k-cell type
 */
unsigned int NewCellNumb_R(unsigned int OCellsNumb);

std::vector<int> NewFacesStrip_RW( int iniFaceNumber, int strip_length, int Leap_friquency = 1, double Leap_dist = 1);

/// (1.2) Quasi-random choice of the element with # New2CellNumb from the list of numbers {0 to OCellsNumb}
/*!
 * @param <int> IniFaceNumber // initial number of face
 * @param <int> strip_length // the number of faces in one strip
 * @param <int> Leap_probability = 0 and <double> Leap_dist = 1 // for the future possibility of leaps: probability (hence frequency) and distance ((?)fraction of the complex size in grains)
 * @return vector<int> NewFacesStrip_RW // vector of faces joining in one strip
 */
// strip_sface_distribution - is a vector containing the numbers of strips of each length starting with 1 expressed in # of faces
// Example: vector<int> strip_sface_distribution = {2 4 5 27 8 6 3 1} means 2 strips of length 1 faces each, 4 strips of length 2 faces each,... , 1 strip of length 8 faces each
// int strip_length - is the number of faces in the current strip or RW path length in # of faces

/// (2.1) The Random generation process function
/*!
 * @param S_Vector
 * @param s_faces_sequence
 * @param max_sFaces_fraction
 */
// S_Vector with its non-zero elements set any pre-define structure of special element feeding to the function Processing_Random
std::vector<unsigned int> Processing_Random(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> max_fractions_vectors);

/// (2.2) Lengthy special strips with the distribution of lengths taken from file
/*!
 * @param face_strip_distribution
 * @param S_Vector
 * @param s_faces_sequence
 * @return
 */
 std::vector<std::vector<int>> RStrips_Distribution( std::vector<unsigned int> const &face_strip_distribution, std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence);

#endif //PCC_PROCESSING_DESIGN_PROCESSING_ASSIGNMENT_FUNCTIONS_H
