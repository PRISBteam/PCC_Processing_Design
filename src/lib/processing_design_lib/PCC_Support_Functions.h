#ifndef PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H

#include <Eigen/SparseCore>
#include <set>

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

/// # 1 # Bool check if the file 'fileName' (path) exists in the directory
bool is_file_exists(const std::string fileName);

/// # 2 # Creates int std::vector from the file (containing in the directory) 'FilePath'
std::vector<unsigned int> VectorIReader(const char* FilePath); // creation int std::Vector from file

/// # 3 # Creates double std::vector from the file (containing in the directory) 'FilePath'
std::vector<double> VectorDReader(const char* FilePath);

/// # 4 # Creation Eigen::Sparse_Matrix from file
//Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols);
Eigen::SparseMatrix<double> SMatrixReader(std::string SMpath, unsigned int Rows, unsigned int Cols);

/// # * # Creation Eigen::Sparse_Matrix from file
std::vector<std::tuple<double, double, double>> Tuple3Reader(std::string SMpath);

/// # * # Vector to Set simple converter
std::set<unsigned int> convertToSet(std::vector<unsigned int> &v);

/// # * # Set to Vector simple converter
std::vector<unsigned int> SetToVector(std::set<unsigned int> &v);


/*! ## 3 ##
 * @brief Log-normal distribution generator used for the RStrips_Distribution() processing function.
 * Its output vector contains the numbers of strips of each length starting with 1 expressed in # of faces.
 * Example: std::vector<int> strip_distribution = {2 4 5 27 8 6 3 1} means 2 strips of length 1 faces each, 4 strips of length 2 faces each,... , 1 strip of length 8 faces each.
 * @param mu_f
 * @param sigm_f
 * @param bins_number
 * @return strip_distribution_vector
 */
std::vector<double> Log_normal_distribution (double mu_f, double sigm_f, int bins_number);

/// # 5 # Finding barycenter coordinates as a tuple<double, double, double> for a given 'facenumb' face
std::tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<std::string> const &paths, std::vector<unsigned int> const &CellNumbs, std::vector<std::tuple<double, double, double>> const &AllSeeds_coordinates);

/// # 5 # Finding barycenter coordinates as a tuple<double, double, double> for a given 'edgenumb' edge
std::tuple<double, double, double> find_anEdgeSeed(unsigned int edgenumb, std::vector<std::string> const &paths, std::vector<unsigned int> const &CellNumbs, std::vector<std::tuple<double, double, double>> const &AllSeeds_coordinates);


template <typename TP> inline constexpr
int sign(TP x, std::false_type is_signed);
template <typename TP> inline constexpr
int sign(TP x, std::true_type is_signed);
template <typename TP> inline constexpr
int sign(TP x);

/*!
 * @brief  Face indexing based on the preliminsry indrxing of their edges (vector<double> TJsTypes)
 * @param face_number
 * @param FES
 * @param TJsTypes
 * @return
 */
std::vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, std::vector<double> const& TJsTypes);

/*!
 * @brief  Node indexing based on the preliminsry indrxing of their edges (vector<double> TJsTypes)
 * @param node_number
 * @param ENS
 * @param TJsTypes
 * @return
 */
std::vector<double> NodeIndex(unsigned int node_number, Eigen::SparseMatrix<double> const& ENS, std::vector<double> const& TJsTypes);

/*!
 * @brief Node Types Calculation function
 * @param CellNumbs
 * @param s_faces_sequence
 * @param ENS
 * @return
 */
std::vector<double> NodesTypesCalc(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &ENS);

/*!
 * @brief  Edge Types Calculation function
 * @param CellNumbs
 * @param s_faces_sequence
 * @param FES
 * @return
 */
std::vector<double> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES);

/*!
 * @brief edge Configuration Entropy measure
 * @param special_faces_seq
 * @return
 */
double Get_TJsEntropy(std::vector<unsigned int> &special_faces_seq);
double Get_TJsEntropy(std::vector<unsigned int> &special_faces_seq, std::vector<double> &Jtypes, std::vector<double> &Dtypes);

/*!
 * @brief
 * @param S_Vector
 * @param EdgeTypes
 * @param FES
 * @param cases_to_sfaces
 * @param p_index
 * @return
 */
std::vector<std::vector<int>> Get_cases_list(std::vector<unsigned int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat const &FES, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, double const p_index = 0);

/*!
 * @brief
 * @param face_sequence
 * @return
 */
std::vector<std::tuple<double, double, double>>  face_sequence_barycentre_coordinates(std::vector<unsigned int> &face_sequence);
std::vector<std::tuple<double, double, double>>  face_sequence_barycentre_coordinates(std::vector<unsigned int> &sfaces_set, std::vector<std::tuple<double, double, double>> &all_face_coordinates);

/*!
 * @brief
 * @param sfaces_set
 * @return
 */
std::vector<std::tuple<double, double, double>>  face_sequence_barycentre_coordinates(std::set<unsigned int> &sfaces_set);
std::vector<std::tuple<double, double, double>>  face_sequence_barycentre_coordinates(std::set<unsigned int> &sfaces_set, std::vector<std::tuple<double, double, double>> &all_face_coordinates);

/*!
 * @brief cout for a vector of type 'unsigned int'
 * @param vector
 * @param text
 */
void Vector_ui_cout(std::vector <unsigned int> &vector, std::string text);

/*!
 * @brief cout for a vector of vectors of type 'unsigned int'
 * @param vector
 * @param text
 */
void Vector_of_vectors_ui_cout(std::vector<std::vector <unsigned int>> &vector, std::string text);

/*!
 * @brief Variable get<var> for a std::tuple<double, double, double>
 * @param direction
 * @param t
 * @return
 */
double get_i(int &direction, std::tuple<double, double, double> &t);

std::vector<std::tuple<double, double, double>> kCell_barycentre_coordinates(int k_type, std::vector<unsigned int> &kcell_sequence);
std::vector<std::tuple<double, double, double>> kSequence_barycentre_coordinates(int k_type, std::vector<unsigned int> &kcell_sequence);
std::vector<std::tuple<double, double, double>> kFaceSeq_barycentre_coordinates(int k_type, std::vector<unsigned int> &kFace_sequence);

/*!
 * @brief Shuffle (coordinates) data rows inside the file
 * @param name_infile
 * @param name_outfile
 * @param rows_number
 */
void shuffle_coord_vector(std::string name_infile, std::string name_outfile, unsigned int rows_number);

#endif //PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H