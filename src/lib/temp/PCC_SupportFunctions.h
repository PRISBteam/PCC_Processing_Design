/// Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols);

/// Function for reading triplet list from file
vector<Triplet<double>> TripletsReader(char* SMpath);

/// Function for reading tuples list from file
vector<vector<double>> VectorVectors4Reader(char* SMpath);

vector<tuple<double, double, double>> TuplesReader(char* SMpath);

vector<tuple<double, double, double>> dTuplesReader(char* SMpath, unsigned int &size);

/// Creation int std::Vector from file
std::vector<unsigned int> VectorIReader(char* FilePath);

/// Creation double std::Vector from file
std::vector<double> VectorDReader(char* FilePath);

/// * The function count the number of Edges possessing with types J0, J1, J2, J3 and J4 for every junction * ///
/// * Calculation Face-Edge index ::                                                                                                     * ///
vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, vector<double> const& TJsTypes);

/// * Function calculates the vector<int> "TJsTypes" of types TJs in the PCC using its FES incidence matrix and special faces sequence (s_faces_sequence) * ///
/// *                                                                                                                                                    * ///
vector<double> EdgesTypesCalc(std::vector<unsigned int> const &CellNumbs, vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const &FES);

/// EntropyIncreaseList
vector<double> Get_EntropyIncreaseList(std::vector<unsigned int> &S_Vector, vector<double> const &TJsTypes, SpMat const &FES);

std::vector<vector<int>> Get_cases_list(std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat const &FES, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, double const &p_index);

std::vector<double> dq_upleft(std::vector<double> const &grain_q_coord, double dq_step);

std::vector<double> dq_upright(vector<double> const &grain_q_coord, double dq_step);

std::vector<double> dq_downleft(vector<double> const &grain_q_coord, double dq_step);

std::vector<double> dq_downright(vector<double> const &grain_q_coord, double dq_step);

/// Crystallographic BCC disorientations
double Get_2grains_FCCdisorientation(std::vector<double> &grain_quaternion1, std::vector<double> &grain_quaternion2);

/// isHAGB checker
bool isHAGB(std::vector<double> &grain_quaternion, std::vector<double> &new_grain_quaternion, double &threshold);

std::vector<vector<int>> Get_crystallographic_cases_list(std::vector<vector<double>> &grain_quaternions_list, std::map<unsigned int, std::vector<unsigned int>> &g_gbs_map, double &dq, double &HAGBs_threshold, std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat &GFS, SpMat &FES, std::map<unsigned int, unsigned int> &cases_to_grains, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, std::map<unsigned int, std::vector<double>> &cases_to_new_quaternions, double const &p_index);

std::vector<vector<int>> Get_crystallographic_cases_random_list(std::vector<vector<double>> &grain_quaternions_list, std::map<unsigned int, std::vector<unsigned int>> &g_gbs_map, double &dq, double &HAGBs_threshold, std::vector<int> const &S_Vector, std::vector<int> const &EdgeTypes, SpMat &GFS, SpMat &FES, std::map<unsigned int, unsigned int> &cases_to_grains, std::map<unsigned int, std::vector<unsigned int>> &cases_to_sfaces, std::map<unsigned int, std::vector<double>> &cases_to_new_quaternions, double const &p_index);

/// DDRX support function :: GFS matrix reading and calculation of new seeds at the centres of GBs
tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates);

bool is_file_exists(const string fileName);

/// Configuration TJs entropy
double Get_TJsEntropy(vector<unsigned int> special_faces_seq);

std::vector<int> state_vector_by_sequence(std::vector<unsigned int> const &cell_sequence, int cell_type);

/// ARCHIVE ///
//Erase First Occurrence of given  substring from main string
void eraseSubStr(std::string & mainStr, const std::string & toErase);