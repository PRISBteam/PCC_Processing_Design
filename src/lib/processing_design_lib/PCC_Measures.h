
std::vector<int> Edge_types_byFaces(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &special_face_sequence, std::vector<double> &j_fractions, std::vector<double> &d_fractions);

std::vector<double> j_fractions_vector(std::vector<int> const &TJsTypes); // based on Edges vector

double Face_edge_index(std::vector<unsigned int> &special_face_sequence, Eigen::SparseMatrix<double> const& FES, double norm_const = 1.0); // based on Edges vector

double Node_edge_index(std::vector<unsigned int> &special_face_sequence, Eigen::SparseMatrix<double> const& ENS, double norm_const = 1.0); // based on Edges vector

std::vector<double> d_fractions_vector(std::vector<double> const &TJsTypes); // based on Edges vector

double Configuration_Entropy(std::vector<double> const &TJsTypes);

double Configuration_Entropy(std::vector<double> const &j_fractions); // 2nd  overloaded function

std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions);

double Configuration_Entropy_change(std::vector<int> const &TJsTypes, std::vector<int> const &TJsTypesNew);

double Configuration_Entropy_change(std::vector<double> const &j_fractions, std::vector<int> const &TJsTypesNew); // 2nd overloaded function