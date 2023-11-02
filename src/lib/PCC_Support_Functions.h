#ifndef PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H

#include <Eigen/SparseCore>

/// # 1 # Bool check if the file 'fileName' (path) exists in the directory
bool is_file_exists(const std::string fileName);

/// # 2 # Creates int std::vector from the file (containing in the directory) 'FilePath'
std::vector<unsigned int> VectorIReader(const char* FilePath); // creation int std::Vector from file

/// # 3 # Creates double std::vector from the file (containing in the directory) 'FilePath'
std::vector<double> VectorDReader(const char* FilePath);

/// # 4 # Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols);

/// # 5 # Finding barycenter coordinates as a tuple<double, double, double> for a given 'facenumb' face
std::tuple<double, double, double> find_aGBseed(unsigned int facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, std::vector<std::tuple<double, double, double>> & AllSeeds_coordinates);

#endif //PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H