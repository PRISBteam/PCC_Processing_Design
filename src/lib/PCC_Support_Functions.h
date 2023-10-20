#ifndef PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H

/// # 1 # Bool check if the file 'fileName' (path) exists in the directory
bool is_file_exists(const std::string fileName);

/// # 2 # Creates int std::vector from the file (containing in the directory) 'FilePath'
std::vector<unsigned int> VectorIReader(const char* FilePath); // creation int std::Vector from file

/// # 3 # Creates double std::vector from the file (containing in the directory) 'FilePath'
std::vector<double> VectorDReader(const char* FilePath);

/// # 4 # Creation Eigen::Sparse_Matrix from file
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols);

#endif //PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H