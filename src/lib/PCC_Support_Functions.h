#ifndef PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H

/// Bool check if the file 'fileName' (path) exists in the directory
bool is_file_exists(const std::string fileName);

/// Creates int std::vector from the file (containing in the directory) 'FilePath'
std::vector<unsigned int> VectorIReader(const char* FilePath); // creation int std::Vector from file

/// Creates double std::vector from the file (containing in the directory) 'FilePath'
std::vector<double> VectorDReader(const char* FilePath);

#endif //PCC_PROCESSING_DESIGN_PCC_SUPPORT_FUNCTIONS_H