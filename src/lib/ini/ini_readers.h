#ifndef PCC_PROCESSING_DESIGN_INI_READERS_H
#define PCC_PROCESSING_DESIGN_INI_READERS_H

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project
std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string &main_type, std::string &e_type);

void config_reader_processing(std::string &source_path, std::vector<std::string> &sequence_source_paths, std::vector<std::vector<double>> &max_fractions_vectors, std::vector<std::vector<double>> &max_cfractions_vectors, double &mu, double &sigma, std::vector<std::string> &ptype_vector, std::vector<std::string> &ctype_vector, std::vector<double> &pindex_vector, std::ofstream &Out_logfile_stream);

void config_reader_writer(std::string &source_path, std::vector<int> &writer_specifications, std::ofstream &Out_logfile_stream);

#endif //PCC_PROCESSING_DESIGN_INI_READERS_H
