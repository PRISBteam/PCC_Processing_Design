#ifndef PCC_PROCESSING_DESIGN_INI_READERS_H
#define PCC_PROCESSING_DESIGN_INI_READERS_H

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project
// === # 1 # === //
std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string & cell_complex_standard, std::string &main_type);
// === # 2 # === //
void config_reader_processing(std::string &source_path, std::vector<std::string> &sequence_source_paths, std::vector<std::vector<double>> &max_fractions_vectors, std::vector<std::vector<double>> &max_cfractions_vectors, double &mu, double &sigma, unsigned int &bins_numb, std::vector<std::string> &ptype_vector, std::vector<std::string> &ctype_vector, std::vector<double> &pindex_vector, std::ofstream &Out_logfile_stream);
// === # 3 # === //
std::vector<double> config_reader_characterisation(std::string const &source_path, std::vector<int> &charlabs_polyhedrons, std::vector<int> &charlabs_faces, std::vector<int> &charlabs_edges, std::vector<int> &charlabs_nodes, std::vector<int> &charlabs_laplacians, std::ofstream &Out_logfile_stream);
// === # 4 # === //
void config_reader_writer(std::string &source_path, std::vector<int> &writer_specifications, std::ofstream &Out_logfile_stream);
// === # 5 # === //
void config_reader_subcomplex(std::string &source_path, std::string &sctype, std::vector<double> &plane_orientation, double &cut_length, unsigned int &grain_neighbour_orders, std::ofstream &Out_logfile_stream);
// === # 6 # === //
void config_reader_multiphysics(std::string &source_path, std::string &Mid_matrix, std::string &Mid_inclusion, std::tuple<double, double, double> &sample_dimensions, double &tau, Eigen::MatrixXd &external_stress, std::vector<double> &macrocrack_ini, std::ofstream &Out_logfile_stream);

#endif //PCC_PROCESSING_DESIGN_INI_READERS_H
