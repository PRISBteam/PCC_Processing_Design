#ifndef PCC_PROCESSING_DESIGN_INI_READERS_H
#define PCC_PROCESSING_DESIGN_INI_READERS_H

/// Tailored Reader for the main.ini file in the ../config/ subdirectory of the project

#include "../../../src/lib/processing_design_lib/PCC_Objects.h"

// === # 0 # === //
/*!
 * @brief config_reader_main :: read input parameters from the project file config/main.ini necessary for the code execution.
 * @param main_ini_data
 * @return std::vector<int> contained (0) PCC dimension and all the modules ID: 1 - if ON, and 0 if OFF in the config/main.ini file;
 * (1) isSubcomplexON, (2) isMultiphysicsON, (3) isProcessingON, (4) isCharacterisationON, (5) isDesignON, (6) isWriterON
 */
std::vector<int> config_reader_main(Config &main_config);
// === # 1 # === //
/*!
 * @brief config_reader_subcomplex :: read input parameters from the project file config/subcomplex.ini necessary for the code execution.
 * @param sctype
 * @param plane_orientation
 * @param cut_length
 * @param grain_neighbour_orders
 * @return void
 */
//void config_reader_subcomplex(std::string &sctype, std::vector<double> &plane_orientation, double &cut_length, unsigned int &grain_neighbour_orders, bool &is_log_file);
void config_reader_subcomplex(Config &subcomplex_config);
// === # 2 # === //
/*!
 * @brief config_reader_multiphysics :: read input parameters from the project file config/multiphysics.ini necessary for the code execution.
 * @param Mid_matrix
 * @param Mid_inclusion
 * @param sample_dimensions
 * @param tau
 * @param external_stress
 * @param macrocrack_ini
 * @return void
 */
void config_reader_multiphysics(Config &multiphysics_config);

// === # 3 # === //
/*!
 * @brief config_reader_processing :: read input parameters from the project file config/processing.ini necessary for the code execution.
 * @param sequence_source_paths
 * @param max_fractions_vectors
 * @param max_cfractions_vectors
 * @param mu
 * @param sigma
 * @param bins_numb
 * @param ptype_vector
 * @param ctype_vector
 * @param pindex_vector
 * @return void
 */
void config_reader_processing(Config &processing_config);

// === # 4 # === //
/*!
 * @brief config_reader_characterisation :: read input parameters from the project file config/characterisation.ini necessary for the code execution.
 * @param charlabs_polyhedrons
 * @param charlabs_faces
 * @param charlabs_edges
 * @param charlabs_nodes
 * @param charlabs_laplacians
 * @return std::vector<double> config_characterisation_vector contained specifications for calculation & output of various characteristics (1 - ON, 0 - OFF)
 */
std::vector<double> config_reader_characterisation(std::vector<int> &charlabs_polyhedrons, std::vector<int> &charlabs_faces, std::vector<int> &charlabs_edges, std::vector<int> &charlabs_nodes, std::vector<int> &charlabs_laplacians, bool &is_log_file);

// === # 5 # === //
/*!
 * @brief config_reader_design :: read input parameters from the project file config/design.ini necessary for the code execution.
 * @param
 * @return void
*/
void config_reader_design(bool &is_log_file);

// === # 6 # === //
/*!
 * @brief config_reader_writer :: read input parameters from the project file config/writer.ini necessary for the code execution.
 * @param writer_specifications
 * @return void
 */
void config_reader_writer(std::vector<int> &writer_specifications, bool &is_log_file);

#endif //PCC_PROCESSING_DESIGN_INI_READERS_H
