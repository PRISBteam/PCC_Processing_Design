///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new "secondary" identifications of k-Cells   **/
/**  in the PCC Processing module. It makes them "fractured" and can go alongside their special or ordinary primal types.         **/
///==============================================================================================================================///

#ifndef PCC_PROCESSING_DESIGN_PROCESSING_INDUCED_LABELLING_H
#define PCC_PROCESSING_DESIGN_PROCESSING_INDUCED_LABELLING_H

/*!
 * @brief Kinematic (with adhesion energies, but without forces) generation algorithm of the induced ('fractured') k-cell indexing
 * @param cell_type
 * @param s_faces_sequence
 * @param Configuration_cState
 * @param max_cfractions_vectors
 * @return
 */
std::vector <unsigned int> PCC_Kinematic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<unsigned int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors);

/*!
 * @brief
 * @param cell_type
 * @param s_faces_sequence
 * @param Configuration_cState
 * @param max_cfractions_vectors
 * @param inclusion_sintensity_factor
 * @param crack_sintensity_factor
 * @return
 */
 std::vector <unsigned int> PCC_Kineic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<unsigned int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors, double inclusion_sintensity_factor, double crack_sintensity_factor);
// std::vector <unsigned int> PCC_Kineic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors);

/*!
 * @brief
 * @param cell_type
 * @param s_faces_sequence
 * @param Configuration_cState
 * @param max_cfractions_vectors
 * @param inclusion_sintensity_factor
 * @param crack_sintensity_factor
 * @param new_cells_energies
 * @param Mat_id
 * @return
 */
std::vector <unsigned int> PCC_Kinematic2_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<unsigned int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors, double inclusion_sintensity_factor, double crack_sintensity_factor, CellEnergies &new_cells_energy, std::vector<Agglomeration> &agglomerations_vector, Material &Mat_id);

//std::vector <unsigned int> PCC_Kineic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<unsigned int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors, double inclusion_sintensity_factor, double crack_sintensity_factor, std::vector<CellEnergies> &new_cells_energies, Material &Mat_id);

#endif // PCC_PROCESSING_DESIGN_PROCESSING_INDUCED_LABELLING_H
