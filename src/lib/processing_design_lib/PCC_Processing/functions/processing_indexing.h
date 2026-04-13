
#ifndef PCC_PROCESSING_DESIGN_PROCESSING_INDEXING_H
#define PCC_PROCESSING_DESIGN_PROCESSING_INDEXING_H

/*!
 * @brief
 * @param k_cell_type
 * @param Configuration_State
 * @return
 */
std::vector<unsigned int> TopDown_cell_indexing(int k_cell_type, std::vector<std::vector<unsigned int>> &Configuration_State);

std::vector<unsigned int> TopDown_cell_indexing(int k_cell_type, std::vector<unsigned int> &higher_order_cells_state_vector);

#endif // PCC_PROCESSING_DESIGN_PROCESSING_INDEXING_H
