#ifndef PCC_PROCESSING_DESIGN_IRRADIATION_LAB_H
#define PCC_PROCESSING_DESIGN_IRRADIATION_LAB_H

/// ======# 1 #================= void interface_irradiation_damage() function ==============================================================///

/*!
 * @brief simulation of interface damage due to particle beam irradiation
 * @param p_cells_history
 */
std::vector<std::vector<double>>  interface_irradiation_damage(Config &config, Material &material, CellDesign &processing_cells_design);

#endif //PCC_PROCESSING_DESIGN_IRRADIATION_LAB_H
