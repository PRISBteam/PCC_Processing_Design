#ifndef PCC_PROCESSING_DESIGN_CORROSION_LAB_H
#define PCC_PROCESSING_DESIGN_CORROSION_LAB_H

/// ======# 1 #================= void surface_interface_corrosion() function ==============================================================///

/*!
 * @brief simulation of interface degradation starting from one of the PCC edges
 * @param p_cells_history
 */
std::vector<std::vector<double>>  surface_interface_corrosion(Config &kinetics_configuration, Material &material, CellDesign &processing_cells_design);

#endif //PCC_PROCESSING_DESIGN_CORROSION_LAB_H
