#ifndef PCC_PROCESSING_DESIGN_CORROSION_LAB_H
#define PCC_PROCESSING_DESIGN_CORROSION_LAB_H

/// ======# 1 #================= void surface_interface_corrosion() function ==============================================================///

/*!
 * @brief simulation of the degradation of interfaces/grain boundariesstarting from one of the PCC edges
 * @param p_cells_history
 */
std::vector<std::vector<double>> surface_interface_corrosion(Config &kinetics_configuration, CellDesign &processing_cells_design, std::vector<CellEnergies> &gb_energies);

/*!
 * @brief simulation of the degradation of interfaces/grain boundaries starting from the macrocrack represented as a half-plane cut Subcomplex
 * @param p_cells_history
 */
std::vector<std::vector<double>> macrocrack_interface_corrosion(Config &kinetics_configuration, CellDesign &processing_cells_design, std::vector<Subcomplex> &plane_crack_subs, std::vector<CellEnergies> &gb_energies);


#endif //PCC_PROCESSING_DESIGN_CORROSION_LAB_H
