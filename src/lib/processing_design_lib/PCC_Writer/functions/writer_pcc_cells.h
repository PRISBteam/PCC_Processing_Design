#ifndef PCC_PROCESSING_DESIGN_WRITER_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_WRITER_FUNCTIONS_H

/// # 1 # Sequences and Designs output
/*!
 * @brief Output of the special cell design sequences - in order of their appearance in a PCC during the processing route.
 * @param new_cells_design
 * @param output_counter
 */
    void PCC_CellSequences_Writer(CellDesign &new_cells_design, int &output_counter);

/// # 2 # Cell energies output
/*!
 * @brief Output of the variety of energies or potentials (elastic, thermal, etc.) associated with a PCC p-cells
 * @param new_cells_design
 * @param output_counter
 */
    void PCC_CellEnergies_Writer(std::vector<CellEnergies> &new_cell_energies, int &output_counter);

#endif //PCC_PROCESSING_DESIGN_WRITER_FUNCTIONS_H
