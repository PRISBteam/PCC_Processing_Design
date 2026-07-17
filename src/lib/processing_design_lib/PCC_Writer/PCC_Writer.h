///================================ PCC Writer module ===================================================================================///
#ifndef PCC_PROCESSING_DESIGN_PCC_WRITER_H
#define PCC_PROCESSING_DESIGN_PCC_WRITER_H

/// PCC Writer
/*! # 1 #
 * @brief PCC_Writer(CellDesign &new_cells_design)
 * @param new_cells_design
 */
void PCC_Writer(CellDesign &new_cells_design);

/*! # 2 #
 * @brief Overloaded PCC_Writer(..) with the input of several classes, such as pcc_subcomplexes, new_cells_energies, new_cells_design, pcc_processed
 * @param pcc_subcomplexes
 * @param new_cells_energies
 * @param new_cells_design
 * @param pcc_processed
 */
void PCC_Writer(std::vector<Subcomplex> &pcc_subcomplexes, std::vector<CellEnergies> &new_cells_energies, CellDesign &new_cells_design, ProcessedComplex &pcc_processed);

#endif //PCC_PROCESSING_DESIGN_PCC_WRITER_H
