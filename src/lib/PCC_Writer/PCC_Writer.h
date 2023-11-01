///================================ PCC Writer module ===================================================================================///
#ifndef PCC_PROCESSING_DESIGN_PCC_WRITER_H
#define PCC_PROCESSING_DESIGN_PCC_WRITER_H

/// PCC Writer
/*! # 1 #
 * @param new_cells_design
 */
void PCC_Writer(CellsDesign &new_cells_design);

/*! # 2 #
 * @param new_cells_design
 * @param pcc_processed
 */
void PCC_Writer(CellsDesign &new_cells_design, ProcessedComplex &pcc_processed);

#endif //PCC_PROCESSING_DESIGN_PCC_WRITER_H
