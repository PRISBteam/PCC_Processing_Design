///================================ PCC Processing module ===================================================================================///

#ifndef PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
#define PCC_PROCESSING_DESIGN_PCC_PROCESSING_H

/// Processing module :: creates design vectors for all k-cells in a PCC and write them to CellDesign object.
CellsDesign PCC_Processing(std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<int>> &Configuration_cState);

#endif //PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
