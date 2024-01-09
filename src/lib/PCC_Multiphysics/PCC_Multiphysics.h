///================================ PCC Multiphysics module =================================================================================///
///=========================================================================================================================================///

#ifndef PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H
#define PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H

/// Multiphysics module :: creates energy vectors for all k-cells in a PCC and write them to CellEnergies object.
std::vector <double> PCC_Multiphysics(Macrocrack &large_crack, double external_vonMizes_stress);


#endif //PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H
