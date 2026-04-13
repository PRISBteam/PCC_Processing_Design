///================================== PCC Kinetics module =======================================================================================///
///=============================================================================================================================================///
///* The interface use functions from PCC_Kinetics/functions C++ libraries to generate for each p-cell in a PCC the moment of 'time'          *///
///* in the [0,1] range when it changed its special 'generated' type as the result of a 'kinetic' process (e.g. corrosion or irradiation)    *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2025 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///

#ifndef PCC_PROCESSING_DESIGN_PCC_KINETICS_H
#define PCC_PROCESSING_DESIGN_PCC_KINETICS_H

/*!
 * @brief Kinetics module ::generate for each p-cell in a PCC the moment of 'time' in the [0,1] range when it changed its special 'generated' type as the result of a 'kinetic' process (e.g. corrosion or irradiation)
 * @param configuration
 * @return std::vector<double> p_cells_history
 */
std::vector<std::vector<double>> PCC_Kinetics(Config &kinetics_configuration, CellDesign &processing_cells_design, std::vector<Subcomplex> &pcc_sub, std::vector<CellEnergies> &gb_energies);

#endif //PCC_PROCESSING_DESIGN_PCC_KINETICS_H
