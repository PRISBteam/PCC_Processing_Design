///================================ PCC Multiphysics module ====================================================================================///
///============================================================================================================================================///
///* The interface use functions from Multiphysics_<***>_functions.h C++ libraries to set physical scale and the corresponding cell energies.*///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///

#ifndef PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H
#define PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H

/*!
 * @brief Multiphysics module :: set physical scale and the corresponding k-cell energies in a PCC saving them into CellEnergies object.
 * @param configuration
 * @param pcc_subcomplexes
 * @return CellEnergies object
 */
std::vector<CellEnergies> PCC_Multiphysics(Config &configuration, std::vector<Subcomplex> &pcc_subcomplexes);
std::vector<CellEnergies> PCC_Multiphysics(Config &configuration, std::vector<Subcomplex> &pcc_subcomplexes, std::vector<Macrocrack> &crack_growth_series);

#endif //PCC_PROCESSING_DESIGN_PCC_MULTIPHYSICS_H