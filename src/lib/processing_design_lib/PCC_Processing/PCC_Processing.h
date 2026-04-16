///================================ PCC Processing module ======================================================================================///
///============================================================================================================================================///
///* The interface use functions from PCC_Processing/functions C++ libraries to generate quasi-random or non-random processes.           *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///

#ifndef PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
#define PCC_PROCESSING_DESIGN_PCC_PROCESSING_H

/*!
 * @brief Processing module :: creates design vectors for all k-cells in a PCC and saves them into a CellDesign object.
 * @param configuration // Object contained initial configuration
 * @param pcc_subcomplexes // vector of Subcomplexes
 * @param new_cells_energies // Cell Energies object
 * @return CellDesign object
 */
CellDesign PCC_Processing(Config &configuration, std::vector<Subcomplex> &pcc_subcomplexes, std::vector<CellEnergies> &new_cells_energies);

#endif //PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
