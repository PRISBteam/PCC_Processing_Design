///================================ PCC Processing module ======================================================================================///
///============================================================================================================================================///
///* The interface use functions from Processing_<***>_functions.h C++ libraries to generate quasi-random or non-random processes.           *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///

#ifndef PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
#define PCC_PROCESSING_DESIGN_PCC_PROCESSING_H

/*!
 * @brief Processing module :: creates design vectors for all k-cells in a PCC and save them into CellDesign object.
 * @param configuration
 * @return CellDesign object
 */
CellDesign PCC_Processing(Config &configuration);

#endif //PCC_PROCESSING_DESIGN_PCC_PROCESSING_H
