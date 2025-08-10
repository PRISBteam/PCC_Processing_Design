///================================ PCC Design module ======================================================================================///
///============================================================================================================================================///
///* The interface use functions from Design_<***>_functions.h C++ libraries for optimisation of the State Vectors           *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///

#ifndef PCC_PROCESSING_DESIGN_PCC_DESIGN_H
#define PCC_PROCESSING_DESIGN_PCC_DESIGN_H
/*!
 * @brief Design module :: creates list of State Vectors optimised according to a goal function
 * @param configuration
 * @return std::vector<std::vector<int>> design_list_of_vectors
 */
std::vector<std::vector<int>> PCC_Design(Config &configuration);

#endif //PCC_PROCESSING_DESIGN_PCC_DESIGN_H
