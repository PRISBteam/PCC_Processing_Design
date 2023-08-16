///================================ PCC Processing module HEADER file ===================================================================================///
///=========================================================================================================================================///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design tool (PDT code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1б 2022-2023 years   *///
///===================================================================================================================================///

/// Attached user-defined C++ libraries:
// Here a Processing_Functions.h library of C++ functions for advanced random and non-random generations of special chains of PCC elements
///-----------------------------------------------------

///* ========================================================= PCC PROCESSING FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
// special_cells_design - vector of vectors containing :: (1) snodes_sequence, (2) sedges_sequence, (3) sfaces_sequence, and (4) svolumes_sequence
/*! processing.ini
 * @param P_type
 * @param Configuration_State [not mandatory - only is there is some "non-zero" initial state]
 * @return cells_design
 */

CellsDesign PCC_Processing(std::vector<std::vector<int>> &Configuration_State,std::vector<std::vector<int>> &Configuration_cState);