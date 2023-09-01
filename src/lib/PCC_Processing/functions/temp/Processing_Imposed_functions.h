///================================ A part of the PCC Processing module HEADER file =============================================================///
///=================================================================================================================================///
/** The library contains functions for generation of new identifiers for m-Cells, m={0,1,2} based on the already assigned types    **/
/**  of the k-Cells with k > m, in the PCC Processing module. It created an "imposed" structure of m-Cells in a PCC.              **/
///==============================================================================================================================///

vector<int> Edge_types_byFaces(std::vector<unsigned int> const &CellNumbs, std::vector<unsigned int> &special_face_sequence, std::vector<double> &j_fractions, std::vector<double> &d_fractions);