#ifndef PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H
#define PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H

/// Creation a subcomplex of a PCC
Subcomplex PCC_Subcomplex(Subcomplex &new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> &sub_faces_sequence, std::vector<unsigned int> const &c_faces_sequence, double a_coeff, double b_coeff, double c_coeff, double D_coeff);

#endif //PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H
