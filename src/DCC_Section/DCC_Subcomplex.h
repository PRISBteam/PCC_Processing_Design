/// Attached user defined C++ libraries:
///-------------------------------------
#include "Planecut_Functions.h"
///-------------------------------------

std::vector <unsigned int> DCC_Subcomplex(double a_coeff, double b_coeff, double c_coeff, double D_coeff, std::vector<unsigned int> const &s_faces_sequence) {
std::vector <unsigned int>  sub_faces_sequence, s_sub_faces_sequence;

sub_faces_sequence = DCC_Plane_cut (a_coeff, b_coeff, c_coeff, D_coeff);
//s_faces_sequence

return s_sub_faces_sequence;
}

/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence