#ifndef PCC_PROCESSING_DESIGN_SUBCOMPLEX_PLANECUT_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_SUBCOMPLEX_PLANECUT_FUNCTIONS_H

/// Creation plane cut as a PCC Subcomplex class
std::vector<unsigned int> PCC_Plane_cut (double a_coeff, double b_coeff, double c_coeff, double D_coeff);

/// Cutting PART defined in the config/subcomplex.ini file of the initial plane cut
Subcomplex Get_half_plane(Subcomplex new_sub, double crack_length, std::vector<unsigned int> const &half_sub_sfaces_sequence);

#endif //PCC_PROCESSING_DESIGN_SUBCOMPLEX_PLANECUT_FUNCTIONS_H
