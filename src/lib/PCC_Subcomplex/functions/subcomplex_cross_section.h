#ifndef PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H
#define PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H

/// Creation plane cut as a PCC Subcomplex class
std::set<unsigned int> PCC_Plane_cut_grains (std::vector<double> &plane_orientation);

// overloaded version
std::set<unsigned int> PCC_Plane_cut_grains (double a_coeff, double b_coeff, double c_coeff, double D_coeff);

/// Cutting PART defined in the config/subcomplex.ini file of the initial plane cut
Subcomplex Get_half_plane(Subcomplex &new_sub, double crack_length);

#endif //PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H
