#ifndef PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H
#define PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H

// creation plane cut as a PCC Subcomplex class
std::set<unsigned int> PCC_Subcomplex_plane_cut_grains (std::vector<double> &plane_orientation);
// overloaded version
std::set<unsigned int> PCC_Subcomplex_plane_cut_grains (double a_coeff, double b_coeff, double c_coeff, double D_coeff);

/// Creation of a vector of k-order grain neighbours as a PCC Subcomplex class
std::set<unsigned int> PCC_Subcomplex_k_order_grain_neighbours(unsigned int grain_number, int k_neighbours_order);
/// Cutting PART defined in the config/subcomplex.ini file of the initial plane cut
Subcomplex Get_half_plane(Subcomplex &plane_subcomplex, double crack_length, int macrocrack_grow_direction);

#endif //PCC_PROCESSING_DESIGN_SUBCOMPLEX_CROSS_SECTION_H
