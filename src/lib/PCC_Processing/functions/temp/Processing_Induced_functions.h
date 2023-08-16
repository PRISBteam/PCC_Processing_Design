///================================ A part of the PCC Processing module HEADER file =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new "secondary" identifications of k-Cells   **/
/**  in the PCC Processing module. It makes them "fractured" and can go alongside their special or ordinary primal types.         **/
///==============================================================================================================================///

/// #1# Kinematic function for multiple cracking
std::vector <unsigned int> PCC_Kinematic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<vector<int>> &Configuration_cState, std::vector<vector<double>> const &max_cfractions_vectors);

/// #2# Kinetic function for multiple cracking
// std::vector <unsigned int> PCC_Kinetic_cracking(std::vector <double> &face_elastic_energies, std::vector<unsigned int> &frac_sfaces_sequence, macrocrack &large_crack, Eigen::SparseMatrix<double> const& AFS, Eigen::SparseMatrix<double> const& FES);