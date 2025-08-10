///================================ PCC Subcomplex module ===================================================================================///
///=========================================================================================================================================///

#ifndef PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H
#define PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H

/*!
 * @brief Create a vector of PCC complexes with their special and induced labels taken from the initial PCC
 * @param configuration
 * @return std::vector<Subcomplex>
 */
std::vector<Subcomplex> PCC_Subcomplex(Config &configuration);
std::vector<Subcomplex> PCC_Subcomplex(Config &configuration, std::vector <unsigned int> &sub_sfaces_sequence);
std::vector<Subcomplex> PCC_Subcomplex(Config &configuration, std::vector <unsigned int> &sub_sfaces_sequence, std::vector <unsigned int> &sub_cfaces_sequence);

#endif //PCC_PROCESSING_DESIGN_PCC_SUBCOMPLEX_H
