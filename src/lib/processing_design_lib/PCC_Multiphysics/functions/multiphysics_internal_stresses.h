#ifndef PCC_PROCESSING_DESIGN_MULTIPHYSICS_INTERNAL_STRESS_H
#define PCC_PROCESSING_DESIGN_MULTIPHYSICS_INTERNAL_STRESS_H

/// Setting internal energies, such as local elastic energies of stress concentrators
/*!
 * @brief Calculate stress field and energy density in a vicinity of a macrocrack.
 * @param new_crack
 * @param matrix_material
 * @param external_stress
 * @param sample_dimensions
 * @return local energy density of faces (2-cells)
 */
std::vector<double> Multiphysics_crack_stress_field(Macrocrack &new_crack, Material &matrix_material, Eigen::MatrixXd &external_stress, std::tuple<double, double, double> &sample_dimensions);

#endif //PCC_PROCESSING_DESIGN_MULTIPHYSICS_INTERNAL_STRESS_H
