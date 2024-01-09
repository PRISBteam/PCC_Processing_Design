
#ifndef PCC_PROCESSING_DESIGN_MULTIPHYSICS_FUNCTIONS_H
#define PCC_PROCESSING_DESIGN_MULTIPHYSICS_FUNCTIONS_H

/// Set stress field for a Micro-crack (class Macrocrack)
void crack_modes_stress_field(std::vector <double> &face_energies, int crack_mode, double new_crack_length, double external_vonMizes_stress, double Puasson_coeff = 0.3, double Young_modulus = 215.0*pow(10.0,9));

#endif //PCC_PROCESSING_DESIGN_MULTIPHYSICS_FUNCTIONS_H
