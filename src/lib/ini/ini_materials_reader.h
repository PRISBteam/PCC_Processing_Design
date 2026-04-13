#ifndef PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H
#define PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H

/// Tailored Reader for materials.ini files in the ../CPD_material_database/ subdirectory of the project
// === # 1 # === //
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &Burgers_vector, double &gb_width, double &gb_inclusion1_adh_energy, double &lagbs_corrosion_current, double &hagbs_corrosion_current, double &sigma3_corrosion_current);

void material_irradiation_reader(std::string &Mid, double &austenite_interface_param, double &martensite_interface_param, double &ferrite_interface_param, double &austenite_martensite_interface_param, double &austenite_ferrite_interface_param, double &martensite_ferrite_interface_param);

// === # 2 # === //
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &gb_width, double &gb_inclusion1_adh_energy, std::string &I1id, std::string &inclusion1_type, double &inclusion_inclusion1_coh_energy, double &inclusion1_mass_density);

#endif //PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H
