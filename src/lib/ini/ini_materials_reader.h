#ifndef PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H
#define PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H

/// Tailored Reader for materials.ini files in the ../CPD_material_database/ subdirectory of the project
// === # 1 # === //
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &gb_width, double &gb_inclusion1_adh_energy);

// === # 2 # === //
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &gb_width, double &gb_inclusion1_adh_energy, std::string &I1id, std::string &inclusion1_type, double &inclusion_inclusion1_coh_energy, double &inclusion1_mass_density);

#endif //PCC_PROCESSING_DESIGN_INI_MATERIALS_READER_H
