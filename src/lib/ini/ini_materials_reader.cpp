/// Author: Dr Elijah Borodin (2024)
/// Manchester, UK
/// Library of specific functions related to the PCC Processing Design code for reading materials.ini database files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

/// Simple reader for *.ini files and a specific CPD code-related library for reading its particular *.ini files ( downloaded from https://github.com/pulzed/mINI )
#include "../ini/ini.h"

using namespace std; // standard namespace

extern std::string source_path;
extern std::string output_dir;
extern std::ofstream main_logfile_stream, kinetics_logfile_stream;

/// ================== # 1 # Materials - reading and output ==================
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &Burgers_vector, double &gb_width, double &gb_inclusion1_adh_energy, double &lagbs_corrosion_current, double &hagbs_corrosion_current, double &sigma3_corrosion_current) {

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "CPD_material_database/"s + Mid + ".ini"s);
    mINI::INIStructure materials_ini;
    file.read(materials_ini);

/// [structural]
// I. material type
    if (materials_ini.has("structural")) {
        auto &collection = materials_ini["structural"];
        if (collection.has("material_type")) {
            material_type = materials_ini.get("structural").get("material_type");
        }
        if (collection.has("Burgers_vector")) {
            Burgers_vector = stod(materials_ini.get("structural").get("Burgers_vector"));
        }
        // Grain Boundary width
        if (collection.has("gb_width")) {
            gb_width = stod(materials_ini.get("structural").get("gb_width"));
            gb_width *= std::pow(10, (-9));
        }
    }
    /// [corrosion]
    if (materials_ini.has("corrosion_current")) {
        auto& collection = materials_ini["corrosion_current"];
        if (collection.has("low_angle_gb")) {
            lagbs_corrosion_current = stod(materials_ini.get("corrosion_current").get("low_angle_gb"));
        }
        if (collection.has("high_angle_gb")) {
            hagbs_corrosion_current = stod(materials_ini.get("corrosion_current").get("high_angle_gb"));
        }
        if (collection.has("sigma3_gb")) {
            sigma3_corrosion_current = stod(materials_ini.get("corrosion_current").get("sigma3_gb"));
        }
    }

/// [thermodynamical]
// II. mass density
    if (materials_ini.has("thermodynamical")) {
        auto& collection = materials_ini["thermodynamical"];
        if (collection.has("mass_density")) {
            mass_density = stod(materials_ini.get("thermodynamical").get("mass_density"));
        }
// III. melting_point
        if (collection.has("melting_point")) {
            melting_point = stod(materials_ini.get("thermodynamical").get("melting_point"));
        }
// Grain boundary energy
        if (collection.has("cohesion_energy")) {
            cohesion_energy = stod(materials_ini.get("thermodynamical").get("cohesion_energy"));
        }
    }

/// [mechanical]
// IV. Young modulus
    if (materials_ini.has("mechanical")) {
        auto& collection = materials_ini["mechanical"];
        if (collection.has("Young_modulus")) {
            Young_modulus = stod(materials_ini.get("mechanical").get("Young_modulus"));
            Young_modulus = Young_modulus*std::pow(10,9);
        }
// V. Poisson ratio
        if (collection.has("Poisson_ratio")) {
            Poisson_ratio = stod(materials_ini.get("mechanical").get("Poisson_ratio"));
        }
// VI. yield strength
        if (collection.has("yield_strength")) {
            yield_strength = stod(materials_ini.get("mechanical").get("yield_strength"));
            yield_strength = yield_strength*std::pow(10,6);
        }
// VII. strength
        if (collection.has("strength")) {
            strength = stod(materials_ini.get("mechanical").get("strength"));
            strength = strength*std::pow(10,6);
        }
// VIII. fracture_toughness
        if (collection.has("fracture_toughness")) {
            fracture_toughness = stod(materials_ini.get("mechanical").get("fracture_toughness"));
        }
    }

/// [inclusions]
// IX. gb-inclusion1 adhesion energy
    if (materials_ini.has("inclusions")) {
        auto& collection = materials_ini["inclusions"];
        if (collection.has("gb_incl1_adhesion_energy")) {
            gb_inclusion1_adh_energy = stod(materials_ini.get("inclusions").get("gb_incl1_adhesion_energy"));
        }
    }

/// Console Output
    cout << endl << ".............................*  NEW MATERIAL  *.................................. " << endl << endl;
    cout << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    cout << "Material ID:.................... " << Mid << endl;
    cout << "Material type:.................. " << material_type << endl;
    if (mass_density > std::pow(10,-100)) cout << "Density [kg/m^3] .................... " << mass_density << endl;
    if (melting_point > std::pow(10,-100)) cout << "Melting point [K] .................... " << melting_point << endl;
    if (gb_width > std::pow(10,-100)) cout << "GB width [nm] .................. " << gb_width*pow(10,9) << endl;
    if (cohesion_energy > std::pow(10,-100)) cout << "Grain boundary energy density [J/m^2] .................... " << cohesion_energy << endl;
    if (Young_modulus > std::pow(10,-100)) cout << "Young modulus [GPa] .................... " << Young_modulus/pow(10,9) << endl;
    if (Poisson_ratio > std::pow(10,-100)) cout << "Poisson ratio  .................... " << Poisson_ratio << endl;
    if (yield_strength > std::pow(10,-100)) cout << "Yield strength [MPa]  .................... " << yield_strength/pow(10,6) << endl;
    if (strength > std::pow(10,-100)) cout << "Strength [MPa] .................... " << strength/pow(10,6) << endl;
    if (fracture_toughness > std::pow(10,-100)) cout << "Fracture toughness [MPa*sqrt(metre)] .................... " << fracture_toughness << endl;
    if (lagbs_corrosion_current > std::pow(10,-100)) cout << "LAGBs corrosion current [ ] .................... " << lagbs_corrosion_current << endl;
    if (hagbs_corrosion_current > std::pow(10,-100)) cout << "HAGBs corrosion current [ ] .................... " << hagbs_corrosion_current << endl;
    if (sigma3_corrosion_current > std::pow(10,-100)) cout << "Sigma3 corrosion current [ ] .................... " << sigma3_corrosion_current << endl;

/// Output into .log file
    main_logfile_stream << endl << ".............................*  NEW MATERIAL  *.................................. " << endl << endl;
    main_logfile_stream << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    main_logfile_stream << "Material ID:.................... " << Mid << endl;
    main_logfile_stream << "Material type:.................. " << material_type << endl;
    if (mass_density > std::pow(10,-100)) main_logfile_stream << "Density [kg/m^3] .................... " << mass_density << endl;
    if (melting_point > std::pow(10,-100)) main_logfile_stream << "Melting point [K] .................... " << melting_point << endl;
    if (gb_width > std::pow(10,-100)) main_logfile_stream << "GB width [nm] ................. " << gb_width*pow(10,9) << endl;
    if (cohesion_energy > std::pow(10,-100)) main_logfile_stream << "Grain boundary energy density [J/m^2] .................... " << cohesion_energy << endl;
    if (Young_modulus > std::pow(10,-100)) main_logfile_stream << "Young modulus [GPa] .................... " << Young_modulus/pow(10,9) << endl;
    if (Poisson_ratio > std::pow(10,-100)) main_logfile_stream << "Poisson ratio  .................... " << Poisson_ratio << endl;
    if (yield_strength > std::pow(10,-100)) main_logfile_stream << "Yield strength [MPa]  .................... " << yield_strength/pow(10,6) << endl;
    if (strength > std::pow(10,-100)) main_logfile_stream << "Strength [MPa] .................... " << strength/pow(10,6) << endl;
    if (fracture_toughness > std::pow(10,-100)) main_logfile_stream << "Fracture toughness [MPa*sqrt(metre)] .................... " << fracture_toughness << endl;
    if (lagbs_corrosion_current > std::pow(10,-100)) main_logfile_stream << "LAGBs corrosion current [ ] .................... " << lagbs_corrosion_current << endl;
    if (hagbs_corrosion_current > std::pow(10,-100)) main_logfile_stream << "HAGBs corrosion current [ ] .................... " << hagbs_corrosion_current << endl;
    if (sigma3_corrosion_current > std::pow(10,-100)) main_logfile_stream << "Sigma3 corrosion current [ ] .................... " << sigma3_corrosion_current << endl;

} // END of material_database_reader() function

void material_irradiation_reader(std::string &Mid, double &austenite_interface_param, double &martensite_interface_param, double &ferrite_interface_param, double &austenite_martensite_interface_param, double &austenite_ferrite_interface_param, double &martensite_ferrite_interface_param){
    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "CPD_material_database/"s + Mid + ".ini"s);
    mINI::INIStructure materials_ini;
    file.read(materials_ini);

    /// irradiation
    if (materials_ini.has("irradiation")) {
        auto& collection = materials_ini["irradiation"];
        if (collection.has("austenite_interface_param")) {
            austenite_interface_param = stod(materials_ini.get("irradiation").get("austenite_interface_param"));
        }

        if (collection.has("martensite_interface_param")) {
            martensite_interface_param = stod(materials_ini.get("irradiation").get("martensite_interface_param"));
        }
        if (collection.has("austenite_martensite_interface_param")) {
            austenite_martensite_interface_param = stod(materials_ini.get("irradiation").get("austenite_martensite_interface_param"));
        }
        if (collection.has("ferrite_interface_param")) {
            ferrite_interface_param = stod(materials_ini.get("irradiation").get("ferrite_interface_param"));
        }
        if (collection.has("austenite_ferrite_interface_param")) {
            austenite_ferrite_interface_param = stod(materials_ini.get("irradiation").get("austenite_ferrite_interface_param"));
        }
        if (collection.has("martensite_ferrite_interface_param")) {
            martensite_ferrite_interface_param = stod(materials_ini.get("irradiation").get("martensite_ferrite_interface_param"));
        }
    }
    /// Console Output
    cout << endl << ".............................*  Irradiation Material Characteristics  *.................................. " << endl << endl;
    cout << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    if (austenite_interface_param > std::pow(10,-100))
    cout << "Austenite interface:..................... " << austenite_interface_param << endl;
    if (martensite_interface_param > std::pow(10,-100))
    cout << "Martensite interface:.................... " << martensite_interface_param << endl;
    if (ferrite_interface_param > std::pow(10,-100))
    cout << "Ferrite interface:....................... " << ferrite_interface_param << endl;
    if (austenite_martensite_interface_param > std::pow(10,-100))
    cout << "Austenite/Martensite interface:.......... " << austenite_martensite_interface_param << endl;
    if (austenite_ferrite_interface_param > std::pow(10,-100))
    cout << "Austenite/Ferrite interface:............. " << austenite_ferrite_interface_param << endl;
    if (martensite_ferrite_interface_param > std::pow(10,-100))
    cout << "Martensite/Ferrite interface:............ " << martensite_ferrite_interface_param << endl;

    /// Console Output
    kinetics_logfile_stream << endl << ".............................*  Irradiation Material Characteristics  *.................................. " << endl << endl;
    kinetics_logfile_stream << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    if (austenite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Austenite interface:..................... " << austenite_interface_param << endl;
    if (martensite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Martensite interface:.................... " << martensite_interface_param << endl;
    if (ferrite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Ferrite interface:....................... " << ferrite_interface_param << endl;
    if (austenite_martensite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Austenite/Martensite interface:.......... " << austenite_martensite_interface_param << endl;
    if (austenite_ferrite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Austenite/Ferrite interface:............. " << austenite_ferrite_interface_param << endl;
    if (martensite_ferrite_interface_param > std::pow(10,-100))
        kinetics_logfile_stream << "Martensite/Ferrite interface:............ " << martensite_ferrite_interface_param << endl;

} // END of material_irradiation_reader()

/// ================== # 2 # Materials - reading and output ==================
void material_database_reader(std::string &Mid, std::string &material_type, double &mass_density, double &melting_point, double &cohesion_energy, double &Young_modulus, double &Poisson_ratio, double &yield_strength, double &strength, double &fracture_toughness, double &gb_width, double &gb_inclusion1_adh_energy, std::string &I1id, std::string &inclusion1_type, double &inclusion1_coh_energy, double &inclusion1_mass_density) {

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "CPD_material_database/"s + Mid + ".ini"s);
    mINI::INIStructure materials_ini;
    file.read(materials_ini);

/// [structural]
// I. material type
    if (materials_ini.has("structural")) {
        auto& collection = materials_ini["structural"];
        if (collection.has("material_type")) {
            material_type = materials_ini.get("structural").get("material_type");
        }

    // Grain Boundary width
        if (collection.has("gb_width")) {
            gb_width = stod(materials_ini.get("structural").get("gb_width"));
            gb_width *= std::pow(10,(-9));
        }
    }

/// [thermodynamical]
// II. mass density
    if (materials_ini.has("thermodynamical")) {
        auto& collection = materials_ini["thermodynamical"];
        if (collection.has("mass_density")) {
            mass_density = stod(materials_ini.get("thermodynamical").get("mass_density"));
        }
// III. melting_point
        if (collection.has("melting_point")) {
            melting_point = stod(materials_ini.get("thermodynamical").get("melting_point"));
        }
// Grain boundary energy
        if (collection.has("gb_cohesion_energy")) {
            cohesion_energy = stod(materials_ini.get("thermodynamical").get("gb_cohesion_energy"));
        }
    }

/// [mechanical]
// IV. Young modulus
    if (materials_ini.has("mechanical")) {
        auto& collection = materials_ini["mechanical"];
        if (collection.has("Young_modulus")) {
            Young_modulus = stod(materials_ini.get("mechanical").get("Young_modulus"));
            Young_modulus = Young_modulus*std::pow(10,9);
        }
// V. Poisson ratio
        if (collection.has("Poisson_ratio")) {
            Poisson_ratio = stod(materials_ini.get("mechanical").get("Poisson_ratio"));
        }
// VI. yield strength
        if (collection.has("yield_strength")) {
            yield_strength = stod(materials_ini.get("mechanical").get("yield_strength"));
            yield_strength = yield_strength*std::pow(10,6);
        }
// VII. strength
        if (collection.has("strength")) {
            strength = stod(materials_ini.get("mechanical").get("strength"));
            strength = strength*std::pow(10,6);
        }
// VIII. fracture_toughness
        if (collection.has("fracture_toughness")) {
            fracture_toughness = stod(materials_ini.get("mechanical").get("fracture_toughness"));
        }
    }

/// [inclusions]
// IX. gb-inclusion1 adhesion energy
    if (materials_ini.has("inclusions")) {
        auto& collection = materials_ini["inclusions"];
        if (collection.has("gb_adhesion_energy")) {
            gb_inclusion1_adh_energy = stod(materials_ini.get("inclusions").get("gb_adhesion_energy"));
        }
    }

    /// ============================== Inclusion materials reader ==========================
    mINI::INIFile file2(source_path + "CPD_material_database/"s + I1id + ".ini"s);
    mINI::INIStructure inclusions_ini;
    file2.read(inclusions_ini);

/// [structural]
// I. Inclusion material type
    if (inclusions_ini.has("structural")) {
        auto& collection = inclusions_ini["structural"];
        if (inclusions_ini.has("material_type")) {
            inclusion1_type = inclusions_ini.get("structural").get("material_type");
        } }

    // Inclusion mass density
    if (inclusions_ini.has("thermodynamical")) {
        auto& collection = inclusions_ini["thermodynamical"];
        if (collection.has("mass_density")) {
            inclusion1_mass_density = stod(inclusions_ini.get("thermodynamical").get("mass_density"));
        } }

    // Inclusion cohesion energy density
    if (inclusions_ini.has("thermodynamical")) {
        auto& collection = inclusions_ini["thermodynamical"];
        if (collection.has("aggl_cohesion")) {
            inclusion1_coh_energy = stod(inclusions_ini.get("thermodynamical").get("aggl_cohesion"));
        } }


/// Console Output
    cout << endl << ".............................*  MATRIX MATERIAL  *.................................. " << endl << endl;
    cout << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    cout << "Material ID: " << Mid << endl;
    cout << "Material type: " << material_type << endl;
    if (mass_density > std::pow(10,-100)) cout << "Density [kg/m^3] .................... " << mass_density << endl;
    if (melting_point > std::pow(10,-100)) cout << "Melting point [K] .................... " << melting_point << endl;
    if (gb_width > std::pow(10,-100)) cout << "GB width [nm] .................... " << gb_width*pow(10,9) << endl;
    if (cohesion_energy > std::pow(10,-100)) cout << "Matrix-matrix interface energy density [J/m^2] .................... " << cohesion_energy << endl;
    if (inclusion1_coh_energy > std::pow(10,-100)) cout << "Matrix-inclusion interface energy density [J/m^2] .................... " << gb_inclusion1_adh_energy << endl;
    if (Young_modulus > std::pow(10,-100)) cout << "Young modulus [GPa] .................... " << Young_modulus/pow(10,9) << endl;
    if (Poisson_ratio > std::pow(10,-100)) cout << "Poisson ratio  .................... " << Poisson_ratio << endl;
    if (yield_strength > std::pow(10,-100)) cout << "Yield strength [MPa]  .................... " << yield_strength/pow(10,6) << endl;
    if (strength > std::pow(10,-100)) cout << "Strength [MPa] .................... " << strength/pow(10,6) << endl;
    if (fracture_toughness > std::pow(10,-100)) cout << "Fracture toughness [MPa*sqrt(metre)] .................... " << fracture_toughness << endl;

    cout << endl << ".............................*  INCLUSION MATERIAL  *.................................. " << endl;
    cout << "Inclusion material source: " << source_path + "CPD_material_database/"s + I1id + ".ini"s << endl;
    cout << "Inclusion ID: " << I1id << endl;
    cout << "Inclusion type: " << inclusion1_type << endl;
    if (inclusion1_mass_density > std::pow(10,-100)) cout << "Inclusion density [kg/m^3] .............. " << inclusion1_mass_density << endl;
    if (inclusion1_coh_energy > std::pow(10,-100)) cout << "Agglomeration cohesion energy density [J/m^2] .............. " << inclusion1_coh_energy << endl;

/// Output into .log file
    main_logfile_stream << endl << ".............................*  MATRIX MATERIAL  *.................................. " << endl << endl;
    main_logfile_stream << "Material source: " << source_path + "CPD_material_database/"s + Mid + ".ini"s << endl;
    main_logfile_stream << "Material ID: " << Mid << endl;
    main_logfile_stream << "Material type: " << material_type << endl;
    if (mass_density > std::pow(10,-100)) main_logfile_stream << "Density [kg/m^3] .................... " << mass_density << endl;
    if (melting_point > std::pow(10,-100)) main_logfile_stream << "Melting point [K] .................... " << melting_point << endl;
    if (gb_width > std::pow(10,-100)) main_logfile_stream << "GB width [nm] .................... " << gb_width*pow(10,9) << endl;
    if (cohesion_energy > std::pow(10,-100)) main_logfile_stream << "Grain boundary energy density [J/m^2] .................... " << cohesion_energy << endl;
    if (Young_modulus > std::pow(10,-100)) main_logfile_stream << "Young modulus [GPa] .................... " << Young_modulus/pow(10,9) << endl;
    if (Poisson_ratio > std::pow(10,-100)) main_logfile_stream << "Poisson ratio  .................... " << Poisson_ratio << endl;
    if (yield_strength > std::pow(10,-100)) main_logfile_stream << "Yield strength [MPa]  .................... " << yield_strength/pow(10,6) << endl;
    if (strength > std::pow(10,-100)) main_logfile_stream << "Strength [MPa] .................... " << strength/pow(10,6) << endl;
    if (fracture_toughness > std::pow(10,-100)) main_logfile_stream << "Fracture toughness [MPa*sqrt(metre)] .................... " << fracture_toughness << endl;
    main_logfile_stream << endl << ".............................*  INCLUSION MATERIAL  *.................................. " << endl;
    main_logfile_stream << "Inclusion material source: " << source_path + "CPD_material_database/"s + I1id + ".ini"s << endl;
    main_logfile_stream << "Inclusion ID: " << I1id << endl;
    main_logfile_stream << "Inclusion type: " << inclusion1_type << endl;
    if (inclusion1_mass_density > std::pow(10,-100)) main_logfile_stream << "Inclusion density [kg/m^3] .............. " << inclusion1_mass_density << endl;
    if (inclusion1_coh_energy > std::pow(10,-100)) main_logfile_stream << "Agglomeration cohesion energy density [J/m^2] .............. " << inclusion1_coh_energy << endl;

//    Out_logfile_stream.close();
} // END of material_database_reader() function