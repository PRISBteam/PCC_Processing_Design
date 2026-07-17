/// Author: Dr Elijah Borodin (2023)
/// Manchester, UK
/// Library of specific functions related to the PCC Processing Design code for reading its *.ini files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2024)
#include "../../../src/lib/external/Eigen-5.0/Dense"
#include "../../../src/lib/external/Eigen-5.0/SparseCore"

#include "../../../src/lib/processing_design_lib/PCC_Objects.h"

/// Simple reader for *.ini files and a specific CPD code-related library for reading its particular *.ini files ( downloaded from https://github.com/pulzed/mINI )
#include "../ini/ini.h"
#include "../ini/ini_readers.h"

using namespace std; // standard namespace

extern std::string source_path;
extern std::string output_dir;
extern std::ofstream main_logfile_stream, subcomplex_logfile_stream, multiphysics_logfile_stream, processing_logfile_stream, kinetics_logfile_stream, characterisation_logfile_stream, design_logfile_stream, writer_logfile_stream;

/// ================== # 1 # Initial configuration - reading and output ==================
//std::vector<int> config_reader_main(std::string &pcc_source_dir, std::string &output_dir, std::string &cell_complex_standard, std::string &main_type) {
/*!
 * @details Read input parameters from the project file config/main.ini necessary for the code execution. Print the read values to a screen and main_logfile_stream --> cpdlog_main.log file.
 * @param main_ini_data
 * @return std::vector<int> res
 */
//Config &configuration
//std::vector<int> config_reader_main(main_config &main_ini_data) {
std::vector<int> config_reader_main(Config &configuration) {
    std::vector<int> res(8,0);

/// [0] - > dim, [1] -> isSubcomplex, [2] -> isMultiphysics, [3] -> isProcessing, [4] -> isCharacterisation, [5] -> isDesign, [6] -> isWriter, [7] -> isKinetics
        bool isSubcomplexON = 0, isProcessingON = 0, isCharacterisationON = 0, isKineticsON = 0, isMultiphysicsON = 0, isDesignON = 0, isWriterON = 0;
        std::string isSubcomplex, isProcessing, isKinetics, isMultiphysics, isCharacterisation, isDesign, isWriter;

        // ini files reader - external (MIT license) library
        mINI::INIFile file(source_path + "main.ini"s);
        mINI::INIStructure main_ini;
        file.read(main_ini);

// I
        if (main_ini.has("execution_mode")) {
            auto &collection = main_ini["execution_mode"];
            if (collection.has("mode")) {
                std::string main_type_str;
                if (main_ini.get("execution_mode").get("mode") == "LIST" ||
                    main_ini.get("execution_mode").get("mode") == "TUTORIAL" ||
                    main_ini.get("execution_mode").get("mode") == "PERFORMANCE_TEST" ||
                    main_ini.get("execution_mode").get("mode") == "TASK") {
                    main_type_str = main_ini.get("execution_mode").get("mode");
                    configuration.Set_main_type(main_type_str);
                } else
                    throw std::invalid_argument(
                            "ERROR in ../src/ini/ini_readers.cpp: WRONG TYPE OF THE 'execution_mode' IN ../config/main.ini FILE; Please change the mode to one of the allowed: 'LIST', 'TUTORIAL', 'PERFORMANCE_TEST' or 'TASK' "s);
            }
        }
// II
        std::string problem_dimension;
        if (main_ini.has("general")) {
            auto &collection = main_ini["general"];
            if (collection.has("PCC_dimension")) {
                if (stoi(main_ini.get("general").get("PCC_dimension")) == 1 ||
                    stoi(main_ini.get("general").get("PCC_dimension")) == 2 ||
                    stoi(main_ini.get("general").get("PCC_dimension")) == 3) {
                    problem_dimension = main_ini.get("general").get("PCC_dimension");
                } else
                    throw std::invalid_argument(
                            "ERROR in ../src/ini/ini_readers.cpp: WRONG TYPE OF DIMENSION 'PCC_dimension' IN ../config/main.ini FILE; Please change the [general] PCC_dimension parameter to one of the allowed: 1, 2 or 3 "s);
            }
        }
        res.at(0) = stoi(problem_dimension); // res[0]

        // III Modules ON/OFF
        if (main_ini.has("modules")) {
            auto &collection = main_ini["modules"];
            if (collection.has("PCC_Subcomplex"))
                isSubcomplex = main_ini.get("modules").get("PCC_Subcomplex");

            if (collection.has("PCC_Multiphysics"))
                isMultiphysics = main_ini.get("modules").get("PCC_Multiphysics");

            if (collection.has("PCC_Processing"))
                isProcessing = main_ini.get("modules").get("PCC_Processing");

            if (collection.has("PCC_Kinetics"))
                isKinetics = main_ini.get("modules").get("PCC_Kinetics");

            if (collection.has("PCC_Characterisation"))
                isCharacterisation = main_ini.get("modules").get("PCC_Characterisation");

            if (collection.has("PCC_Design"))
                isDesign = main_ini.get("modules").get("PCC_Design");

            if (collection.has("PCC_Writer"))
                isWriter = main_ini.get("modules").get("PCC_Writer");

        } // end if(main_ini.has("modules"))

/// Forming the output RES vector
        if (isSubcomplex == "ON") {
            isSubcomplexON = 1;
            res.at(1) = 1;
        }
        else res.at(1) = 0; // res[1] - Section -> ConfigVector.at(1) in main.cpp
        if (isMultiphysics == "ON") {
            isMultiphysicsON = 1;
            res.at(2) = 1;
        }
        else res.at(2) = 0; // res[2] - Multiphysics -> ConfigVector.at(2) in main.cpp
        if (isProcessing == "ON") {
            isProcessingON = 1;
            res.at(3) = 1;
        }
        else res.at(3) = 0; // res[3] - Processing -> ConfigVector.at(3) in main.cpp
        if (isKinetics == "ON") {
            isKineticsON = 1;
            res.at(7) = 1;
        }
        else res.at(7) = 0; // res[7] - Kinetics -> ConfigVector.at(7) in main.cpp
        if (isCharacterisation == "ON") {
            isCharacterisationON = 1;
            res.at(4) = 1;
        }
        else res.at(4) = 0; // res[4] - Characterisation -> ConfigVector.at(4) in main.cpp
        if (isDesign == "ON") {
            isDesignON = 1;
            res.at(5) = 1;
        }
        else res.at(5) = 0; // res[5] - Design -> ConfigVector.at(5) in main.cpp
        if (isWriter == "ON") {
            isWriterON = 1;
            res.at(6) = 1;
        }
        else res.at(6) = 0; // res[6] - Writer -> ConfigVector.at(6) in main.cpp

        // additional parameters
        if (main_ini.has("general")) {
            auto &collection = main_ini["general"];
            std::string pcc_source_directory;
            if (collection.has("pcc_source_dir"))
                pcc_source_directory = main_ini.get("general").get("pcc_source_dir");
            configuration.Set_pcc_source_dir(pcc_source_directory);
        }

        if (main_ini.has("general")) {
            auto &collection = main_ini["general"];
            std::string pcc_standard_id;
            if (collection.has("pcc_standard"))
                pcc_standard_id = main_ini.get("general").get("pcc_standard");
            configuration.Set_pcc_standard_id(pcc_standard_id);
        } // like 'pcc1s' - the first computational standard (2024)

        if (main_ini.has("general")) {
            auto &collection = main_ini["general"];
            if (collection.has("output_dir"))
                output_dir = main_ini.get("general").get("output_dir");
            configuration.Set_output_dir(output_dir);
        }

/// Output to the screen/console
    if (configuration.main_reader_switch) {
        main_logfile_stream.open(output_dir + "cpd_main.log"s, ios::trunc); // the main_logfile_stream.log stream will be closed at the end of the main function

        cout << "The problem dimension that is the maximum value k_max of k-cells in the PCC\t\t|\t\t"s << "dim = "
             << res.at(0) << endl;
        cout << "Execution mode:\t"s << "\t" << configuration.Get_main_config().main_type << endl;
        cout << "Output directory:\t"s << "\t" << output_dir << endl;
        cout << "PCC source directory:\t"s << configuration.Get_main_config().pcc_source_dir << endl;
        cout << "PCC standard ID:\t\t"s << configuration.Get_main_config().pcc_standard << endl;
        cout << endl;

        if (isSubcomplexON == 1) cout << "ON    | PCC_Subcomplex"s << endl;
        else cout << "OFF    | PCC_Subcomplex"s << endl;
        if (isMultiphysicsON == 1) cout << "ON    | PCC_Multiphysics"s << endl;
        else cout << "OFF    | PCC_Multiphysics"s << endl;
        if (isProcessingON == 1) cout << "ON    | PCC_Processing"s << endl;
        else cout << "OFF    | PCC_Processing"s << endl;
        if (isKineticsON == 1) cout << "ON    | PCC_Kinetics"s << endl;
        else cout << "OFF    | PCC_Kinetics"s << endl;
        if (isCharacterisationON == 1) cout << "ON    | PCC_Characterisation"s << endl;
        else cout << "OFF    | PCC_Characterisation"s << endl;
        if (isDesignON == 1) cout << "ON    | PCC_Design"s << endl;
        else cout << "OFF    | PCC_Design"s << endl;
        if (isWriterON == 1) cout << "ON    | PCC_Writer"s << endl;
        else cout << "OFF    | PCC_Writer"s << endl;
        cout << endl;

/// Output into main.log file
        main_logfile_stream << endl;
        main_logfile_stream << "The problem dimension that is the maximum value k_max of k-cells in the PCC:\t\t|\t\t"s
                            << "dim = " << res.at(0) << endl;
        main_logfile_stream << "Execution mode:\t"s << "\t" << configuration.Get_main_config().main_type << endl;
        main_logfile_stream << "Output directory:\t"s << "\t" << output_dir << endl;
        main_logfile_stream << "PCC source directory:\t"s << configuration.Get_main_config().pcc_source_dir << endl;
        main_logfile_stream << "PCC standard ID:\t\t"s << configuration.Get_main_config().pcc_standard << endl;

        main_logfile_stream << endl;
        if (isSubcomplexON == 1) main_logfile_stream << "ON    | PCC_Subcomplex"s << endl;
        else main_logfile_stream << "OFF    | PCC_Subcomplex"s << endl;
        if (isMultiphysicsON == 1) main_logfile_stream << "ON    | PCC_Multiphysics"s << endl;
        else main_logfile_stream << "OFF    | PCC_Multiphysics"s << endl;
        if (isProcessingON == 1) main_logfile_stream << "ON    | PCC_Processing"s << endl;
        else main_logfile_stream << "OFF    | PCC_Processing"s << endl;
        if (isKineticsON == 1) main_logfile_stream << "ON    | PCC_Kinetics"s << endl;
        else main_logfile_stream << "OFF    | PCC_Kinetics"s << endl;
        if (isCharacterisationON == 1) main_logfile_stream << "ON    | PCC_Characterisation"s << endl;
        else main_logfile_stream << "OFF    | PCC_Characterisation"s << endl;
        if (isDesignON == 1) main_logfile_stream << "ON    | PCC_Design"s << endl;
        else main_logfile_stream << "OFF    | PCC_Design"s << endl;
        if (isWriterON == 1) main_logfile_stream << "ON    | PCC_Writer"s << endl;
        else main_logfile_stream << "OFF    | PCC_Writer"s << endl;
        main_logfile_stream << endl;

        main_logfile_stream.close();
        configuration.main_reader_switch = false;
    } //  END if (configuration.main_reader_switch)

    return res;
} /// END of the 'config_reader_main()' function

/// ================== # 2 # Initial SUBCOMPLEX module configuration - reading and output ==================
/*!
 * @details Read input parameters from the project file config/subcomplex.ini necessary for the code execution. Print the read values to a screen and subcomplex_logfile_stream --> cpdlog_subcomplex.log file.
 * @param sctype
 * @param plane_orientation
 * @param cut_length
 * @param grain_neighbour_orders
 * @return void
 */
//void config_reader_subcomplex(std::string &sctype, std::vector<double> &plane_orientation, double &cut_length, unsigned int &grain_neighbour_orders, bool &is_log_file) {
void config_reader_subcomplex(Config &configuration) {

    std::string log_file_output;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "subcomplex.ini"s);
    mINI::INIStructure subcomplex_ini;
    file.read(subcomplex_ini);

    std::string subcomplex_mode;
    std::vector<double> subcomplex_plane;
    double cut_length;
    unsigned int grain_neighbour_orders;

//subcomplex type
    if (subcomplex_ini.has("subcomplex_type")) {
        auto& collection = subcomplex_ini["subcomplex_type"];
        if (collection.has("subPCC_type"))
        {
            subcomplex_mode = subcomplex_ini.get("subcomplex_type").get("subPCC_type");
            configuration.Set_subcomplex_mode(subcomplex_mode);
        } }

    //plane orientation
    double a_coeff = 0, b_coeff = 0, c_coeff = 0, D_coeff = 0;
    if (subcomplex_ini.has("plane_section")) {
        auto& collection = subcomplex_ini["plane_section"];
        if (collection.has("a_coeff")) {
            a_coeff = stod(subcomplex_ini.get("plane_section").get("a_coeff"));
            subcomplex_plane.push_back(a_coeff);
        }
        if (collection.has("b_coeff")) {
            b_coeff = stod(subcomplex_ini.get("plane_section").get("b_coeff"));
            subcomplex_plane.push_back(b_coeff);
        }
        if (collection.has("c_coeff")) {
            c_coeff = stod(subcomplex_ini.get("plane_section").get("c_coeff"));
            subcomplex_plane.push_back(c_coeff);
        }
        if (collection.has("D_coeff")) {
            D_coeff = stod(subcomplex_ini.get("plane_section").get("D_coeff"));
            subcomplex_plane.push_back(D_coeff);
        }
        configuration.Set_subcomplex_plane(subcomplex_plane);
    } // end of if (subcomplex_ini.has("plane_section"))

    //half-plane cut length
    if (subcomplex_ini.has("half_plane_section")) {
        auto& collection = subcomplex_ini["half_plane_section"];
        if (collection.has("half_plane_length"))
        {
            cut_length = stod(subcomplex_ini.get("half_plane_section").get("half_plane_length"));
            configuration.Set_subcomplex_cut_length(cut_length);
        } }

    // k_order_neighbours
    if (subcomplex_ini.has("k_order_neighbours")) {
        auto& collection = subcomplex_ini["k_order_neighbours"];
        if (collection.has("neighbours_order"))
        {
            grain_neighbour_orders = stoi(subcomplex_ini.get("k_order_neighbours").get("neighbours_order"));
            configuration.Set_subcomplex_grain_neighbour_orders(grain_neighbour_orders);
        } }

    // module output
        if (subcomplex_ini.has("module_output")) {
        auto& collection = subcomplex_ini["module_output"];
        if (collection.has("module_log_file"))
        {
            log_file_output = subcomplex_ini.get("module_output").get("module_log_file");
        } }
    if (log_file_output == "ON") configuration.Set_is_subcomplex_log_file(true);

/// Output to the screen/console
    if (configuration.subcomplex_reader_switch) {
        cout << "The Subcomplex module type and initial parameters:\t\t" << endl << endl;
        cout << "Subcomplex type:\t"s << configuration.Get_subcomplex_mode()<< endl;
        if (configuration.Get_subcomplex_mode() == "H"s) {
            cout << "Half-plane length:\t"s << configuration.Get_subcomplex_cut_length() << endl << endl;
        }
        if (configuration.Get_subcomplex_mode() == "P"s ||
            configuration.Get_subcomplex_mode() == "H"s) {
            cout << "Plane orientation:\ta_coeff*X + b_coeff*Y + c_coeff*Z = D"s << endl << "Plane normal vector\t"s
                 << "\ta =\t" << configuration.Get_subcomplex_plane().at(0) << "\tb =\t"
                 << configuration.Get_subcomplex_plane().at(1) << "\tc =\t"
                 << configuration.Get_subcomplex_plane().at(2) << "\t\tPlane position D =\t"s
                 << configuration.Get_subcomplex_plane().at(3) << endl;
        } else if (configuration.Get_subcomplex_mode() == "N"s) {
            cout << "Grain k-neighbours order:\t"s << configuration.Get_subcomplex_grain_neighbour_orders() << endl << endl;
        }
        cout << "Subcomplex module cpdlog_subcomplex.log file output:\t"s << log_file_output << endl;

        if (configuration.Get_is_subcomplex_log_file()) {
            subcomplex_logfile_stream << "The Subcomplex module type and initial parameters:\t\t" << endl << endl;
            subcomplex_logfile_stream << "Subcomplex type:\t"s << configuration.Get_subcomplex_mode() << endl;
            if (configuration.Get_subcomplex_mode() == "H"s) {
                subcomplex_logfile_stream << "Half-plane length:\t"s << configuration.Get_subcomplex_cut_length() << endl << endl;
            }
            if (configuration.Get_subcomplex_mode() == "P"s ||
                configuration.Get_subcomplex_mode() == "H"s) {
                subcomplex_logfile_stream << "Plane orientation:\ta_coeff*X + b_coeff*Y + c_coeff*Z = D"s << endl
                                          << "Plane normal vector:\t"s << " a: "
                                          << configuration.Get_subcomplex_plane().at(0) << " b: "
                                          << configuration.Get_subcomplex_plane().at(1) << " c: "
                                          << configuration.Get_subcomplex_plane().at(2)
                                          << "\tPlane position\t" << " D: "
                                          << configuration.Get_subcomplex_plane().at(3) << endl;
            } else if (configuration.Get_subcomplex_mode() == "N"s) {
                subcomplex_logfile_stream << "Grain k-neighbours order:\t"s
                                          << configuration.Get_subcomplex_grain_neighbour_orders() << endl
                                          << endl;
            }
        }
        configuration.subcomplex_reader_switch = false;
    } // end if (configuration.subcomplex_reader_switch) {

    return;
} /// end of the 'config_reader_subcomplex() function


/// ================== # 3 # Initial MULTIPFYSICS module configuration - physical dimesions and all ==================
/*!
 * @details Read input parameters from the project file config/multiphysics.ini necessary for the code execution. Print the read values to a screen and multiphysics_logfile_stream --> cpdlog_processing.log file.
 * @param Mid_matrix
 * @param Mid_inclusion1
 * @param sample_dimensions
 * @param tau
 * @param ext_stress_tensor
 * @param macrocrack_ini
 * @return void
 */
//void config_reader_multiphysics(std::string &Mid_matrix, std::string &Mid_inclusion1, std::tuple<double, double, double> &sample_dimensions, double &tau, Eigen::MatrixXd &ext_stress_tensor, std::vector<double> &macrocrack_ini, bool &is_log_file) {
void config_reader_multiphysics(Config &configuration) {

    double lx_size = 0.0, ly_size = 0.0, lz_size = 0.0; // sample dimensions
    double sxx = 0.0, sxy = 0.0, sxz = 0.0, syx = 0.0, syy = 0.0, syz = 0.0, szx = 0.0, szy = 0.0, szz = 0.0; // external stress tensor components [homogeneous stress state]
    double new_ambient_temperature = 300.0;
    std::string log_file_output;

// ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "multiphysics.ini"s);
    mINI::INIStructure multiphysics_ini;
    file.read(multiphysics_ini);

// I
// Material ID for the CPD code Database
    std::string Mid_matrix, Mid_inclusion1;
    if (multiphysics_ini.has("material_id")) {
        auto &collection = multiphysics_ini["material_id"];
        if (collection.has("Mid_matrix")) {
            Mid_matrix = multiphysics_ini.get("material_id").get("Mid_matrix");
            configuration.Set_multiphysics_matrixMaterial_id(Mid_matrix);
        }
// Material ID for the CPD code Database
        if (collection.has("Mid_inclusion1")) {
            Mid_inclusion1 = multiphysics_ini.get("material_id").get("Mid_inclusion1");
            configuration.Set_multiphysics_inclusionMaterial_id(Mid_inclusion1);
        }
    }

// II
// sequences and designs output
    if (multiphysics_ini.has("sample_dimensions")) {
        auto &collection = multiphysics_ini["sample_dimensions"];
        if (collection.has("lx"))
            lx_size = stod(multiphysics_ini.get("sample_dimensions").get("lx"));
        if (collection.has("ly"))
            ly_size = stod(multiphysics_ini.get("sample_dimensions").get("ly"));
        if (collection.has("lz"))
            lz_size = stod(multiphysics_ini.get("sample_dimensions").get("lz"));
    } // end of  if (multiphysics_ini.has("physical_dimensions"))
    std::tuple<double, double, double> new_sample_dimensions = make_tuple(lx_size,ly_size,lz_size);
    configuration.Set_multiphysics_sample_dimensions(new_sample_dimensions);

// III
    if (multiphysics_ini.has("time_scale")) {
        auto &collection = multiphysics_ini["time_scale"];
        double tau_parameter = 0;
        if (collection.has("multiphysics_time_scale"))
            tau_parameter = stod(multiphysics_ini.get("time_scale").get("multiphysics_time_scale"));
        configuration.Set_multiphysics_time_scale(tau_parameter);
    }

// IV
// sExternal stress state
    if (multiphysics_ini.has("stress_tensor")) {
        auto &collection = multiphysics_ini["stress_tensor"];
        if (collection.has("sxx"))
            sxx = stod(multiphysics_ini.get("stress_tensor").get("sxx"));
        if (collection.has("sxy"))
            sxy = stod(multiphysics_ini.get("stress_tensor").get("sxy"));
        if (collection.has("sxz"))
            sxz = stod(multiphysics_ini.get("stress_tensor").get("sxz"));

        if (collection.has("syx"))
            syx = stod(multiphysics_ini.get("stress_tensor").get("syx"));
        if (collection.has("syy"))
            syy = stod(multiphysics_ini.get("stress_tensor").get("syy"));
        if (collection.has("syz"))
            syz = stod(multiphysics_ini.get("stress_tensor").get("syz"));

        if (collection.has("szx"))
            szx = stod(multiphysics_ini.get("stress_tensor").get("szx"));
        if (collection.has("szy"))
            szy = stod(multiphysics_ini.get("stress_tensor").get("szy"));
        if (collection.has("szz"))
            szz = stod(multiphysics_ini.get("stress_tensor").get("szz"));
        //    Eigen::MatrixXd &new_external_stress_tensor
        Eigen::MatrixXd stress_tensor(3, 3); // Declare a 3x3 dynamic matrix
        stress_tensor << sxx, sxy, sxz,
                          syx, syy, syz,
                          szx, szy, szz; // dense matrix elements row by row
        configuration.Set_multiphysics_external_stress_tensor(stress_tensor);
        configuration.Get_multiphysics_vonMises_stress_in_MPa(stress_tensor);
        configuration.Get_multiphysics_pressure_in_MPa(stress_tensor);
    } // end of  if (multiphysics_ini.has("stress_tensor"))

 // V
 // Temperature
    if (multiphysics_ini.has("temperature")) {
        auto &collection = multiphysics_ini["temperature"];
        double tau_parameter = 0;
        if (collection.has("ambient_temperature"))
            new_ambient_temperature = stod(multiphysics_ini.get("temperature").get("ambient_temperature"));
        configuration.Set_multiphysics_temperature(new_ambient_temperature);
    }

    // module output
    if (multiphysics_ini.has("module_output"))
    {
        auto& collection = multiphysics_ini["module_output"];
        if (collection.has("module_log_file"))
        {
            log_file_output = multiphysics_ini.get("module_output").get("module_log_file");
        }
    }
    if (log_file_output == "ON") configuration.Set_is_multiphysics_log_file(true);

    double inclusion_stress_intensity_factor = 0.0, crack_stress_intensity_factor = 0.0;
    if (multiphysics_ini.has("microcracks")) {
        auto &collection = multiphysics_ini["microcracks"];

        if (collection.has("inclusion_stress_intensity_factor"))
            inclusion_stress_intensity_factor = stod(multiphysics_ini.get("microcracks").get("inclusion_stress_intensity_factor"));
            configuration.Set_multiphysics_inclusion_stress_intensity_factor(inclusion_stress_intensity_factor);

        if (collection.has("crack_stress_intensity_factor"))
            crack_stress_intensity_factor = stod(multiphysics_ini.get("microcracks").get("crack_stress_intensity_factor"));
            configuration.Set_multiphysics_crack_stress_intensity_factor(crack_stress_intensity_factor);
    }

    std::string grow_direction;
    int macrocrack_number = 0, number_of_crack_sizes = 0, grow_direction_id;
    double max_crack_lenghts = 0.0, min_crack_lenghts = 0.0, crack_stress_mode = 0.0;
    if (multiphysics_ini.has("macrocracks")) {
        auto &collection = multiphysics_ini["macrocracks"];

        if (collection.has("number_of_macrocracks"))
            macrocrack_number = stoi(multiphysics_ini.get("macrocracks").get("number_of_macrocracks"));
            configuration.Set_multiphysics_number_of_cracks(macrocrack_number);

        if (macrocrack_number > 0) {
            if (collection.has("crack_stress_mode"))
                crack_stress_mode = stod(multiphysics_ini.get("macrocracks").get("crack_stress_mode"));
            configuration.Set_multiphysics_crack_stress_mode(crack_stress_mode);

            if (collection.has("grow_direction"))
                grow_direction = multiphysics_ini.get("macrocracks").get("grow_direction");

            //grow_direction: '0' - for x-axis, '1' - for y-axis, '2' - for z-axis
            if(grow_direction == "xx"s) grow_direction_id = 0;
            else if (grow_direction == "yy"s) grow_direction_id = 1;
            else if (grow_direction == "zz"s) grow_direction_id = 2;
            else {
                cout << "ERROR: 'grow_direction' parameter in config/multiphysics.ini must be 'x', 'y' or 'z'. Please change accordingly!"<< endl;
                multiphysics_logfile_stream << "ERROR: 'grow_direction' parameter in config/multiphysics.ini must be 'x', 'y' or 'z'. Please change accordingly!"<< endl;
            }
            configuration.Set_multiphysics_crack_grow_direction(grow_direction_id); //; '1' - for x-axis, '2' - for y-axis, '3' - for z-axis

            if (collection.has("min_crack_lenghts"))
                min_crack_lenghts = stod(multiphysics_ini.get("macrocracks").get("min_crack_lenghts"));
            configuration.Set_multiphysics_min_crack_lenghts(min_crack_lenghts);

            if (collection.has("max_crack_lenghts"))
                max_crack_lenghts = stod(multiphysics_ini.get("macrocracks").get("max_crack_lenghts"));
            configuration.Set_multiphysics_max_crack_lenghts(max_crack_lenghts);

            if (collection.has("number_of_crack_sizes"))
                number_of_crack_sizes = stoi(multiphysics_ini.get("macrocracks").get("number_of_crack_sizes"));
            configuration.Set_multiphysics_number_of_crack_sizes(number_of_crack_sizes);

        } // end of if (macrocrack_number > 0)
    } // of  if (multiphysics_ini.has("macrocracks"))

/// Output to the screen/console
    cout << "______________________________________________________________________________________" << endl;
    cout << "The Multiphysics module specifications:\t\t" << endl;
    cout << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << ", y: "
         << std::get<1>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << ", z: " << std::get<2>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << endl;
    cout << "Characteristic time is \t\t\t"s << " tau: " << configuration.Get_multiphysics_time_scale() * pow(10, 6) << " [microseconds] "s << endl;

    // Homogeneous External Stress State
    double eq_stress_val = configuration.Get_multiphysics_pressure_in_MPa(configuration.Get_multiphysics_external_stress_tensor());
    double pressure_val = configuration.Get_multiphysics_vonMises_stress_in_MPa(configuration.Get_multiphysics_external_stress_tensor());
    cout << "Pressure is equal to \t\t\t"s << " P: "<< eq_stress_val << " [MPa] "s  << endl;
    cout << "Von Mises stress is equal to \t"s << " Sv: " << pressure_val << " [MPa] "s << endl;

    if (eq_stress_val || pressure_val > 0) {
        cout << "External Stress [MPa]: "s << endl;
        cout << configuration.Get_multiphysics_external_stress_tensor() << endl;
    }
    cout << endl;
    if(macrocrack_number > 0) {
        cout << "Number of macrocracks  \t\t\t\t"s << configuration.Get_multiphysics_number_of_cracks() << endl;
        cout << "Crack grow_direction (0->xx,1->yy,2->zz): \t\t"s << configuration.Get_multiphysics_crack_grow_direction() << endl;
        cout << "Crack mode  \t\t\t\t\t"s << configuration.Get_multiphysics_crack_stress_mode() << endl;
        cout << "MIN crack lenghts (fraction)  \t\t"s << configuration.Get_multiphysics_min_crack_lenghts() << endl;
        cout << "MAX crack lenghts (fraction)  \t\t"s << configuration.Get_multiphysics_max_crack_lenghts() << endl;
       if(configuration.Get_multiphysics_number_of_crack_sizes() > 1)
        cout << "Series of crack sizes (number)  \t"s << configuration.Get_multiphysics_number_of_crack_sizes() << endl;
    }
    cout << endl;
    cout << "Multiphysics module cpdlog_multiphysics.log file output:\t"s << configuration.Get_is_multiphysics_log_file() << endl;
    cout << endl;

// Output into .log file
    if(configuration.Get_is_multiphysics_log_file()) {
        multiphysics_logfile_stream
                << "______________________________________________________________________________________"
                << endl;
        multiphysics_logfile_stream << "The Multiphysics module specifications:\t\t" << endl;
        multiphysics_logfile_stream << "______________________________________________________________________________________" << endl;
        multiphysics_logfile_stream << "The Multiphysics module specifications:\t\t" << endl;
        multiphysics_logfile_stream << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << ", y: "
             << std::get<1>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << ", z: " << std::get<2>(configuration.Get_multiphysics_sample_dimensions()) << " [m] "s << endl;
        multiphysics_logfile_stream << "Characteristic time is \t\t\t"s << " tau: " << configuration.Get_multiphysics_time_scale() * pow(10, 6) << " [microseconds] "s << endl;

        // Homogeneous External Stress State
        double eq_stress_val = configuration.Get_multiphysics_pressure_in_MPa(configuration.Get_multiphysics_external_stress_tensor());
        double pressure_val = configuration.Get_multiphysics_vonMises_stress_in_MPa(configuration.Get_multiphysics_external_stress_tensor());
        multiphysics_logfile_stream << "Pressure is equal to \t\t\t"s << " P: "<< eq_stress_val << " [MPa] "s  << endl;
        multiphysics_logfile_stream << "Von Mises stress is equal to \t"s << " Sv: " << pressure_val << " [MPa] "s << endl;

        if (eq_stress_val || pressure_val > 0) {
            multiphysics_logfile_stream << "External Stress [MPa]: "s << endl;
            multiphysics_logfile_stream << configuration.Get_multiphysics_external_stress_tensor() << endl;
        }
        multiphysics_logfile_stream << endl;
        if(macrocrack_number > 0) {
            multiphysics_logfile_stream << "Number of macrocracks  \t\t\t\t"s << configuration.Get_multiphysics_number_of_cracks() << endl;
            multiphysics_logfile_stream << "Crack grow_direction (0->xx,1->yy,2->zz): \t\t"s << configuration.Get_multiphysics_crack_grow_direction() << endl;
            multiphysics_logfile_stream << "Crack mode  \t\t\t\t\t"s << configuration.Get_multiphysics_crack_stress_mode() << endl;
            multiphysics_logfile_stream << "MIN crack lenghts (fraction)  \t\t"s << configuration.Get_multiphysics_min_crack_lenghts() << endl;
            multiphysics_logfile_stream << "MAX crack lenghts (fraction)  \t\t"s << configuration.Get_multiphysics_max_crack_lenghts() << endl;
            if(configuration.Get_multiphysics_number_of_crack_sizes() > 1)
                multiphysics_logfile_stream << "Series of crack sizes (number)  \t"s << configuration.Get_multiphysics_number_of_crack_sizes() << endl;
        }
        multiphysics_logfile_stream << endl;
        multiphysics_logfile_stream << "Multiphysics module cpdlog_multiphysics.log file output:\t"s << configuration.Get_is_multiphysics_log_file() << endl;
        multiphysics_logfile_stream << endl;
        multiphysics_logfile_stream << endl;
    }

    return;
} /// end of the 'config_reader_multiphysics()' function

/// ================== # 3 # Initial PROCESSING module configuration - reading and output ==================
/*!
 * @details Read input parameters from the project file config/processing.ini necessary for the code execution. Print the read values to a screen and processing_logfile_stream --> cpdlog_processing.log file.
 * @param sequence_source_paths
 * @param max_fractions_vectors
 * @param max_cfractions_vectors
 * @param mu
 * @param sigma
 * @param bins_numb
 * @param ptype_vector
 * @param ctype_vector
 * @param pindex_vector
 * @return void
 */
//void config_reader_processing(std::vector<string> &sequence_source_paths, std::vector<vector<double>> &max_fractions_vectors, std::vector<vector<double>> &max_cfractions_vectors, double &mu, double &sigma, unsigned int &bins_numb, std::vector<string> &ptype_vector, std::vector<string> &ctype_vector, std::vector<double> &pindex_vector, bool &is_log_file) {
void config_reader_processing(Config &configuration) {
    std::string log_file_output;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "processing.ini"s);
    mINI::INIStructure processing_ini;
    file.read(processing_ini);

    int pindex, findex, eindex, nindex;
    std::string new_pp_mode, new_pf_mode, new_pe_mode, new_pn_mode;
    std::string new_ip_mode, new_if_mode, new_ie_mode, new_in_mode;
    std::string pseq_source, fseq_source, eseq_source, nseq_source;
    double p1_max = 0, p2_max = 0, p3_max = 0, f1_max = 0, f2_max = 0, f3_max = 0, e1_max = 0, e2_max = 0, e3_max = 0, n1_max = 0, n2_max = 0, n3_max = 0;
    int ptypes_number, ftypes_number, etypes_number, ntypes_number;
    double cp_max, cf_max, ce_max, cn_max;
    std::string cftypes_number_string, cetypes_number_string;
    bool is_processing_log_file = 0;
    std::vector<double> pp_max_vector, if_max_vector, ie_max_vector, pn_max_vector;

// I: cell types and max fractions and processing modes
//if (dim == 3) {
/// Polyhedrons
//processing_mode
    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];

        if (collection.has("polyhedron_types_number")){
            ptypes_number = stoi(processing_ini.get("polyhedrons").get("polyhedron_types_number")); // [2]
        }

        if (collection.has("pp_mode")) {
            new_pp_mode = processing_ini.get("polyhedrons").get("pp_mode");
            configuration.Set_processing_pp_mode(new_pp_mode);
        }

        if (collection.has("source")) {
            pseq_source = processing_ini.get("polyhedrons").get("source");
            configuration.Set_processing_p_source_path(pseq_source);
        }

        if (collection.has("p_multiplexity")) {
            pindex = stoi(processing_ini.get("polyhedrons").get("p_multiplexity"));
            configuration.Set_processing_p_multiplexity(pindex);
        }

// fractions
        if (collection.has("pmax_fraction1"))
            p1_max = stod(processing_ini.get("polyhedrons").get("pmax_fraction1"));
            pp_max_vector.push_back(p1_max);

        if (collection.has("pmax_fraction2"))
            p2_max = stod(processing_ini.get("polyhedrons").get("pmax_fraction2"));
            pp_max_vector.push_back(p2_max);

        if (collection.has("pmax_fraction3"))
            p3_max = stod(processing_ini.get("polyhedrons").get("pmax_fraction3"));
            pp_max_vector.push_back(p3_max);

    configuration.Set_processing_pp_max_fractions(pp_max_vector);
    }
    // } // end of dim == 3

/// Faces
//processing_mode
    if (processing_ini.has("faces")) {
        auto &collection = processing_ini["faces"];
        if (collection.has("face_types_number"))
            ftypes_number = stoi(processing_ini.get("faces").get("face_types_number"));

        if (collection.has("pf_mode"))
            new_pf_mode = processing_ini.get("faces").get("pf_mode");
        configuration.Set_processing_pf_mode(new_pf_mode);

        if (collection.has("f_multiplexity")) {
            findex = stoi(processing_ini.get("faces").get("f_multiplexity"));
            configuration.Set_processing_f_multiplexity(findex);
        }        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode

        if (collection.has("source")) {
            fseq_source = processing_ini.get("faces").get("source");
            configuration.Set_processing_f_source_path(fseq_source);
        }

// fractions
        std::vector<double> pf_max_vector;
        if (collection.has("fmax_fraction1"))
            f1_max = stod(processing_ini.get("faces").get("fmax_fraction1"));
            pf_max_vector.push_back(f1_max); // 2 - faces

        if (collection.has("fmax_fraction2"))
            f2_max =  stod(processing_ini.get("faces").get("fmax_fraction2"));
            pf_max_vector.push_back(f2_max); // 2 - faces

        if (collection.has("fmax_fraction3"))
            f3_max =  stod(processing_ini.get("faces").get("fmax_fraction3"));
            pf_max_vector.push_back(f3_max); // 2 - faces

        configuration.Set_processing_pf_max_fractions(pf_max_vector);

// induced structure
        if (collection.has("crack_types_number"))
            cftypes_number_string = processing_ini.get("faces").get("crack_types_number");

        if (collection.has("cf_mode")) {
            new_if_mode = processing_ini.get("faces").get("cf_mode");
            configuration.Set_processing_if_mode(new_if_mode);
        }

        if (collection.has("cfmax_fraction"))
            cf_max = stod(processing_ini.get("faces").get("cfmax_fraction"));
        if (cf_max > 0.0) if_max_vector.push_back(cf_max); // 2 - generated faces

        configuration.Set_processing_if_max_fractions(if_max_vector);

    }
/// Edges
//processing_mode
    if (processing_ini.has("edges")) {
        auto &collection = processing_ini["edges"];

        if (collection.has("edge_types_number"))
            etypes_number = stoi(processing_ini.get("edges").get("edge_types_number"));

        if (collection.has("pe_mode")) {
            new_pe_mode = processing_ini.get("edges").get("pe_mode");
            configuration.Set_processing_pe_mode(new_pe_mode);
        }

        if (collection.has("e_multiplexity")) {
            eindex = stoi(processing_ini.get("edges").get("e_multiplexity"));
            configuration.Set_processing_e_multiplexity(eindex);
        } // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode

        if (collection.has("source"))
            eseq_source = processing_ini.get("edges").get("source");
        configuration.Set_processing_e_source_path(eseq_source);

// fractions
        std::vector<double> pe_max_vector;
        if (collection.has("emax_fraction1"))
            e1_max = stod(processing_ini.get("edges").get("emax_fraction1"));
            pe_max_vector.push_back(e1_max); // 1 - edges

        if (collection.has("emax_fraction2"))
            e2_max = stod(processing_ini.get("edges").get("emax_fraction2"));
            pe_max_vector.push_back(e2_max); // 1 - edges

        if (collection.has("emax_fraction3"))
            e3_max = stod(processing_ini.get("edges").get("emax_fraction3"));
            pe_max_vector.push_back(e3_max); // 1 - edges

        configuration.Set_processing_pe_max_fractions(pe_max_vector);

/// Fracture for edges
        // induced structure
        if (collection.has("crack_types_number"))
            cetypes_number_string = processing_ini.get("edges").get("crack_types_number");

        if (collection.has("ce_mode")) {
            new_ie_mode = processing_ini.get("edges").get("ce_mode");
            configuration.Set_processing_ie_mode(new_ie_mode);
        }

        if (collection.has("cemax_fraction"))
            ce_max = stod(processing_ini.get("edges").get("cemax_fraction"));
        if (ce_max > 0.0) ie_max_vector.push_back(ce_max); // 1 - generated edges

        configuration.Set_processing_ie_max_fractions(ie_max_vector);
    }

/// Nodes
//processing_mode
    if (processing_ini.has("nodes")) {
        auto &collection = processing_ini["nodes"];

        if (collection.has("node_types_number"))
            ntypes_number = stoi(processing_ini.get("nodes").get("node_types_number"));

        if (collection.has("pn_mode")) {
            new_pn_mode = processing_ini.get("nodes").get("pn_mode");
            configuration.Set_processing_pn_mode(new_pn_mode);
        }

        if (collection.has("n_multiplexity")) {
            nindex = stoi(processing_ini.get("nodes").get("n_multiplexity"));
            configuration.Set_processing_n_multiplexity(nindex);
        } // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode

        if (collection.has("source"))
            nseq_source = processing_ini.get("nodes").get("source");
        configuration.Set_processing_n_source_path(nseq_source);

// fractions
        if (collection.has("nmax_fraction1"))
            n1_max = stod(processing_ini.get("nodes").get("nmax_fraction1"));
            pn_max_vector.push_back(n1_max); // 0 - nodes

        if (collection.has("nmax_fraction2"))
            n2_max = stod(processing_ini.get("nodes").get("nmax_fraction2"));
            pn_max_vector.push_back(n2_max); // 0 - nodes

        if (collection.has("nmax_fraction3"))
            n3_max = stod(processing_ini.get("nodes").get("nmax_fraction3"));
            pn_max_vector.push_back(n3_max); // 0 - nodes

        configuration.Set_processing_pn_max_fractions(pn_max_vector);
    }

// III: distribution
    if (processing_ini.has("distribution")) {
        auto& collection = processing_ini["distribution"];

        double distribution_mu = 0;
        if (collection.has("mu"))
            distribution_mu = stod(processing_ini.get("distribution").get("mu"));
        configuration.Set_processing_mu(distribution_mu);

        double distribution_sigma = 0;
        if (collection.has("sigma"))
            distribution_sigma = stod(processing_ini.get("distribution").get("sigma"));
        configuration.Set_processing_sigma(distribution_sigma);

        unsigned int distribution_bins_numb = 0;
        if (collection.has("bins_number"))
            distribution_bins_numb = stod(processing_ini.get("distribution").get("bins_number"));
        configuration.Set_processing_bins_number(distribution_bins_numb);
    }

    // Module output
    if (processing_ini.has("output")) {
        auto &collection = processing_ini["output"];
        std::string is_processing_log_out = "OFF";
        if (collection.has("polyhedron_types_number")) {
            is_processing_log_out =processing_ini.get("output").get("module_log_file");
        }
        if (is_processing_log_out == "ON") is_processing_log_file = 1; else is_processing_log_file = 0;
    }

    /// sequences
    std::vector<std::string> ptype_vector = {new_pn_mode, new_pe_mode, new_pf_mode, new_pp_mode};
//    configuration.Get_processing_config().sequence_source_paths = {nseq_source, eseq_source, fseq_source, pseq_source};
    vector<double> max_fractions_output(3, 0); // temporary vector serving as an output template for max fractions

/// Output to the screen/console
    if (configuration.processing_reader_switch) {

        cout << "The Processing module simulation type and initial parameters:\t\t" << endl;
        cout << endl;
        if (ptypes_number != 0) {
            // polyhedrons
            cout << "Processing p_type:\t"s << configuration.Get_processing_pp_mode() << "\t with p_index:\t"s
                 << configuration.Get_processing_p_multiplexity() << endl;
            if (configuration.Get_processing_pp_mode() == "L")
                cout << "mu = \t"s << configuration.Get_processing_config().mu << " and " << "sigma = \t"s
                     << configuration.Get_processing_config().sigma << endl;
            if (configuration.Get_processing_pp_mode() == "S")
                cout << "polyhedron sequence source: "s << pseq_source << endl;
            cout << "Number of polyhedron types:\t"s << ptypes_number << endl;
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pp_max_fractions().size() > 0 &&
                    configuration.Get_processing_pp_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pp_max_fractions()[i];
            cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
                 << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
            cout << endl;
        }
        if (ftypes_number != 0) {
            // faces
            cout << "Processing f_type:\t"s << configuration.Get_processing_pf_mode() << "\t with f_index:\t"s
                 << configuration.Get_processing_f_multiplexity() << endl;
            if (configuration.Get_processing_pf_mode() == "L")
                cout << "mu = \t"s << configuration.Get_processing_config().mu << " and " << "sigma = \t"s
                     << configuration.Get_processing_config().sigma << endl;
            if (configuration.Get_processing_pf_mode() == "S") cout << "face sequence source: "s << fseq_source << endl;
            cout << "Number of face types:\t"s << ftypes_number << endl;
// refill 0s
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pf_max_fractions().size() > 0 &&
                    configuration.Get_processing_pf_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pf_max_fractions()[i];
            cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
                 << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
            cout << endl;
        }

        if (if_max_vector.size() > 0) {
            cout << "Processing cf_mode:\t"s << configuration.Get_processing_if_mode() << endl;
            cout << "Number of generated face types:\t"s << cftypes_number_string << endl;
            if (configuration.Get_processing_if_mode() == "Km")
                cout << "Their maximum fractions:\t"s << configuration.Get_processing_if_max_fractions()[0] << endl;
        }

        if (etypes_number != 0) {
            //edges
            cout << "Processing e_type:\t"s << configuration.Get_processing_pe_mode() << "\twith e_index:\t"s
                 << configuration.Get_processing_e_multiplexity() << endl;
            if (configuration.Get_processing_pe_mode() == "L")
                cout << "mu = \t"s << configuration.Get_processing_config().mu << " and " << "sigma = \t"s
                     << configuration.Get_processing_config().sigma << endl;
            if (configuration.Get_processing_pe_mode() == "S")
                cout << "edges sequence source: "s << eseq_source << endl;
            cout << "Number of edge types:\t"s << etypes_number << endl;
// refill 0s
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pe_max_fractions().size() > 0 &&
                    configuration.Get_processing_pe_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pe_max_fractions()[i];
            cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
                 << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
            cout << endl;
        }

        if (ie_max_vector.size() > 0) {
            cout << "Processing c_edge_type:\t"s << configuration.Get_processing_in_mode() << endl;
            cout << "Number of edge crack types:\t"s << cetypes_number_string << endl;
            if (configuration.Get_processing_in_mode() == "Km")
                cout << "Their maximum fractions:\t"s << configuration.Get_processing_ie_max_fractions()[0] << endl;
        }

        if (ntypes_number != 0) {
            // nodes
            cout << "Processing n_type:\t"s << configuration.Get_processing_pn_mode() << "\twith n_index:\t"s
                 << configuration.Get_processing_n_multiplexity() << endl;
            if (configuration.Get_processing_pn_mode() == "L")
                cout << "mu = \t"s << configuration.Get_processing_config().mu << " and " << "sigma = \t"s
                     << configuration.Get_processing_config().sigma << endl;
            if (configuration.Get_processing_pn_mode() == "S")
                cout << "nodes sequence source: "s << nseq_source << endl;
            cout << "Number of node types:\t"s << ntypes_number << endl;

            // refill 0s
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pn_max_fractions().size() > 0 &&
                    configuration.Get_processing_pn_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pn_max_fractions()[i];
            cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
                 << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
            cout << "_________________________________________________" << endl << endl;
// Output to the screen/console
            cout << "Processing module cpdlog_processing.log file output:\t"s << log_file_output << endl << endl;
        }

// Output into .log file
        if (configuration.Get_is_processing_log_file()) {
            processing_logfile_stream << "The Processing module simulation type and initial parameters:\t\t" << endl;
            processing_logfile_stream << endl;
            if (ptypes_number != 0) {
                // polyhedrons
                processing_logfile_stream << "Processing p_type:\t"s << configuration.Get_processing_pp_mode()
                                          << "\t with p_index:\t"s
                                          << configuration.Get_processing_p_multiplexity() << endl;
                if (configuration.Get_processing_pp_mode() == "L")
                    processing_logfile_stream << "mu = \t"s << configuration.Get_processing_config().mu << " and "
                                              << "sigma = \t"s << configuration.Get_processing_config().sigma << endl;
                if (configuration.Get_processing_pp_mode() == "S")
                    processing_logfile_stream << "polyhedron sequence source: "s << pseq_source << endl;
                processing_logfile_stream << "Number of polyhedron types:\t"s << ptypes_number << endl;
                std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
                for (int i = 0; i < 3; ++i) {
                    if (configuration.Get_processing_pp_max_fractions().size() > 0 &&
                        configuration.Get_processing_pp_max_fractions()[i] > 0)
                        max_fractions_output.at(i) = configuration.Get_processing_pp_max_fractions()[i];
                }
                processing_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t"
                                          << max_fractions_output.at(1) << "\t\t" << max_fractions_output.at(2)
                                          << "\t\t" << endl;
                processing_logfile_stream << endl;
            }
            if (ftypes_number != 0) {
                // faces
                processing_logfile_stream << "Processing f_type:\t"s << configuration.Get_processing_pf_mode()
                                          << "\t with f_index:\t"s
                                          << configuration.Get_processing_f_multiplexity() << endl;
                if (configuration.Get_processing_pf_mode() == "L")
                    processing_logfile_stream << "mu = \t"s << configuration.Get_processing_config().mu << " and "
                                              << "sigma = \t"s << configuration.Get_processing_config().sigma << endl;
                if (configuration.Get_processing_pf_mode() == "S")
                    processing_logfile_stream << "face sequence source: "s << fseq_source << endl;
                processing_logfile_stream << "Number of face types:\t"s << ftypes_number << endl;
// refill 0s
                std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
                for (int i = 0; i < 3; ++i)
                    if (configuration.Get_processing_pf_max_fractions().size() > 0 &&
                        configuration.Get_processing_pf_max_fractions()[i] > 0)
                        max_fractions_output.at(i) = configuration.Get_processing_pf_max_fractions()[i];
                processing_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t"
                                          << max_fractions_output.at(1) << "\t\t" << max_fractions_output.at(2)
                                          << "\t\t" << endl;
                processing_logfile_stream << endl;
            }
            if (cftypes_number_string != "0") {
                processing_logfile_stream << "Processing c_face_type:\t"s << configuration.Get_processing_if_mode()
                                          << endl;
                processing_logfile_stream << "Number of face crack types:\t"s << cftypes_number_string << endl;
                if (configuration.Get_processing_if_mode() == "Km")
                    processing_logfile_stream << "Their maximum fractions:\t"s
                                              << configuration.Get_processing_if_max_fractions()[0] << endl;
            }

            //edges
            if (etypes_number != 0) {
                processing_logfile_stream << "Processing e_type:\t"s << configuration.Get_processing_pe_mode()
                                          << "\twith e_index:\t"s
                                          << configuration.Get_processing_e_multiplexity() << endl;
                if (configuration.Get_processing_pe_mode() == "L")
                    processing_logfile_stream << "mu = \t"s << configuration.Get_processing_config().mu << " and "
                                              << "sigma = \t"s << configuration.Get_processing_config().sigma << endl;
                if (configuration.Get_processing_pe_mode() == "S")
                    processing_logfile_stream << "edges sequence source: "s << eseq_source << endl;
                processing_logfile_stream << "Number of edge types:\t"s << etypes_number << endl;
// refill 0s
                std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
                for (int i = 0; i < 3; ++i)
                    if (configuration.Get_processing_pe_max_fractions().size() > 0 &&
                        configuration.Get_processing_pe_max_fractions()[i] > 0)
                        max_fractions_output.at(i) = configuration.Get_processing_pe_max_fractions()[i];
                processing_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t"
                                          << max_fractions_output.at(1)
                                          << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
                processing_logfile_stream << endl;
            }

            if (ie_max_vector.size() > 0) {
                processing_logfile_stream << "Processing c_edge_type:\t"s << configuration.Get_processing_ie_mode()
                                          << endl;
                processing_logfile_stream << "Number of edge crack types:\t"s << cetypes_number_string << endl;
                if (configuration.Get_processing_ie_mode() == "Km")
                    processing_logfile_stream << "Their maximum fractions:\t"s
                                              << configuration.Get_processing_pn_max_fractions()[0] << endl;
            }

            if (ntypes_number != 0) {
                // nodes
                processing_logfile_stream << "Processing n_type:\t"s << configuration.Get_processing_pn_mode()
                                          << "\twith n_index:\t"s
                                          << configuration.Get_processing_n_multiplexity() << endl;
                if (configuration.Get_processing_pn_mode() == "L")
                    processing_logfile_stream << "mu = \t"s << configuration.Get_processing_config().mu << " and "
                                              << "sigma = \t"s << configuration.Get_processing_config().sigma << endl;
                if (configuration.Get_processing_pn_mode() == "S")
                    processing_logfile_stream << "nodes sequence source: "s << nseq_source << endl;
                processing_logfile_stream << "Number of node types:\t"s << ntypes_number << endl;
            }
// refill 0s
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pe_max_fractions().size() > 0 &&
                    configuration.Get_processing_pe_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pe_max_fractions()[i];
            processing_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t"
                                      << max_fractions_output.at(1) << "\t\t" << max_fractions_output.at(2) << "\t\t"
                                      << endl;
            processing_logfile_stream << endl;
            processing_logfile_stream << "Processing n_type:\t"s << configuration.Get_processing_pn_mode()
                                      << "\twith n_index:\t"s
                                      << configuration.Get_processing_n_multiplexity() << endl;
            if (configuration.Get_processing_pn_mode() == "L")
                processing_logfile_stream << "mu = \t"s << configuration.Get_processing_config().mu << " and "
                                          << "sigma = \t"s << configuration.Get_processing_config().sigma << endl;
            processing_logfile_stream << "Number of node types:\t"s << ntypes_number << endl;
// refill 0s
            std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
            for (int i = 0; i < 3; ++i)
                if (configuration.Get_processing_pn_max_fractions().size() > 0 &&
                    configuration.Get_processing_pn_max_fractions()[i] > 0)
                    max_fractions_output.at(i) = configuration.Get_processing_pn_max_fractions()[i];
            processing_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t"
                                      << max_fractions_output.at(1) << "\t\t" << max_fractions_output.at(2) << "\t\t"
                                      << endl;
            processing_logfile_stream << "_________________________________________________" << endl << endl;
        }

        configuration.processing_reader_switch = false;
        } //     if (configuration.processing_reader_switch) {

    return;
} /// END of the 'config_reader_processing' function


/// ================== # 4 # Initial CHARACTERISATION module configuration - reading and output ==================
std::vector<double> config_reader_characterisation(std::vector<int> &charlabs_polyhedrons, std::vector<int> &charlabs_faces, std::vector<int> &charlabs_edges, std::vector<int> &charlabs_nodes, std::vector<int> &charlabs_laplacians, bool &is_log_file) {
    std::vector<double> config_characterisation_vector;
    std::string log_file_output;

// ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "characterisation.ini"s);
    mINI::INIStructure char_ini;
    file.read(char_ini);

/// Polyhedrons
    if (char_ini.has("polyhedrons_lab")) {
        auto& collection = char_ini["polyhedrons_lab"];
        if (collection.has("pl_active"))
        {
            charlabs_polyhedrons.push_back(stoi(char_ini.get("polyhedrons_lab").get("pl_active")));
        } }

    if (char_ini.has("polyhedrons_lab")) {
        auto& collection = char_ini["polyhedrons_lab"];
        if (collection.has("config_entropy"))
        {
            charlabs_polyhedrons.push_back(stoi(char_ini.get("polyhedrons_lab").get("config_entropy")));
        } }

    if (char_ini.has("polyhedrons_lab")) {
        auto& collection = char_ini["polyhedrons_lab"];
        if (collection.has("S_mean"))
        {
            charlabs_polyhedrons.push_back(stoi(char_ini.get("polyhedrons_lab").get("S_mean")));
        } }

    if (char_ini.has("polyhedrons_lab")) {
        auto& collection = char_ini["polyhedrons_lab"];
        if (collection.has("S_skew"))
        {
            charlabs_polyhedrons.push_back(stoi(char_ini.get("polyhedrons_lab").get("S_skew")));
        } }

/// Faces
    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("fl_active"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("fl_active")));
        } }

    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("config_entropy"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("config_entropy")));
        } }

    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("S_mean"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("S_mean")));
        } }

    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("S_skew"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("S_skew")));
        } }

    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("j_fractions"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("j_fractions")));
        } }

    if (char_ini.has("faces_lab")) {
        auto& collection = char_ini["faces_lab"];
        if (collection.has("d_fractions"))
        {
            charlabs_faces.push_back(stoi(char_ini.get("faces_lab").get("d_fractions")));
        } }

/// Edges
    if (char_ini.has("edges_lab")) {
        auto& collection = char_ini["edges_lab"];
        if (collection.has("el_active"))
        {
            charlabs_edges.push_back(stoi(char_ini.get("edges_lab").get("el_active"))); // [0]
        } }

    if (char_ini.has("edges_lab")) {
        auto& collection = char_ini["edges_lab"];
        if (collection.has("config_entropy"))
        {
            charlabs_edges.push_back(stoi(char_ini.get("edges_lab").get("config_entropy"))); // [1]
        } }

    if (char_ini.has("edges_lab")) {
        auto& collection = char_ini["edges_lab"];
        if (collection.has("S_mean"))
        {
            charlabs_edges.push_back(stoi(char_ini.get("edges_lab").get("S_mean"))); // [2]
        } }

    if (char_ini.has("edges_lab")) {
        auto& collection = char_ini["edges_lab"];
        if (collection.has("S_skew"))
        {
            charlabs_edges.push_back(stoi(char_ini.get("edges_lab").get("S_skew"))); // [3]
        } }

    if (char_ini.has("edges_lab")) {
        auto& collection = char_ini["edges_lab"];
        if (collection.has("analytical"))
        {
            charlabs_edges.push_back(stoi(char_ini.get("edges_lab").get("analytical"))); // [4]
        } }

/// Nodes
    if (char_ini.has("nodes_lab")) {
        auto& collection = char_ini["nodes_lab"];
        if (collection.has("nl_active"))
        {
            charlabs_nodes.push_back(stoi(char_ini.get("nodes_lab").get("nl_active")));
        } }

    if (char_ini.has("nodes_lab")) {
        auto& collection = char_ini["nodes_lab"];
        if (collection.has("config_entropy"))
        {
            charlabs_nodes.push_back(stoi(char_ini.get("nodes_lab").get("config_entropy")));
        } }

    if (char_ini.has("nodes_lab")) {
        auto& collection = char_ini["nodes_lab"];
        if (collection.has("S_mean"))
        {
            charlabs_nodes.push_back(stoi(char_ini.get("nodes_lab").get("S_mean")));
        } }

    if (char_ini.has("nodes_lab")) {
        auto& collection = char_ini["nodes_lab"];
        if (collection.has("S_skew"))
        {
            charlabs_nodes.push_back(stoi(char_ini.get("nodes_lab").get("S_skew")));
        } }

/// Laplacians
    if (char_ini.has("spectra_lab")) {
        auto& collection = char_ini["spectra_lab"];
        if (collection.has("calc_steps_numb"))
        {
            charlabs_laplacians.push_back(stoi(char_ini.get("spectra_lab").get("calc_steps_numb"))); // 0
        } }

    if (char_ini.has("spectra_lab")) {
        auto& collection = char_ini["spectra_lab"];
        if (collection.has("laplacians"))
        {
            charlabs_laplacians.push_back(stoi(char_ini.get("spectra_lab").get("laplacians"))); // 1
        } }

// module output
    if (char_ini.has("module_output")) {
        auto& collection = char_ini["module_output"];
        if (collection.has("module_log_file"))
        {
            log_file_output = char_ini.get("module_output").get("module_log_file");
        } }
    if (log_file_output == "ON") is_log_file = 1;

/// Console output
///    if (configuration.characterisation_reader_switch) {
        cout << "The Characterisation module simulation type and initial parameters:\t\t" << endl << endl;
        cout << "Polyhedrons lab ON/OFF:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
        cout << "Faces lab       ON/OFF:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
        cout << "Edges lab       ON/OFF:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
        cout << "Nodes lab       ON/OFF:\t"s << charlabs_nodes.at(0) << "\t\t" << endl << endl;

        if (charlabs_polyhedrons.at(0) == 1) { // Polyhedrons
            cout << "Polyhedrons configuration entropy:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
            cout << "Conf entropy mean part:\t"s << charlabs_polyhedrons.at(1) << "\t\t" << endl;
            cout << "Conf entropy skew part:\t"s << charlabs_polyhedrons.at(2) << "\t\t" << endl << endl;
        } // if(charlabs_polyhedrons.at(0) == 1)
        if (charlabs_faces.at(0) == 1) { // Faces
            cout << "Faces configuration entropy:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
            cout << "Conf entropy mean part:\t"s << charlabs_faces.at(1) << "\t\t" << endl;
            cout << "Conf entropy skew part:\t"s << charlabs_faces.at(2) << "\t\t" << endl;
            cout << "Edges fractions:\t"s << charlabs_faces.at(3) << "\t\t" << endl;
            cout << "Edges degree fractions:\t"s << charlabs_faces.at(4) << "\t\t" << endl << endl;
        } // if(charlabs_faces.at(0) == 1)
        if (charlabs_edges.at(0) == 1) { // Edges
            cout << "Edges configuration entropy:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
            cout << "Conf entropy Mean (-) Skew:\t"s << charlabs_edges.at(1) << "\t\t" << endl;
            cout << "Conf entropy mean part:\t"s << charlabs_edges.at(2) << "\t\t" << endl;
            cout << "Conf entropy skew part:\t"s << charlabs_edges.at(3) << "\t\t" << endl;
            cout << "Analytical solutions :\t"s << charlabs_edges.at(4) << "\t\t" << endl << endl;
        } // if(charlabs_edges.at(0) == 1)
        if (charlabs_nodes.at(0) == 1) { // Nodes
            cout << "Node configuration entropy:\t"s << charlabs_nodes.at(0) << "\t\t" << endl;
            cout << "Conf entropy mean part:\t"s << charlabs_nodes.at(1) << "\t\t" << endl;
            cout << "Conf entropy skew part:\t"s << charlabs_nodes.at(2) << "\t\t" << endl << endl;
        } // if(charlabs_polyhedrons.at(0) == 1)

        if (charlabs_laplacians.at(0) > 0) { // Laplacians lab
            cout << "Laplacians: number of calculation steps \t"s << charlabs_laplacians.at(0) << "\t\t" << endl;
            cout << endl;
        } // if(charlabs_laplacians.at(0) == 1)

        if (charlabs_laplacians.at(1) == 1) { // Laplacians lab
            cout << "Special cell Laplacians:\t"s << charlabs_laplacians.at(1) << "\t\t" << endl;
            cout << endl;
        } // if(charlabs_laplacians.at(0) == 1)

        cout << "Characterisation module cpdlog_characterisation.log file output:\t"s << log_file_output << endl
             << endl;

/// Output into Processing_Design.log file
        if (is_log_file) {
            characterisation_logfile_stream << "The Characterisation module simulation type and initial parameters:\t\t"
                                            << endl << endl;
            characterisation_logfile_stream << "Polyhedrons lab ON/OFF:\t"s << charlabs_polyhedrons.at(0) << "\t\t"
                                            << endl;
            characterisation_logfile_stream << "Faces lab       ON/OFF:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
            characterisation_logfile_stream << "Edges lab       ON/OFF:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
            characterisation_logfile_stream << "Nodes lab       ON/OFF:\t"s << charlabs_nodes.at(0) << "\t\t" << endl
                                            << endl;

            if (charlabs_polyhedrons.at(0) == 1) { // Polyhedrons
                characterisation_logfile_stream << "Polyhedrons configuration entropy:\t"s << charlabs_polyhedrons.at(0)
                                                << "\t\t" << endl;
                characterisation_logfile_stream << "Conf entropy mean part:\t"s << charlabs_polyhedrons.at(1) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy skew part:\t"s << charlabs_polyhedrons.at(2) << "\t\t"
                                                << endl << endl;
            } // if(charlabs_polyhedrons.at(0) == 1)
            if (charlabs_faces.at(0) == 1) { // Faces
                characterisation_logfile_stream << "Faces configuration entropy:\t"s << charlabs_faces.at(0) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy mean part:\t"s << charlabs_faces.at(1) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy skew part:\t"s << charlabs_faces.at(2) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Edges fractions:\t"s << charlabs_faces.at(3) << "\t\t" << endl;
                characterisation_logfile_stream << "Edges degree fractions:\t"s << charlabs_faces.at(4) << "\t\t"
                                                << endl << endl;
            } // if(charlabs_faces.at(0) == 1)
            if (charlabs_edges.at(0) == 1) { // Edges
                characterisation_logfile_stream << "Edges configuration entropy:\t"s << charlabs_edges.at(0) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy mean part:\t"s << charlabs_edges.at(1) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy skew part:\t"s << charlabs_edges.at(2) << "\t\t"
                                                << endl << endl;
            } // if(charlabs_edges.at(0) == 1)
            if (charlabs_nodes.at(0) == 1) { // Nodes
                characterisation_logfile_stream << "Node configuration entropy:\t"s << charlabs_nodes.at(0) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy mean part:\t"s << charlabs_nodes.at(1) << "\t\t"
                                                << endl;
                characterisation_logfile_stream << "Conf entropy skew part:\t"s << charlabs_nodes.at(2) << "\t\t"
                                                << endl << endl;
            } // if(charlabs_polyhedrons.at(0) == 1)

            if (charlabs_laplacians.at(0) > 0) { // Laplacians lab
                characterisation_logfile_stream << "Laplacians: number of calculation steps \t"s
                                                << charlabs_laplacians.at(0) << "\t\t" << endl << endl;
            } // if(charlabs_laplacians.at(0) == 1)

            if (charlabs_laplacians.at(1) == 1) { // Laplacians lab
                characterisation_logfile_stream << "Special cell Laplacians:\t"s << charlabs_laplacians.at(1) << "\t\t"
                                                << endl << endl;
            } // if(charlabs_laplacians.at(0) == 1)
        }

///        configuration.characterisation_reader_switch = false;
///    } // if(configuration.characterisation_reader_switch)

    return config_characterisation_vector;
} // END of config characterisation reader function

/// ================== # 5 # Initial DESIGN module configuration - reading and output ==================
void config_reader_design(Config &configuration) {

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "design.ini"s);
    mINI::INIStructure design_ini;
    file.read(design_ini);

    int cell_type;
    //'0' - nodes, '1' - edges, '2' - faces, '3' - grains
    std:: string design_goal;
    // 'min' for minimisation, and 'max' for xaimisation of the goal function output
    std:: string PCCDesign_type; // Design type ('G' for 'genetic')

    unsigned int population_size; // Number of 'creatures' in the initial and maybe future populations
    double mutation_rate; // fraction of mutation in the population
    double crossover_rate; // probability of acceptance
    double survival_rate; // The ratio (population_size/ survival_rate) gives the number of considered newly created State Vectors at each calculation step
    int max_generation_number; // maximal generation as a computation limit
    int genes_diversity; // number of different kinds of genes {0, 1, 2, 3, ..}
    bool is_design_log_file = false;
    std::string log_file_output; // output of results/design.log file

    if (design_ini.has("design_type")) {
        auto &collection = design_ini["design_type"];
        if (collection.has("design_mode")) {
            PCCDesign_type = design_ini.get("design_type").get("design_mode");
            configuration.Set_design_mode(PCCDesign_type);
        } }

    if (design_ini.has("genetic_algorithm")) {
        auto &collection = design_ini["genetic_algorithm"];

        if (collection.has("cell_type")) {
            cell_type = stoi(design_ini.get("genetic_algorithm").get("cell_type"));
            configuration.Set_design_cell_type(cell_type);
        }
        if (collection.has("design_goal")) {
            design_goal = design_ini.get("genetic_algorithm").get("design_goal");
            configuration.Set_design_goal(design_goal);
        }

        if (collection.has("genes_diversity")) {
            genes_diversity = stoi(design_ini.get("genetic_algorithm").get("genes_diversity"));
            configuration.Set_design_genes_diversity(genes_diversity);
        }
        if (collection.has("population_size")) {
            population_size = stoi(design_ini.get("genetic_algorithm").get("population_size"));
            configuration.Set_design_population_size(population_size);
        }
        if (collection.has("mutation_rate")) {
            mutation_rate = stod(design_ini.get("genetic_algorithm").get("mutation_rate"));
            configuration.Set_design_mutation_rate(mutation_rate);
        }
        if (collection.has("crossover_rate")) {
            crossover_rate = stod(design_ini.get("genetic_algorithm").get("crossover_rate"));
            configuration.Set_design_crossover_rate(crossover_rate);
        }
        if (collection.has("survival_rate")) {
            survival_rate = stod(design_ini.get("genetic_algorithm").get("survival_rate"));
            configuration.Set_design_survival_rate(survival_rate);
        }
        if (collection.has("max_generation_number")) {
            max_generation_number = stoi(design_ini.get("genetic_algorithm").get("max_generation_number"));
            configuration.Set_design_max_generation_number(max_generation_number);
        }

        // design module output
        if (design_ini.has("module_output")) {
            auto &collection = design_ini["module_output"];
            if (collection.has("module_log_file")) {
                log_file_output = design_ini.get("module_output").get("module_log_file");
            }
        }
        if (log_file_output == "ON") {
            configuration.Set_is_design_log_file(true);
            is_design_log_file = true;
        }

/// Console output
        if (configuration.design_reader_switch) {
            cout << "The Design module simulation type and initial parameters:\t\t" << endl << endl;
            cout << "Design mode:\t\t\t\t"s << PCCDesign_type << "\t\t" << endl;
            cout << "Design cell type:\t\t\t"s << cell_type << "\t\t" << endl;
            cout << "Diversity of gene types:\t"s << genes_diversity << "\t\t" << endl;
            cout << "Population size:\t\t\t"s << population_size << "\t\t" << endl;
            cout << "Mutation rate:\t\t\t\t"s << mutation_rate << "\t\t" << endl;
            cout << "Crossover rate\t\t\t\t"s << crossover_rate << "\t\t" << endl;
            cout << "Survival rate:\t\t\t\t"s << survival_rate << "\t\t" << endl;
            cout << "Max generation number:\t\t"s << max_generation_number << "\t\t" << endl;

            // Design logfile output
            if (is_design_log_file) {
                design_logfile_stream << "The Design module simulation type and initial parameters:\t\t" << endl
                                      << endl;
                design_logfile_stream << "Design mode:\t\t\t\t"s << PCCDesign_type << "\t\t" << endl;
                design_logfile_stream << "Design cell type:\t\t\t"s << cell_type << "\t\t" << endl;
                design_logfile_stream << "Diversity of gene types:\t"s << genes_diversity << "\t\t" << endl;
                design_logfile_stream << "Population size:\t\t\t"s << population_size << "\t\t" << endl;
                design_logfile_stream << "Mutation rate:\t\t\t\t"s << mutation_rate << "\t\t" << endl;
                design_logfile_stream << "Crossover rate\t\t\t\t"s << crossover_rate << "\t\t" << endl;
                design_logfile_stream << "Survival rate:\t\t\t\t"s << survival_rate << "\t\t" << endl;
                design_logfile_stream << "Max generation number:\t\t"s << max_generation_number << "\t\t" << endl;
            }
        } // end of 'if (design_ini.has("genetic_algorithm"))'

        configuration.design_reader_switch = false;
    } // if(configuration.design_reader_switch)

} // END of config_reader_design(Config &configuration)

/// ================== # 6 # Initial WRITER module configuration - reading and output ==================
void config_reader_writer(std::vector<int> &writer_specifications, bool &is_log_file) {
/// writer_specifications vector ::
    int    isSequencesOutput = 0;      // - >     [0]
    int    isDesignvectorsOutput = 0;  // - >     [1]
    int    isEnergiesOutput = 0;       // - >     [8]
    int    isSubPCCOutput = 0;        // - >     [9]
    int isEdgeConfEntropy = 0, isEdgeFractions = 0, isDegreeEdgeFractions = 0, isEdgeAnFractions = 0, isEdgeAnConfEntropies = 0; // [2], [3], [4], [5], [6]
    int isBetti = 0; // Laplacians lab  // - >     [7]
    std::string log_file_output;

// ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "writer.ini"s);
    mINI::INIStructure writer_ini;
    file.read(writer_ini);

// I
// sequences and designs output
    if (writer_ini.has("sequences")) {
        auto& collection = writer_ini["sequences"];
        if (collection.has("isSequencesOutput"))
        {
            isSequencesOutput = stoi(writer_ini.get("sequences").get("isSequencesOutput"));
        } }
        writer_specifications.push_back(isSequencesOutput); // [0]

    if (writer_ini.has("sequences")) {
        auto& collection = writer_ini["sequences"];
        if (collection.has("isDesignvectorsOutput"))
        {
            isDesignvectorsOutput = stoi(writer_ini.get("sequences").get("isDesignvectorsOutput"));
        } }
        writer_specifications.push_back(isDesignvectorsOutput); // [1]

// II Entropic
    if (writer_ini.has("entropic_edges")) {
        auto& collection = writer_ini["entropic_edges"];
        if (collection.has("isConfEntropy"))
        {
            isEdgeConfEntropy = stoi(writer_ini.get("entropic_edges").get("isConfEntropy"));
        } }
    writer_specifications.push_back(isEdgeConfEntropy); // [2]

    if (writer_ini.has("entropic_edges")) {
        auto& collection = writer_ini["entropic_edges"];
        if (collection.has("isFractions"))
        {
            isEdgeFractions = stoi(writer_ini.get("entropic_edges").get("isFractions"));
        } }
    writer_specifications.push_back(isEdgeFractions); // [3]

    if (writer_ini.has("entropic_edges")) {
        auto& collection = writer_ini["entropic_edges"];
        if (collection.has("isDegreeFractions"))
        {
            isDegreeEdgeFractions = stoi(writer_ini.get("entropic_edges").get("isDegreeFractions"));
        } }
    writer_specifications.push_back(isDegreeEdgeFractions); // [4]

    if (writer_ini.has("entropic_analytical")) {
        auto& collection = writer_ini["entropic_analytical"];
        if (collection.has("isEdgeFractions"))
        {
            isEdgeAnFractions = stoi(writer_ini.get("entropic_analytical").get("isEdgeFractions"));
        } }
    writer_specifications.push_back(isEdgeAnFractions); // [5]

    if (writer_ini.has("entropic_analytical")) {
        auto& collection = writer_ini["entropic_analytical"];
        if (collection.has("isEdgeConfEntropies"))
        {
            isEdgeAnConfEntropies = stoi(writer_ini.get("entropic_analytical").get("isEdgeConfEntropies"));
        } }
    writer_specifications.push_back(isEdgeAnConfEntropies); // [6]

// III Laplacians
    if (writer_ini.has("component_analysis")) {
        auto& collection = writer_ini["component_analysis"];
        if (collection.has("isBetti"))
        {
            isBetti = stoi(writer_ini.get("component_analysis").get("isBetti"));
        } }
    writer_specifications.push_back(isBetti); // [7]

// IV
// cell energies output
    if (writer_ini.has("energies")) {
        auto& collection = writer_ini["energies"];
        if (collection.has("isEnergiesOutput"))
        {
            isEnergiesOutput = stoi(writer_ini.get("energies").get("isEnergiesOutput"));
        } }
    writer_specifications.push_back(isEnergiesOutput); // [8]

    if (writer_ini.has("sequences")) {
        auto& collection = writer_ini["sequences"];
        if (collection.has("isSubcompexesOutput"))
        {
            isSubPCCOutput = stoi(writer_ini.get("sequences").get("isSubcompexesOutput"));
        } }
    writer_specifications.push_back(isSubPCCOutput); // [9]

    // module output
     if (writer_ini.has("module_output")) {
         auto& collection = writer_ini["module_output"];
         if (collection.has("module_log_file"))
         {
             log_file_output = writer_ini.get("module_output").get("module_log_file");
         } }
     if (log_file_output == "ON") is_log_file = 1;

/// Output to the screen/console
///    if(configuration.writer_reader_switch) {
        cout << "The Writer module specifications:\t\t" << endl;
        cout << "Sequences output \t\t\t\t\t"s << writer_specifications.at(0) << endl;
        cout << "Design vectors output \t\t\t\t"s << writer_specifications.at(1) << endl;
//        cout << "Configuration Edges entropy \t\t"s << writer_specifications.at(2) << endl;
//        cout << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(3) << endl;
//        cout << "Special Edge degree fractions \t\t"s << writer_specifications.at(4) << endl;
//        cout << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
//        cout << "Analytical Edge degree fractions \t"s << writer_specifications.at(5) << endl;
//        cout << "Analytical Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
        cout << "Laplacians and Betti numbers \t\t"s << writer_specifications.at(7) << endl;
        cout << "Cell Energies \t\t\t\t\t\t"s << writer_specifications.at(8) << endl;
        cout << "Subcomplexes \t\t\t\t\t\t"s << writer_specifications.at(9) << endl << endl;
        cout << "Writer module cpdlog_writer.log file output:\t"s << log_file_output << endl;

/// Output into .log file
        if (is_log_file) {
            writer_logfile_stream << "The Writer module specifications:\t\t" << endl;
            writer_logfile_stream << "Sequences output \t\t\t\t\t"s << writer_specifications.at(0) << endl;
            writer_logfile_stream << "Design vectors output \t\t\t\t"s << writer_specifications.at(1) << endl;
//            writer_logfile_stream << "Configuration Edges entropy \t\t\t"s << writer_specifications.at(2) << endl;
//            writer_logfile_stream << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(3) << endl;
//            writer_logfile_stream << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(4) << endl;
//            writer_logfile_stream << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
//            writer_logfile_stream << "Analytical Edge degree fractions \t"s << writer_specifications.at(5) << endl;
//            writer_logfile_stream << "Analytical Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
            writer_logfile_stream << "Laplacians and Betti numbers \t\t"s << writer_specifications.at(7) << endl;
            writer_logfile_stream << "Cell Energies \t\t\t\t\t\t"s << writer_specifications.at(8) << endl;
            writer_logfile_stream << "Subcomplexes \t\t\t\t\t\t"s << writer_specifications.at(9) << endl << endl;
        }

///        configuration.writer_reader_switch = false;
///    } // if(configuration.writer_reader_switch)

    return;
} /// END of config_reader_writer function

/// ================== # 7 # Initial KINETICS module configuration - kinetic (time-dependent) processes ==================
/*!
 * @details config_reader_kinetics :: read input parameters from the project file config/kinetics.ini necessary for the Kinetics module execution.
 * @param Config configuration
 * @return void
 */
void config_reader_kinetics(Config &configuration) {
    std::string log_file_output;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "kinetics.ini"s);
    mINI::INIStructure kinetics_ini;
    file.read(kinetics_ini);

    std::string new_nk_mode, new_ek_mode, new_fk_mode, new_pk_mode;
//kinetics type
    if (kinetics_ini.has("kinetic_mode")) {
        auto& collection = kinetics_ini["kinetic_mode"];
        if (collection.has("node_kinetic_mode")) {
            new_nk_mode = kinetics_ini.get("kinetic_mode").get("node_kinetic_mode");
            configuration.Set_kinetics_nk_mode(new_nk_mode);
        }
        if (collection.has("edge_kinetic_mode")) {
            new_ek_mode = kinetics_ini.get("kinetic_mode").get("edge_kinetic_mode");
            configuration.Set_kinetics_ek_mode(new_ek_mode);
        }
        if (collection.has("face_kinetic_mode")) {
            new_fk_mode = kinetics_ini.get("kinetic_mode").get("face_kinetic_mode");
            configuration.Set_kinetics_fk_mode(new_fk_mode);
        }
        if (collection.has("polyhedra_kinetic_mode")) {
            new_pk_mode = kinetics_ini.get("kinetic_mode").get("polyhedra_kinetic_mode");
            configuration.Set_kinetics_pk_mode(new_pk_mode);
        }
    }
    // Material ID from the CPD code Database
    std::string mat_id;
    if (kinetics_ini.has("general")) {
        auto &collection = kinetics_ini["general"];
        if (collection.has("mat_id")) {
            mat_id = kinetics_ini.get("general").get("mat_id");
            configuration.Set_kinetics_material_id(mat_id);
        }
        double time_constant;
        if (collection.has("kinetic_time_scale")) {
            time_constant = stod(kinetics_ini.get("general").get("kinetic_time_scale"));
            configuration.Set_kinetics_time_scale(time_constant);
        }
    }
        //corrosion
    if (kinetics_ini.has("corrosion")) {
        auto &collection = kinetics_ini["corrosion"];
        double corrosion_rate, corrosion_activation_volume;
        if (collection.has("corrosion_rate_coeff")) {
            corrosion_rate = stod(kinetics_ini.get("corrosion").get("corrosion_rate_coeff"));
            configuration.Set_kinetics_corrosion_rate_scale(corrosion_rate);
        }
        if (collection.has("corrosion_activation_volume")) {
            corrosion_activation_volume = stod(kinetics_ini.get("corrosion").get("corrosion_activation_volume"));
            configuration.Set_kinetics_corrosion_activation_volume(corrosion_activation_volume);
        }
    }

    // irradiation
    double rx = 0.0, rz = 0.0, ry = 0.0;
    double beam_energy_flux = 0; //energy flux in [J/s]
    double beam_current = 0; //beam current in [particles/s]
    double energy_dissipation_rate = 0; // energy dissipation of a beam with material depth
    double irradiation_damage_rate = 0; // an additional coefficient of the irradition damage model
    double observation_time = 0; // time of the module output

    if (kinetics_ini.has("irradiation")) {
        auto &collection = kinetics_ini["irradiation"];
        if (collection.has("rx"))
            rx = stod(kinetics_ini.get("irradiation").get("rx"));
        if (collection.has("rz"))
            rz = stod(kinetics_ini.get("irradiation").get("rz"));
        ry = std::sqrt(1.0 - std::pow(rx,2.0) - std::pow(rz,2.0));
        std::tuple<double, double, double> new_beam_direction = make_tuple(rx,ry,rz);
                configuration.Set_kinetics_beam_direction(new_beam_direction);

        if (collection.has("irradiation_damage_rate")) {
            irradiation_damage_rate = stod(kinetics_ini.get("irradiation").get("irradiation_damage_rate"));
            configuration.Set_kinetics_irradiation_damage_rate(irradiation_damage_rate);
        }

        if (collection.has("beam_energy_flux")) {
            beam_energy_flux = stod(kinetics_ini.get("irradiation").get("beam_energy_flux"));
            configuration.Set_kinetics_beam_energy_flux(beam_energy_flux);
        }
        if (collection.has("beam_current")) {
            beam_current = stod(kinetics_ini.get("irradiation").get("beam_current"));
            configuration.Set_kinetics_beam_current(beam_current);
        }
        if (collection.has("energy_dissipation_rate")) {
            energy_dissipation_rate = stod(kinetics_ini.get("irradiation").get("energy_dissipation_rate"));
            configuration.Set_kinetics_energy_dissipation_rate(energy_dissipation_rate);
        }
        if (collection.has("observation_time")) {
            observation_time = stod(kinetics_ini.get("irradiation").get("observation_time"));
            configuration.Set_kinetics_observation_time(observation_time);
        }
    }

    // module output
    if (kinetics_ini.has("module_output")) {
        auto& collection = kinetics_ini["module_output"];
        if (collection.has("module_log_file"))
        {
            log_file_output = kinetics_ini.get("module_output").get("module_log_file");
        } }
    if (log_file_output == "ON") configuration.Set_is_kinetics_log_file(true);

/// Output to the screen/console
    if(configuration.kinetics_reader_switch) {
        cout << "The Kinetics module type and initial parameters:\t\t" << endl;
        if (configuration.Get_kinetics_nk_mode() != "N"s)
            cout << "Node kinetics type:\t\t\t\t"s << configuration.Get_kinetics_nk_mode() << endl;
        if (configuration.Get_kinetics_ek_mode() != "N"s)
            cout << "Edge kinetics type:\t\t\t\t"s << configuration.Get_kinetics_ek_mode() << endl;
        if (configuration.Get_kinetics_fk_mode() != "N"s)
            cout << "Face kinetics type:\t\t\t\t"s << configuration.Get_kinetics_fk_mode() << endl;
        if (configuration.Get_kinetics_pk_mode() != "N"s)
            cout << "Volume kinetics type:\t\t\t\t"s << configuration.Get_kinetics_pk_mode() << endl;
        cout << "Material's ID:\t\t\t\t\t"s << configuration.Get_kinetics_material_id() << endl;
        cout << "Kinetic time scale:\t\t\t\t"s << configuration.Get_kinetics_time_scale() << endl;
        if (configuration.Get_kinetics_nk_mode() == "C"s || configuration.Get_kinetics_ek_mode() == "C"s ||
            configuration.Get_kinetics_fk_mode() == "C"s || configuration.Get_kinetics_pk_mode() == "C"s)
            cout << "Corrosion rate:\t\t\t\t\t"s << configuration.Get_kinetics_corrosion_rate_scale() << endl;
        if (configuration.Get_kinetics_nk_mode() == "I"s || configuration.Get_kinetics_ek_mode() == "I"s ||
            configuration.Get_kinetics_fk_mode() == "I"s || configuration.Get_kinetics_pk_mode() == "I"s) {
            cout << "Energy flux of the beam:\t\t"s << configuration.Get_kinetics_beam_energy_flux() << endl;
            cout << "Beam current:\t\t\t\t\t"s << configuration.Get_kinetics_beam_current() << endl;
            cout << "Energy dissipation rate:\t\t"s << configuration.Get_kinetics_energy_dissipation_rate() << endl;
        }
        if (configuration.Get_is_kinetics_log_file()) {
            cout << "Kinetics *.log file output is \tON"s << endl;

            kinetics_logfile_stream << "The Kinetics module type and initial parameters:\t\t" << endl;
            if (configuration.Get_kinetics_nk_mode() != "N"s)
                kinetics_logfile_stream << "Node kinetics type:\t\t\t\t"s << configuration.Get_kinetics_nk_mode()
                                        << endl;
            if (configuration.Get_kinetics_ek_mode() != "N"s)
                kinetics_logfile_stream << "Edge kinetics type:\t\t\t\t"s << configuration.Get_kinetics_ek_mode()
                                        << endl;
            if (configuration.Get_kinetics_fk_mode() != "N"s)
                kinetics_logfile_stream << "Face kinetics type:\t\t\t\t"s << configuration.Get_kinetics_fk_mode()
                                        << endl;
            if (configuration.Get_kinetics_pk_mode() != "N"s)
                kinetics_logfile_stream << "Volume kinetics type:\t\t\t\t"s << configuration.Get_kinetics_pk_mode()
                                        << endl;
            kinetics_logfile_stream << "Material's ID:\t\t\t\t\t"s << configuration.Get_kinetics_material_id() << endl;
            kinetics_logfile_stream << "Kinetic time scale:\t\t\t\t"s << configuration.Get_kinetics_time_scale()
                                    << endl;
            if (configuration.Get_kinetics_nk_mode() == "C"s || configuration.Get_kinetics_ek_mode() == "C"s ||
                configuration.Get_kinetics_fk_mode() == "C"s || configuration.Get_kinetics_pk_mode() == "C"s)
                kinetics_logfile_stream << "Corrosion rate:\t\t\t\t\t"s
                                        << configuration.Get_kinetics_corrosion_rate_scale() << endl;
            if (configuration.Get_kinetics_nk_mode() == "I"s || configuration.Get_kinetics_ek_mode() == "I"s ||
                configuration.Get_kinetics_fk_mode() == "I"s || configuration.Get_kinetics_pk_mode() == "I"s) {
                kinetics_logfile_stream << "Energy flux of the beam:\t\t"s
                                        << configuration.Get_kinetics_beam_energy_flux() << endl;
                kinetics_logfile_stream << "Beam current:\t\t\t\t\t"s << configuration.Get_kinetics_beam_current()
                                        << endl;
                kinetics_logfile_stream << "Energy dissipation rate:\t\t"s
                                        << configuration.Get_kinetics_energy_dissipation_rate() << endl;
            }

        } else
            cout << "Kinetics *.log file output is \tOFF"s << endl;

        configuration.kinetics_reader_switch = false;
    } // if(configuration.kinetics_reader_switch)

    return;
} /// END of the 'config_reader_kinetics() function
