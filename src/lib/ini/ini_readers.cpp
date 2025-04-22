/// Author: Dr Elijah Borodin (2023)
/// Manchester, UK
/// Library of specific functions related to the PCC Processing Design code for reading its *.ini files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Simple reader for *.ini files and a specific CPD code-related library for reading its particular *.ini files ( downloaded from https://github.com/pulzed/mINI )
#include "../ini/ini.h"

///* ------------------------------------------------------------------------------- *
///* Attached user-defined C++ libraries (must be copied in the directory for STL):
///* ------------------------------------------------------------------------------- *
/// Eigen source: https://eigen.tuxfamily.org/ (2024)
#include <Eigen/Dense>

using namespace std; // standard namespace

extern std::string output_dir;
extern std::ofstream Out_logfile_stream;

/// ================== # 1 # Initial configuration - reading and output ==================
std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string &cell_complex_standard, std::string &main_type) {

    std::vector<int> res(7,0);
/// [0] - > dim, [1] -> isSubcomplex, [2] -> isProcessing, [3] -> isCharacterisation, [4] -> isMultiphysics, [5] -> isKinetic, [6] -> isWriter
    bool isSubcomplexON = 0, isProcessingON = 0, isCharacterisationON = 0, isKineticON = 0, isMultiphysicsON = 0, isWriterON = 0;
    std::string isSubcomplex, isProcessing, isCharacterisation, isKinetic, isMultiphysics, isWriter;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "main.ini"s);
//    mINI::INIFile file(source_path + "main_2D.ini"s);

    mINI::INIStructure main_ini;
    file.read(main_ini);

// 0
//    if (main_ini.has("execution_type")) {
//        auto& collection = main_ini["execution_type"];
//        if (collection.has("e_type"))
//        {
//            e_type = main_ini.get("execution_type").get("e_type");
//        } }

// I
    if (main_ini.has("simulation_mode")) {
        auto& collection = main_ini["simulation_mode"];
        if (collection.has("mode"))
        {
            if (main_ini.get("simulation_mode").get("mode") == "LIST" || main_ini.get("simulation_mode").get("mode") == "TUTORIAL" ||
                    main_ini.get("simulation_mode").get("mode") == "PERFORMANCE_TEST" || main_ini.get("simulation_mode").get("mode") == "TASK")
            {
                main_type = main_ini.get("simulation_mode").get("mode");
                }
            else throw std::invalid_argument("ERROR in ../src/ini/ini_readers.cpp: WRONG TYPE OF THE 'simulation_mode' IN ../config/main.ini FILE; Please change the mode to one of the allowed: 'LIST', 'TUTORIAL', 'PERFORMANCE_TEST' or 'TASK' "s);
        }
    }
// II
    std::string problem_dimension;
    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("dim")) {
            if(stoi(main_ini.get("general").get("dim")) == 1 || stoi(main_ini.get("general").get("dim")) == 2 || stoi(main_ini.get("general").get("dim")) == 3) {
                problem_dimension = main_ini.get("general").get("dim");
            }
            else throw std::invalid_argument("ERROR in ../src/ini/ini_readers.cpp: WRONG TYPE OF DIMENSION 'dim' IN ../config/main.ini FILE; Please change the [general] dim parameter to one of the allowed: 1, 2 or 3 "s);
        }    }
    res.at(0) = stoi(problem_dimension); // res[0]

    // III
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Subcomplex"))
            isSubcomplex = main_ini.get("modules").get("PCC_Subcomplex");
    }

    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Multiphysics"))
            isMultiphysics = main_ini.get("modules").get("PCC_Multiphysics");
    }

    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Processing"))
            isProcessing = main_ini.get("modules").get("PCC_Processing");
    }

    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Characterisation"))
            isCharacterisation = main_ini.get("modules").get("PCC_Characterisation");
    }

    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Writer"))
            isWriter = main_ini.get("modules").get("PCC_Writer");
    }

    /// forming the output RES vector
// ON/OFF IDs
    if (isSubcomplex == "ON") { isSubcomplexON = 1; res.at(1) = 1; } else res.at(1) = 0; // res[1] - Section -> ConfigVector.at(1) in main.cpp
    if (isProcessing == "ON") { isProcessingON = 1; res.at(2) = 1; } else res.at(2) = 0; // res[2] - Processing -> ConfigVector.at(2) in main.cpp
    if (isCharacterisation == "ON") { isCharacterisationON = 1; res.at(3) = 1; } else res.at(3) = 0; // res[3] - Characterisation -> ConfigVector.at(3) in main.cpp
    if (isMultiphysics == "ON") { isMultiphysicsON = 1; res.at(4) = 1; } else res.at(4) = 0; // res[4] - Multiphysics -> ConfigVector.at(4) in main.cpp
    if (isWriter == "ON") { isWriterON = 1; res.at(6) = 1; } else res.at(6) = 0; // res[6] - Writer -> ConfigVector.at(6) in main.cpp

    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("source_dir"))
            source_dir = main_ini.get("general").get("source_dir");
    }

    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("pcc_standard"))
            cell_complex_standard = main_ini.get("general").get("pcc_standard");
    } // like 'pcc1s'

    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("output_dir"))
            output_dir = main_ini.get("general").get("output_dir");
    }

/// Output to the screen/console
    cout << "The problem dimension that is the maximum value k_max of k-cells in the PCC\t\t|\t\t"s << "dim = " << res.at(0) << endl;
    cout << "Simulation mode:\t"s << "\t" << main_type << endl;
    cout << "Output directory:\t"s << "\t" << output_dir << endl;
    cout << "PCC source directory:\t"s << source_dir << endl;
    cout << "PCC standard ID:\t\t"s << cell_complex_standard << endl;
    cout << endl;
    if (isSubcomplexON == 1) cout << "ON    | PCC_Subcomplex"s << endl;
    else cout << "OFF    | PCC_Subcomplex"s << endl;
    if (isMultiphysicsON == 1) cout << "ON    | PCC_Multiphysics"s << endl;
    else cout << "OFF    | PCC_Multiphysics"s << endl;
    if (isProcessingON == 1) cout << "ON    | PCC_Processing"s << endl;
    else cout << "OFF    | PCC_Processing"s << endl;
    if (isCharacterisationON == 1) cout << "ON    | PCC_Characterisation"s << endl;
    else cout << "OFF    | PCC_Characterisation"s << endl;
    if (isWriterON == 1) cout << "ON    | PCC_Writer"s << endl;
    else cout << "OFF    | PCC_Writer"s << endl;
    cout << endl;

/// Output into .log file
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::trunc); // this *.log stream will be closed at the end of the main function

    Out_logfile_stream << endl;
    Out_logfile_stream << "The problem dimension that is the maximum value k_max of k-cells in the PCC:\t\t|\t\t"s << "dim = " << res.at(0) << endl;
    Out_logfile_stream << "PCC standard ID:\t\t"s << cell_complex_standard << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Simulation mode:\t"s << "\t" << main_type << endl;
    Out_logfile_stream << "Output directory:\t"s << "\t" << output_dir << endl;
    Out_logfile_stream << "PCC source directory:\t"s << "\t" << source_dir << endl;

    Out_logfile_stream << endl;
    if (isSubcomplexON == 1) Out_logfile_stream << "ON    | PCC_Subcomplex"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Subcomplex"s << endl;
    if (isMultiphysicsON == 1) cout << "ON    | PCC_Multiphysics"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Multiphysics"s << endl;
    if (isProcessingON == 1) Out_logfile_stream << "ON    | PCC_Processing"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Processing"s << endl;
    if (isCharacterisationON == 1) Out_logfile_stream << "ON    | PCC_Characterisation"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Characterisation"s << endl;
    if (isWriterON == 1) Out_logfile_stream << "ON    | PCC_Writer"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Writer"s << endl;
    Out_logfile_stream << endl;

    Out_logfile_stream.close();

    return res;
} /// END of config_reader_main() function

/// ================== # 2 # Initial PROCESSING module configuration - reading and output ==================
void config_reader_processing(std::string &source_path, std::vector<string> &sequence_source_paths, std::vector<vector<double>> &max_fractions_vectors, std::vector<vector<double>> &max_cfractions_vectors, double &mu, double &sigma, unsigned int &bins_numb, std::vector<string> &ptype_vector, std::vector<string> &ctype_vector, std::vector<double> &pindex_vector, std::ofstream &Out_logfile_stream) {
    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "processing.ini"s);
//    mINI::INIFile file(source_path + "processing_2D.ini"s);
    mINI::INIStructure processing_ini;
    file.read(processing_ini);

// I: cell types and max fractions and processing modes
//if (dim == 3) {
/// Polyhedrons
//processing_mode
    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("pp_mode")) {
            ptype_vector.at(3) = processing_ini.get("polyhedrons").get("pp_mode");
        }
    }

    string pseq_source;
    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("source"))
            pseq_source = processing_ini.get("polyhedrons").get("source");
    }

    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("p_multiplexity"))
            pindex_vector.at(3) = stod(processing_ini.get("polyhedrons").get("p_multiplexity"));
        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode
    }

    string ptypes_number_string;
    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("polyhedron_types_number"))
            ptypes_number_string = processing_ini.get("polyhedrons").get("polyhedron_types_number"); // [2]
    }

// fractions
    string p1_max, p2_max, p3_max;
    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction1"))
            p1_max = processing_ini.get("polyhedrons").get("pmax_fraction1");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p1_max)); // 3 - polyhedra

    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction2"))
            p2_max = processing_ini.get("polyhedrons").get("pmax_fraction2");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p2_max)); // 3 - polyhedra

    if (processing_ini.has("polyhedrons")) {
        auto &collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction3"))
            p3_max = processing_ini.get("polyhedrons").get("pmax_fraction3");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p3_max)); // 3 - polyhedra
// } // end of dim == 3

/// Faces
//processing_mode
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("pf_mode"))
            ptype_vector.at(2) = processing_ini.get("faces").get("pf_mode");
         }

    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("f_multiplexity"))
            pindex_vector.at(2) = stod(processing_ini.get("faces").get("f_multiplexity"));
        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode
    }

    string ftypes_number_string;
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("face_types_number"))
            ftypes_number_string = processing_ini.get("faces").get("face_types_number");
    }

    string fseq_source;
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("source"))
            fseq_source = processing_ini.get("faces").get("source");
    }

// fractions
    string f1_max, f2_max, f3_max;
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("fmax_fraction1"))
            f1_max = processing_ini.get("faces").get("fmax_fraction1");
    }
    if (stoi(ftypes_number_string) > 0) max_fractions_vectors.at(2).push_back(stod(f1_max)); // 2 - faces

    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("fmax_fraction2"))
            f2_max = processing_ini.get("faces").get("fmax_fraction2");
    }
    if (stoi(ftypes_number_string) > 0) max_fractions_vectors.at(2).push_back(stod(f2_max)); // 2 - faces

    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("fmax_fraction3"))
            f3_max = processing_ini.get("faces").get("fmax_fraction3");
    }
    if (stoi(ftypes_number_string) > 0) max_fractions_vectors.at(2).push_back(stod(f3_max)); // 2 - faces

// induced structure
    string cftypes_number_string;
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("crack_types_number"))
            cftypes_number_string = processing_ini.get("faces").get("crack_types_number");
    }

    if (processing_ini.has("faces")) { /// FRACTURE mode
        auto& collection = processing_ini["faces"];
        if (collection.has("cf_mode"))
            ctype_vector.at(2) = processing_ini.get("faces").get("cf_mode");
    }
    string cf_max;
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("cfmax_fraction"))
            cf_max = processing_ini.get("faces").get("cfmax_fraction");
    }
    if (stoi(cftypes_number_string) > 0) max_cfractions_vectors.at(2).push_back(stod(cf_max)); // 2 - faces

/// Edges
//processing_mode
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("pe_mode"))
        {
            ptype_vector.at(1) = processing_ini.get("edges").get("pe_mode");
        } }

    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("e_multiplexity"))
            pindex_vector.at(1) = stod(processing_ini.get("edges").get("e_multiplexity"));
        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode
    }

    string eseq_source;
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("source"))
            eseq_source = processing_ini.get("edges").get("source");
    }

    string etypes_number_string;
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("edge_types_number"))
            etypes_number_string = processing_ini.get("edges").get("edge_types_number");
    }

// fractions
    string e1_max, e2_max, e3_max;
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("emax_fraction1"))
            e1_max = processing_ini.get("edges").get("emax_fraction1");
    }
    if (stoi(etypes_number_string) > 0) max_fractions_vectors.at(1).push_back(stod(e1_max)); // 1 - edges

    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("emax_fraction2"))
            e2_max = processing_ini.get("edges").get("emax_fraction2");
    }
    if (stoi(etypes_number_string) > 0) max_fractions_vectors.at(1).push_back(stod(e2_max)); // 1 - edges

    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("emax_fraction3"))
            e3_max = processing_ini.get("edges").get("emax_fraction3");
    }
    if (stoi(etypes_number_string) > 0) max_fractions_vectors.at(1).push_back(stod(e3_max)); // 1 - edges

/// Fracture for edges
    // induced structure
    string cetypes_number_string;
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("crack_types_number"))
            cetypes_number_string = processing_ini.get("edges").get("crack_types_number");
    }

    if (processing_ini.has("edges")) { /// FRACTURE mode
        auto& collection = processing_ini["edges"];
        if (collection.has("ce_mode"))
            ctype_vector.at(1) = processing_ini.get("edges").get("ce_mode");
    }
    string ce_max;
    if (processing_ini.has("edges")) {
        auto& collection = processing_ini["edges"];
        if (collection.has("cemax_fraction"))
            ce_max = processing_ini.get("edges").get("cemax_fraction");
    }
    if (stoi(cetypes_number_string) > 0) max_cfractions_vectors.at(1).push_back(stod(ce_max)); // 2 - faces

/// Nodes
//processing_mode
    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("pn_mode"))
        {
            ptype_vector.at(0) = processing_ini.get("nodes").get("pn_mode");
        } }

    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("n_multiplexity"))
            pindex_vector.at(0) = stod(processing_ini.get("nodes").get("n_multiplexity"));
        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode
    }

    string nseq_source;
    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("source"))
            nseq_source = processing_ini.get("nodes").get("source");
    }

    string ntypes_number_string;
    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("node_types_number"))
            ntypes_number_string = processing_ini.get("nodes").get("node_types_number");
    }

// fractions
    string n1_max, n2_max, n3_max;
    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("nmax_fraction1"))
            n1_max = processing_ini.get("nodes").get("nmax_fraction1");
    }
    if (stoi(ntypes_number_string) > 0) max_fractions_vectors.at(0).push_back(stod(n1_max)); // 0 - nodes

    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("nmax_fraction2"))
            n2_max = processing_ini.get("nodes").get("nmax_fraction2");
    }
    if (stoi(ntypes_number_string) > 0) max_fractions_vectors.at(0).push_back(stod(n2_max)); // 0 - nodes

    if (processing_ini.has("nodes")) {
        auto& collection = processing_ini["nodes"];
        if (collection.has("nmax_fraction3"))
            n3_max = processing_ini.get("nodes").get("nmax_fraction3");
    }
    if (stoi(ntypes_number_string) > 0) max_fractions_vectors.at(0).push_back(stod(n3_max)); // 0 - nodes

// III: distribution
    if (processing_ini.has("distribution")) {
        auto& collection = processing_ini["distribution"];
        if (collection.has("mu"))
            mu = stod(processing_ini.get("distribution").get("mu"));
    }
    if (processing_ini.has("distribution")) {
        auto& collection = processing_ini["distribution"];
        if (collection.has("sigma"))
            sigma = stod(processing_ini.get("distribution").get("sigma"));
    }
    if (processing_ini.has("distribution")) {
        auto& collection = processing_ini["distribution"];
        if (collection.has("bins_number"))
            bins_numb = stod(processing_ini.get("distribution").get("bins_number"));
    }

    /// sequences
    sequence_source_paths = {nseq_source, eseq_source, fseq_source, pseq_source};

vector<double> max_fractions_output(3, 0); // temporary vector serving as an output template for max fractions
/// Output to the screen/console
    cout << "The Processing module simulation type and initial parameters:\t\t" << endl;
    cout << endl;
    if (ptypes_number_string != "0") {
        // polyhedrons
        cout << "Processing p_type:\t"s << ptype_vector.at(3) << "\t with p_index:\t"s << pindex_vector.at(3) << endl;
        if (ptype_vector.at(3) == "L") cout << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(3) == "S") cout << "polyhedron sequence source: "s << pseq_source << endl;
        cout << "Number of polyhedron types:\t"s << ptypes_number_string << endl;
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[3].size() > 0 && max_fractions_vectors[3][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[3][i];
        cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        cout << endl;
    }
    if (ftypes_number_string != "0") {
        // faces
        cout << "Processing f_type:\t"s << ptype_vector.at(2) << "\t with f_index:\t"s << pindex_vector.at(2) << endl;
        if (ptype_vector.at(2) == "L") cout << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(2) == "S") cout << "face sequence source: "s << fseq_source << endl;
        cout << "Number of face types:\t"s << ftypes_number_string << endl;
// refill 0s
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[2].size() > 0 && max_fractions_vectors[2][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[2][i];
        cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        cout << endl;
    }

    if (cftypes_number_string != "0") {
        cout << "Processing c_face_type:\t"s << ctype_vector.at(2) << endl;
        cout << "Number of face crack types:\t"s << cftypes_number_string << endl;
        if (ctype_vector.at(2) == "Km") cout << "Their maximum fractions:\t"s << max_cfractions_vectors[2][0] << endl;
    }

    if (etypes_number_string != "0") {
        //edges
        cout << "Processing e_type:\t"s << ptype_vector.at(1) << "\twith e_index:\t"s << pindex_vector.at(1) << endl;
        if (ptype_vector.at(1) == "L") cout << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(1) == "S") cout << "edges sequence source: "s << eseq_source << endl;
        cout << "Number of edge types:\t"s << etypes_number_string << endl;
// refill 0s
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[1].size() > 0 && max_fractions_vectors[1][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[1][i];
        cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        cout << endl;
    }

    if (cetypes_number_string != "0") {
        cout << "Processing c_edge_type:\t"s << ctype_vector.at(1) << endl;
        cout << "Number of edge crack types:\t"s << cetypes_number_string << endl;
        if (ctype_vector.at(1) == "Km") cout << "Their maximum fractions:\t"s << max_cfractions_vectors[1][0] << endl;
    }

    if (ntypes_number_string != "0") {
        // nodes
        cout << "Processing n_type:\t"s << ptype_vector.at(0) << "\twith n_index:\t"s << pindex_vector.at(0) << endl;
        if (ptype_vector.at(0) == "L") cout << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(0) == "S") cout << "nodes sequence source: "s << nseq_source << endl;
        cout << "Number of node types:\t"s << ntypes_number_string << endl;
    }
    // refill 0s
    std::fill(max_fractions_output.begin(), max_fractions_output.end(),0);
    for (int i = 0; i < 3; ++i)
        if (max_fractions_vectors[0].size() > 0 && max_fractions_vectors[0][i] > 0) max_fractions_output.at(i) = max_fractions_vectors[0][i];
    cout << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1) << "\t\t"<< max_fractions_output.at(2) << "\t\t" << endl;
    cout<< "_________________________________________________" << endl << endl;

/// Output into .log file
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function

    Out_logfile_stream << "The Processing module simulation type and initial parameters:\t\t" << endl;
    Out_logfile_stream << endl;
    if (ptypes_number_string != "0") {
        // polyhedrons
        Out_logfile_stream << "Processing p_type:\t"s << ptype_vector.at(3) << "\t with p_index:\t"s << pindex_vector.at(3) << endl;
        if (ptype_vector.at(3) == "L") Out_logfile_stream << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(3) == "S") Out_logfile_stream << "polyhedron sequence source: "s << pseq_source << endl;
        Out_logfile_stream << "Number of polyhedron types:\t"s << ptypes_number_string << endl;
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[3].size() > 0 && max_fractions_vectors[3][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[3][i];
        Out_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    }
    if (ftypes_number_string != "0") {
        // faces
        Out_logfile_stream << "Processing f_type:\t"s << ptype_vector.at(2) << "\t with f_index:\t"s << pindex_vector.at(2) << endl;
        if (ptype_vector.at(2) == "L") Out_logfile_stream << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(2) == "S") Out_logfile_stream << "face sequence source: "s << fseq_source << endl;
        Out_logfile_stream << "Number of face types:\t"s << ftypes_number_string << endl;
// refill 0s
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[2].size() > 0 && max_fractions_vectors[2][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[2][i];
        Out_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    }
    if (cftypes_number_string != "0") {
        Out_logfile_stream << "Processing c_face_type:\t"s << ctype_vector.at(2) << endl;
        Out_logfile_stream << "Number of face crack types:\t"s << cftypes_number_string << endl;
        if (ctype_vector.at(2) == "Km") Out_logfile_stream << "Their maximum fractions:\t"s << max_cfractions_vectors[2][0] << endl;
    }

    //edges
    if (etypes_number_string != "0") {
        Out_logfile_stream << "Processing e_type:\t"s << ptype_vector.at(1) << "\twith e_index:\t"s << pindex_vector.at(1) << endl;
        if (ptype_vector.at(1) == "L") Out_logfile_stream << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(1) == "S") Out_logfile_stream << "edges sequence source: "s << eseq_source << endl;
        Out_logfile_stream << "Number of edge types:\t"s << etypes_number_string << endl;
// refill 0s
        std::fill(max_fractions_output.begin(), max_fractions_output.end(), 0);
        for (int i = 0; i < 3; ++i)
            if (max_fractions_vectors[1].size() > 0 && max_fractions_vectors[1][i] > 0)
                max_fractions_output.at(i) = max_fractions_vectors[1][i];
        Out_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1)
             << "\t\t" << max_fractions_output.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    }

    if (cetypes_number_string != "0") {
        Out_logfile_stream << "Processing c_edge_type:\t"s << ctype_vector.at(1) << endl;
        Out_logfile_stream << "Number of edge crack types:\t"s << cetypes_number_string << endl;
        if (ctype_vector.at(1) == "Km") Out_logfile_stream << "Their maximum fractions:\t"s << max_cfractions_vectors[1][0] << endl;
    }

    if (ntypes_number_string != "0") {
        // nodes
        Out_logfile_stream << "Processing n_type:\t"s << ptype_vector.at(0) << "\twith n_index:\t"s << pindex_vector.at(0) << endl;
        if (ptype_vector.at(0) == "L") Out_logfile_stream << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
        if (ptype_vector.at(0) == "S") Out_logfile_stream << "nodes sequence source: "s << nseq_source << endl;
        Out_logfile_stream << "Number of node types:\t"s << ntypes_number_string << endl;
    }
// refill 0s
    std::fill(max_fractions_output.begin(), max_fractions_output.end(),0);
    for (int i = 0; i < 3; ++i)
        if (max_fractions_vectors[1].size() > 0 && max_fractions_vectors[1][i] > 0) max_fractions_output.at(i) = max_fractions_vectors[1][i];
    Out_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1) << "\t\t"<< max_fractions_output.at(2) << "\t\t" << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Processing n_type:\t"s << ptype_vector.at(0) << "\twith n_index:\t"s << pindex_vector.at(0) << endl;
    if (ptype_vector.at(0) == "L") Out_logfile_stream << "mu = \t"s << mu << " and " << "sigma = \t"s << sigma << endl;
    Out_logfile_stream << "Number of node types:\t"s << ntypes_number_string << endl;
// refill 0s
    std::fill(max_fractions_output.begin(), max_fractions_output.end(),0);
    for (int i = 0; i < 3; ++i)
        if (max_fractions_vectors[0].size() > 0 && max_fractions_vectors[0][i] > 0) max_fractions_output.at(i) = max_fractions_vectors[0][i];
    Out_logfile_stream << "Their maximum fractions:\t"s << max_fractions_output.at(0) << "\t\t" << max_fractions_output.at(1) << "\t\t"<< max_fractions_output.at(2) << "\t\t" << endl;
    Out_logfile_stream<< "_________________________________________________" << endl << endl;

    Out_logfile_stream.close(); // this *.log stream will be closed at the end of the main function

    return;
} /// END of config_reader_processing function


/// ================== # 3 # Initial CHARACTERISATION module configuration - reading and output ==================
std::vector<double> config_reader_characterisation(std::string const &source_path, std::vector<int> &charlabs_polyhedrons, std::vector<int> &charlabs_faces, std::vector<int> &charlabs_edges, std::vector<int> &charlabs_nodes, std::vector<int> &charlabs_laplacians, std::ofstream &Out_logfile_stream) {
    std::vector<double> config_characterisation_vector;

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

    /// Console output
//    cout<< "______________________________________________________________________________________" << endl;
    cout << "The Characterisation module simulation type and initial parameters:\t\t" << endl;
    cout << endl;
    cout << "Polyhedrons lab ON/OFF:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
    cout << "Faces lab       ON/OFF:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
    cout << "Edges lab       ON/OFF:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
    cout << "Nodes lab       ON/OFF:\t"s << charlabs_nodes.at(0) << "\t\t" << endl;
    cout << endl;

    if(charlabs_polyhedrons.at(0) == 1) { // Polyhedrons
        cout << "Polyhedrons configuration entropy:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
        cout << "Conf entropy mean part:\t"s << charlabs_polyhedrons.at(1) << "\t\t" << endl;
        cout << "Conf entropy skew part:\t"s << charlabs_polyhedrons.at(2) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_polyhedrons.at(0) == 1)
    if(charlabs_faces.at(0) == 1) { // Faces
        cout << "Faces configuration entropy:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
        cout << "Conf entropy mean part:\t"s << charlabs_faces.at(1) << "\t\t" << endl;
        cout << "Conf entropy skew part:\t"s << charlabs_faces.at(2) << "\t\t" << endl;
        cout << "Edges fractions:\t"s << charlabs_faces.at(3) << "\t\t" << endl;
        cout << "Edges degree fractions:\t"s << charlabs_faces.at(4) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_faces.at(0) == 1)
    if(charlabs_edges.at(0) == 1) { // Edges
        cout << "Edges configuration entropy:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
        cout << "Conf entropy Mean (-) Skew:\t"s << charlabs_edges.at(1) << "\t\t" << endl;
        cout << "Conf entropy mean part:\t"s << charlabs_edges.at(2) << "\t\t" << endl;
        cout << "Conf entropy skew part:\t"s << charlabs_edges.at(3) << "\t\t" << endl;
        cout << "Analytical solutions :\t"s << charlabs_edges.at(4) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_edges.at(0) == 1)
    if(charlabs_nodes.at(0) == 1) { // Nodes
        cout << "Node configuration entropy:\t"s << charlabs_nodes.at(0) << "\t\t" << endl;
        cout << "Conf entropy mean part:\t"s << charlabs_nodes.at(1) << "\t\t" << endl;
        cout << "Conf entropy skew part:\t"s << charlabs_nodes.at(2) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_polyhedrons.at(0) == 1)

    if(charlabs_laplacians.at(0) > 0) { // Laplacians lab
        cout << "Laplacians: number of calculation steps \t"s << charlabs_laplacians.at(0) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_laplacians.at(0) == 1)

    if(charlabs_laplacians.at(1) == 1) { // Laplacians lab
        cout << "Special cell Laplacians:\t"s << charlabs_laplacians.at(1) << "\t\t" << endl;
        cout << endl;
    } // if(charlabs_laplacians.at(0) == 1)
//    cout<< "______________________________________________________________________________________" << endl;

/// Output into Processing_Design.log file
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function

//    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
    Out_logfile_stream << "The Characterisation module simulation type and initial parameters:\t\t" << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Polyhedrons lab ON/OFF:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
    Out_logfile_stream << "Faces lab       ON/OFF:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
    Out_logfile_stream << "Edges lab       ON/OFF:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
    Out_logfile_stream << "Nodes lab       ON/OFF:\t"s << charlabs_nodes.at(0) << "\t\t" << endl;
    Out_logfile_stream << endl;

    if(charlabs_polyhedrons.at(0) == 1) { // Polyhedrons
        Out_logfile_stream << "Polyhedrons configuration entropy:\t"s << charlabs_polyhedrons.at(0) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy mean part:\t"s << charlabs_polyhedrons.at(1) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy skew part:\t"s << charlabs_polyhedrons.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_polyhedrons.at(0) == 1)
    if(charlabs_faces.at(0) == 1) { // Faces
        Out_logfile_stream << "Faces configuration entropy:\t"s << charlabs_faces.at(0) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy mean part:\t"s << charlabs_faces.at(1) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy skew part:\t"s << charlabs_faces.at(2) << "\t\t" << endl;
        Out_logfile_stream << "Edges fractions:\t"s << charlabs_faces.at(3) << "\t\t" << endl;
        Out_logfile_stream << "Edges degree fractions:\t"s << charlabs_faces.at(4) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_faces.at(0) == 1)
    if(charlabs_edges.at(0) == 1) { // Edges
        Out_logfile_stream << "Edges configuration entropy:\t"s << charlabs_edges.at(0) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy mean part:\t"s << charlabs_edges.at(1) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy skew part:\t"s << charlabs_edges.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_edges.at(0) == 1)
    if(charlabs_nodes.at(0) == 1) { // Nodes
        Out_logfile_stream << "Node configuration entropy:\t"s << charlabs_nodes.at(0) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy mean part:\t"s << charlabs_nodes.at(1) << "\t\t" << endl;
        Out_logfile_stream << "Conf entropy skew part:\t"s << charlabs_nodes.at(2) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_polyhedrons.at(0) == 1)

    if(charlabs_laplacians.at(0) > 0) { // Laplacians lab
        Out_logfile_stream << "Laplacians: number of calculation steps \t"s << charlabs_laplacians.at(0) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_laplacians.at(0) == 1)

    if(charlabs_laplacians.at(1) == 1) { // Laplacians lab
        Out_logfile_stream << "Special cell Laplacians:\t"s << charlabs_laplacians.at(1) << "\t\t" << endl;
        Out_logfile_stream << endl;
    } // if(charlabs_laplacians.at(0) == 1)
//    Out_logfile_stream<< "______________________________________________________________________________________" << endl;

    Out_logfile_stream.close(); // this *.log stream will be closed at the end of the main function

    return config_characterisation_vector;
} // END of config characterisation reader function


/// ================== # 4 # Initial WRITER module configuration - reading and output ==================
void config_reader_writer(std::string &source_path, std::vector<int> &writer_specifications, std::ofstream &Out_logfile_stream) {
/// writer_specifications vector ::
int    isSequencesOutput = 0;      // - >     [0]
int    isDesignvectorsOutput = 0;  // - >     [1]
int    isEnergiesOutput = 0;       // - >     [2]
int isEdgeConfEntropy = 0, isEdgeFractions = 0, isDegreeEdgeFractions = 0, isEdgeAnFractions = 0, isEdgeAnConfEntropies = 0; // [2], [3], [4], [5], [6]
int isBetti = 0; // Laplacians lab

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


/// Output to the screen/console
//    cout << endl;
//    cout<< "______________________________________________________________________________________" << endl;
    cout << "The Writer module specifications:\t\t" << endl;
//    cout << endl;
    cout << "Sequences output \t\t\t\t\t"s << writer_specifications.at(0) << endl;
    cout << "Design vectors output \t\t\t\t"s << writer_specifications.at(1) << endl;
    cout << "Configuration Edges entropy \t\t"s << writer_specifications.at(2) << endl;
    cout << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(3) << endl;
    cout << "Special Edge degree fractions \t\t"s << writer_specifications.at(4) << endl;
    cout << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
    cout << "Analytical Edge degree fractions \t"s << writer_specifications.at(5) << endl;
    cout << "Analytical Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
    cout << "Laplacians and Betti numbers \t\t"s << writer_specifications.at(7) << endl;
    cout << "Cell Energies \t\t\t\t\t\t"s << writer_specifications.at(8) << endl;
    cout << endl;

/// Output into .log file
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function

//    Out_logfile_stream << endl;
//    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
    Out_logfile_stream << "The Writer module specifications:\t\t" << endl;
//    Out_logfile_stream << endl;
    Out_logfile_stream << "Sequences output \t\t\t\t\t"s << writer_specifications.at(0) << endl;
    Out_logfile_stream << "Design vectors output \t\t\t\t"s << writer_specifications.at(1) << endl;
    Out_logfile_stream << "Configuration Edges entropy \t\t\t"s << writer_specifications.at(2) << endl;
    Out_logfile_stream << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(3) << endl;
    Out_logfile_stream << "Special Edge fractions \t\t\t\t"s << writer_specifications.at(4) << endl;
    Out_logfile_stream << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
    Out_logfile_stream << "Analytical Edge degree fractions \t"s << writer_specifications.at(5) << endl;
    Out_logfile_stream << "Analytical Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
    Out_logfile_stream << "Laplacians and Betti numbers \t\t"s << writer_specifications.at(7) << endl;
    Out_logfile_stream << "Cell Energies \t\t\t\t\t\t"s << writer_specifications.at(8) << endl;
    Out_logfile_stream << endl;

    Out_logfile_stream.close(); // this *.log stream will be closed at the end of the main function

    return;
} /// END of config_reader_writer function

/// ================== # 5 # Initial SUBCOMPLEX module configuration - reading and output ==================
void config_reader_subcomplex(std::string &source_path, std::string &sctype, std::vector<double> &plane_orientation, double &cut_length, unsigned int &grain_neighbour_orders, std::ofstream &Out_logfile_stream) {

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "subcomplex.ini"s);
    mINI::INIStructure subcomplex_ini;
    file.read(subcomplex_ini);

    cut_length;
//subcomplex type
    if (subcomplex_ini.has("subcomplex_type")) {
        auto& collection = subcomplex_ini["subcomplex_type"];
        if (collection.has("subPCC_type"))
        {
            sctype = subcomplex_ini.get("subcomplex_type").get("subPCC_type");
        } }

    //plane orientation
    if (subcomplex_ini.has("plane_section")) {
        auto& collection = subcomplex_ini["plane_section"];
        if (collection.has("a_coeff")) {
            plane_orientation.at(0) = stod(subcomplex_ini.get("plane_section").get("a_coeff"));
        }
        if (collection.has("b_coeff")) {
            plane_orientation.at(1) = stod(subcomplex_ini.get("plane_section").get("b_coeff"));
        }
        if (collection.has("c_coeff")) {
            plane_orientation.at(2) = stod(subcomplex_ini.get("plane_section").get("c_coeff"));
        }
        if (collection.has("D_coeff")) {
            plane_orientation.at(3) = stod(subcomplex_ini.get("plane_section").get("D_coeff"));
        }
    } // end of if (subcomplex_ini.has("plane_section"))

    //half-plane cut length
    if (subcomplex_ini.has("half_plane_section")) {
        auto& collection = subcomplex_ini["half_plane_section"];
        if (collection.has("half_plane_length"))
        {
            cut_length = stod(subcomplex_ini.get("half_plane_section").get("half_plane_length"));
        } }

    // k_order_neighbours
    if (subcomplex_ini.has("k_order_neighbours")) {
        auto& collection = subcomplex_ini["k_order_neighbours"];
        if (collection.has("neighbours_order"))
        {
            grain_neighbour_orders = stoi(subcomplex_ini.get("k_order_neighbours").get("neighbours_order"));
        } }

    /// Output to the screen/console
    cout << "The Subcomplex module type and initial parameters:\t\t" << endl << endl;
    cout << "Subcomplex type:\t"s << sctype << endl;
    if (sctype == "H"s) {
        cout << "Half-plane length:\t"s << cut_length << endl << endl;
    }
    if (sctype == "P"s || sctype == "H"s) {
        cout << "Plane orientation:\ta_coeff*X + b_coeff*Y + c_coeff*Z = D"s << endl << "Plane normal vector\t"s << "\ta =\t" << plane_orientation.at(0) << "\tb =\t" << plane_orientation.at(1) << "\tc =\t" << plane_orientation.at(2) << "\t\tPlane position D =\t"s << plane_orientation.at(3) << endl;
    }
    else if(sctype == "N"s){
        cout << "Grain k-neighbours order:\t"s << grain_neighbour_orders << endl << endl;
    }
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function

    Out_logfile_stream << "The Subcomplex module type and initial parameters:\t\t" << endl << endl;
    Out_logfile_stream << "Subcomplex type:\t"s << sctype << endl;
    if (sctype == "H"s) {
        Out_logfile_stream << "Half-plane length:\t"s << cut_length << endl << endl;
    }
    if (sctype == "P"s || sctype == "H"s) {
        Out_logfile_stream << "Plane orientation:\ta_coeff*X + b_coeff*Y + c_coeff*Z = D"s << endl << "Plane normal vector:\t"s << " a: " << plane_orientation.at(0) << " b: " << plane_orientation.at(1) << " c: " << plane_orientation.at(2) << "\t\tPlane position:\t" << " D: " << plane_orientation.at(3) << endl;
    }
    else if(sctype == "N"s){
        Out_logfile_stream << "Grain k-neighbours order:\t"s << grain_neighbour_orders << endl << endl;
    }

    Out_logfile_stream.close();

    return;
} /// end of the bool SubcomplexON() function


/// ================== # 6 # Initial MULTIPFYSICS module configuration - physical dimesions and all ==================
void config_reader_multiphysics(std::string &source_path, std::string &Mid_matrix, std::string &Mid_inclusion1, std::tuple<double, double, double> &sample_dimensions, double &tau, Eigen::MatrixXd &ext_stress_tensor, std::vector<double> &macrocrack_ini, std::ofstream &Out_logfile_stream) {

    double lx_size = 0.0, ly_size = 0.0, lz_size = 0.0; // sample dimensions
    double sxx = 0.0, sxy = 0.0, sxz = 0.0, syx = 0.0, syy = 0.0, syz = 0.0, szx = 0.0, szy = 0.0, szz = 0.0; // external stress tensor components [homogeneous stress state]

// ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "multiphysics.ini"s);
    mINI::INIStructure multiphysics_ini;
    file.read(multiphysics_ini);

// I
// Material ID for the CPD code Database
    if (multiphysics_ini.has("material_id")) {
        auto &collection = multiphysics_ini["material_id"];
        if (collection.has("Mid_matrix"))
            Mid_matrix = multiphysics_ini.get("material_id").get("Mid_matrix");
    }
// Material ID for the CPD code Database
    if (multiphysics_ini.has("material_id")) {
        auto &collection = multiphysics_ini["material_id"];
        if (collection.has("Mid_inclusion1"))
            Mid_inclusion1 = multiphysics_ini.get("material_id").get("Mid_inclusion1");
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

    std::get<0>(sample_dimensions) = lx_size;
    std::get<1>(sample_dimensions) = ly_size;
    std::get<2>(sample_dimensions) = lz_size;

// III
    if (multiphysics_ini.has("time_scale")) {
        auto &collection = multiphysics_ini["time_scale"];
        if (collection.has("tau"))
            tau = stod(multiphysics_ini.get("time_scale").get("tau"));
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
    } // end of  if (multiphysics_ini.has("physical_dimensions"))

    double vonMises_stress = 0.0, pressure = 0.0;
    if (sxx > 1000000.0 || syy > 1000000.0 || szz > 1000000.0 || sxy > 1000000.0 || sxz > 1000000.0 || syx > 1000000.0 || syz > 1000000.0) {
        cout << endl
             << "WARNING: Highly likely stress values in the '../config/multiphysics.ini' file are not in [MPa] as requested !!!"
             << endl << endl;
        Out_logfile_stream << endl
                           << "WARNING: Highly likely stress values in the '../config/multiphysics.ini' file are not in [MPa] as requested !!!"
                           << endl << endl;
    }
        ext_stress_tensor(0,0) = sxx*pow(10,6);
        ext_stress_tensor(0,1) = sxy*pow(10,6);
        ext_stress_tensor(0,2) = sxz*pow(10,6);
        ext_stress_tensor(1,0) = syx*pow(10,6);
        ext_stress_tensor(1,1) = syy*pow(10,6);
        ext_stress_tensor(1,2) = syz*pow(10,6);
        ext_stress_tensor(2,0) = szx*pow(10,6);
        ext_stress_tensor(2,1) = szy*pow(10,6);
        ext_stress_tensor(2,2) = szz*pow(10,6);
        vonMises_stress = std::sqrt(0.5 * (pow((ext_stress_tensor(0,0) - ext_stress_tensor(1,1)), 2.0) + pow((ext_stress_tensor(0,0) - ext_stress_tensor(2,2)), 2.0) + pow((ext_stress_tensor(1,1) - ext_stress_tensor(2,2)), 2.0) ) + 3.0 * ( pow(ext_stress_tensor(0,1), 2.0) + pow(ext_stress_tensor(1,2), 2.0) + pow(ext_stress_tensor(2,0), 2.0) ));
        pressure = (ext_stress_tensor(0,0) + ext_stress_tensor(1,1) + ext_stress_tensor(2,2)) / 3.0;


macrocrack_ini.resize(5,0);
double number_of_crack_sizes = 0.0, max_crack_lenghts = 0.0, min_crack_lenghts = 0.0, crack_stress_mode = 0.0, macrocrack_number = 0.0;
    if (multiphysics_ini.has("macrocracks")) {
        auto &collection = multiphysics_ini["macrocracks"];
        if (collection.has("number_of_macrocracks"))
            macrocrack_number = stod(multiphysics_ini.get("macrocracks").get("number_of_macrocracks"));
    }
    if(macrocrack_number > 0){
        if (multiphysics_ini.has("macrocracks")) {
            auto &collection = multiphysics_ini["macrocracks"];
            if (collection.has("crack_stress_mode"))
                crack_stress_mode = stod(multiphysics_ini.get("macrocracks").get("crack_stress_mode"));
        }
        if (multiphysics_ini.has("macrocracks")) {
            auto &collection = multiphysics_ini["macrocracks"];
            if (collection.has("min_crack_lenghts"))
                min_crack_lenghts = stod(multiphysics_ini.get("macrocracks").get("min_crack_lenghts"));
        }
        if (multiphysics_ini.has("macrocracks")) {
            auto &collection = multiphysics_ini["macrocracks"];
            if (collection.has("max_crack_lenghts"))
                max_crack_lenghts = stod(multiphysics_ini.get("macrocracks").get("max_crack_lenghts"));
        }
        if (multiphysics_ini.has("macrocracks")) {
            auto &collection = multiphysics_ini["macrocracks"];
            if (collection.has("number_of_crack_sizes"))
                number_of_crack_sizes = stod(multiphysics_ini.get("macrocracks").get("number_of_crack_sizes"));
        }
    } // end of if(macrocrack_number > 0)

    macrocrack_ini.at(0) = macrocrack_number; macrocrack_ini.at(1) = crack_stress_mode; macrocrack_ini.at(2) = max_crack_lenghts;
    macrocrack_ini.at(3) = min_crack_lenghts; macrocrack_ini.at(4) = number_of_crack_sizes;
/// Output to the screen/console
    cout << "______________________________________________________________________________________" << endl;
    cout << "The Multiphysics module specifications:\t\t" << endl;
    cout << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(sample_dimensions) << " [m] "s << ", y: "
         << std::get<1>(sample_dimensions) << " [m] "s << ", z: " << std::get<2>(sample_dimensions) << " [m] "s << endl;
    cout << "Characteristic time is \t\t\t"s << " tau: " << tau*pow(10,6) << " [microseconds] "s << endl;

    // Homogeneous External Stress State
    cout << "Pressure is equal to \t\t\t"s << " P: "<< pressure/pow(10,6) << " [MPa] "s  << endl;
    cout << "Von Mises stress is equal to \t"s << " Sv: " << vonMises_stress/pow(10,6) << " [MPa] "s << endl;

    if (pressure || vonMises_stress > 0) {
        cout << "External Stress [MPa]: "s << endl;
        cout << ext_stress_tensor << endl;
    }
    cout << endl;
    if(macrocrack_number > 0) {
        cout << "Number of macrocracks  \t\t\t\t"s << macrocrack_ini.at(0) << endl;
        cout << "Crack mode  \t\t\t\t\t\t"s << macrocrack_ini.at(1) << endl;
        cout << "MAX crack lenghts (fraction)  \t\t"s << macrocrack_ini.at(2) << endl;
        cout << "MIN crack lenghts (fraction)  \t\t"s << macrocrack_ini.at(3) << endl;
        cout << "Series of crack sizes (number)  \t"s << macrocrack_ini.at(4) << endl;
    }
    cout << endl;
/// Output into .log file
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this *.log stream will be closed at the end of the main function

    Out_logfile_stream << "______________________________________________________________________________________"
                       << endl;
    Out_logfile_stream << "The Multiphysics module specifications:\t\t" << endl;
    Out_logfile_stream << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(sample_dimensions) << " [m] "s << ", y: "
         << std::get<1>(sample_dimensions) << " [m] "s << ", z: " << std::get<2>(sample_dimensions) << " [m] "s << endl;
    Out_logfile_stream << "Characteristic time is \t\t\t"s << " tau: " << tau*pow(10,6) << " [microseconds] "s << endl;

    // Homogeneous External Stress State
    Out_logfile_stream << "Pressure is equal to \t\t\t"s << " P: "<< pressure/pow(10,6) << " [MPa] "s  << endl;
    Out_logfile_stream << "Von Mises stress is equal to \t"s << " Sv: " << vonMises_stress/pow(10,6) << " [MPa] "s << endl;

    if (pressure || vonMises_stress > 0) {
        Out_logfile_stream << "External Stress [MPa]: "s << endl;
        Out_logfile_stream << ext_stress_tensor << endl;
    }
    Out_logfile_stream << endl;
    if(macrocrack_number > 0) {
        Out_logfile_stream << "Number of macrocracks  \t\t\t\t"s << macrocrack_ini.at(0) << endl;
        Out_logfile_stream << "Crack mode  \t\t\t\t\t\t"s << macrocrack_ini.at(1) << endl;
        Out_logfile_stream << "MAX crack lenghts (fraction)  \t\t"s << macrocrack_ini.at(2) << endl;
        Out_logfile_stream << "MIN crack lenghts (fraction)  \t\t"s << macrocrack_ini.at(3) << endl;
        Out_logfile_stream << "Series of crack sizes (number)  \t"s << macrocrack_ini.at(4) << endl;
    }
    Out_logfile_stream << endl;
    Out_logfile_stream.close();

    return;
} /// end of the void config_reader_multiphysics() function