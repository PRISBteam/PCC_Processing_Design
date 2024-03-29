/// Author: Dr Elijah Borodin (2023)
/// Manchester, UK
/// Library of specific functions related to the PCC Processing Design code for reading its *.ini files
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Simple reader for *.ini files and a specific CPD code-related library for reading its particular *.ini files ( downloaded from https://github.com/pulzed/mINI )
#include "../ini/ini.h"

using namespace std; // standard namespace

/// ================== # 1 # Initial configuration - reading and output ==================
std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string &main_type, std::string &e_type) {

    std::vector<int> res(7,0);
/// [0] - > dim, [1] -> isSection, [2] -> isProcessing, [3] -> isCharacterisation, [4] -> isMultiphysics, [5] -> isKinetic, [6] -> isWriter
    bool isSectionON = 0, isProcessingON = 0, isCharacterisationON = 0, isKineticON = 0, isMultiphysicsON = 0, isWriterON = 0;
    std::string isSection, isProcessing, isCharacterisation, isKinetic, isMultiphysics, isWriter;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "main.ini"s);
//    mINI::INIFile file(source_path + "main_2D.ini"s);

    mINI::INIStructure main_ini;
    file.read(main_ini);

// 0
    if (main_ini.has("execution_type")) {
        auto& collection = main_ini["execution_type"];
        if (collection.has("e_type"))
        {
            e_type = main_ini.get("execution_type").get("e_type");
        } }
// I
    if (main_ini.has("simulation_mode")) {
        auto& collection = main_ini["simulation_mode"];
        if (collection.has("mode"))
        {
            main_type = main_ini.get("simulation_mode").get("mode");
        } }
// II
    std::string problem_dimension;
    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("dim"))
            problem_dimension = main_ini.get("general").get("dim");
    }
    res.at(0) = stoi(problem_dimension); // res[0]

    // III
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Section"))
            isSection = main_ini.get("modules").get("PCC_Section");
    }
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Processing"))
            isProcessing = main_ini.get("modules").get("PCC_Processing");
    }
/*
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Characterisation"))
            isCharacterisation = main_ini.get("modules").get("PCC_Characterisation");
    }
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Multiphysics"))
            isMultiphysics = main_ini.get("modules").get("PCC_Multiphysics");
    }
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Kinetic"))
            isKinetic = main_ini.get("modules").get("PCC_Kinetic");
    }
*/
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Writer"))
            isWriter = main_ini.get("modules").get("PCC_Writer");
    }

    /// forming the output RES vector
// ON/OFF IDs
    if (isSection == "ON") { isSectionON = 1; res.at(1) = 1; } else res.at(1) = 0; // res[1] - Section
    if (isProcessing == "ON") { isProcessingON = 1; res.at(2) = 1; } else res.at(2) = 0; // res[2] - Processing
//    if (isCharacterisation == "ON") { isCharacterisationON = 1; res.at(3) = 1; } else res.at(3) = 0; // res[3] - Characterisation
//    if (isMultiphysics == "ON") { isMultiphysicsON = 1; res.at(4) = 1; } else res.at(4) = 0; // res[4] - Multiphysics
//    if (isKinetic == "ON") { isKineticON = 1; res.at(5) = 1; } else res.at(5) = 0; // res[5] - Kinetic
    if (isWriter == "ON") { isWriterON = 1; res.at(6) = 1; } else res.at(6) = 0; // res[6] - Writer

    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("source_dir"))
            source_dir = main_ini.get("general").get("source_dir");
    }
    if (main_ini.has("general")) {
        auto& collection = main_ini["general"];
        if (collection.has("output_dir"))
            output_dir = main_ini.get("general").get("output_dir");
    }

/// Output to the screen/console
    cout << "The problem dimension that is the maximum value k_max of k-cells in the PCC\t\t|\t\t"s << "dim = " << res.at(0) << endl;
    cout << "Execution type:\t\t"s << e_type << endl;
    cout << "Simulation type:\t"s << main_type << endl;
    cout << "Source directory:\t"s << source_dir << endl;
    cout << "Output directory:\t"s << output_dir << endl;
    cout << endl;
    if (isSectionON == 1) cout << "ON    | PCC_Section"s << endl;
    else cout << "OFF    | PCC_Section"s << endl;
    if (isProcessingON == 1) cout << "ON    | PCC_Processing"s << endl;
    else cout << "OFF    | PCC_Processing"s << endl;
//    if (isCharacterisationON == 1) cout << "ON    | PCC_Characterisation"s << endl;
//    else cout << "OFF    | PCC_Characterisation"s << endl;
//    if (isMultiphysicsON == 1) cout << "ON    | PCC_Multiphysics"s << endl;
//    else cout << "OFF    | PCC_Multiphysics"s << endl;
//    if (isKineticON == 1) cout << "ON    | PCC_Kinetic"s << endl;
//    else cout << "OFF    | PCC_Kinetic"s << endl;
    if (isWriterON == 1) cout << "ON    | PCC_Writer"s << endl;
    else cout << "OFF    | PCC_Writer"s << endl;
    cout << endl;

/// Output into .log file
    std::ofstream Out_logfile_stream;
    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::trunc); // this *.log stream will be closed at the end of the main function
    Out_logfile_stream << endl;
    Out_logfile_stream << "The problem dimension that is the maximum value k_max of k-cells in the PCC:\t\t|\t\t"s << "dim = " << res.at(0) << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Execution type:\t\t"s << e_type << endl;
    Out_logfile_stream << "Simulation type:\t"s << main_type << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Source directory:\t"s << source_dir << endl;
    Out_logfile_stream << "Output directory:\t"s << output_dir << endl;
    Out_logfile_stream << endl;
    if (isSectionON == 1) Out_logfile_stream << "ON    | PCC_Section"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Section"s << endl;
    if (isProcessingON == 1) Out_logfile_stream << "ON    | PCC_Processing"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Processing"s << endl;
//    if (isCharacterisationON == 1) Out_logfile_stream << "ON    | PCC_Characterisation"s << endl;
//    else Out_logfile_stream << "OFF    | PCC_Characterisation"s << endl;
//    if (isMultiphysicsON == 1) Out_logfile_stream << "ON    | PCC_Multiphysics"s << endl;
//    else Out_logfile_stream << "OFF    | PCC_Multiphysics"s << endl;
//    if (isKineticON == 1) Out_logfile_stream << "ON    | PCC_Kinetic"s << endl;
//    else Out_logfile_stream << "OFF    | PCC_Kinetic"s << endl;
    if (isWriterON == 1) Out_logfile_stream << "ON    | PCC_Writer"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Writer"s << endl;
    Out_logfile_stream << endl;

    Out_logfile_stream.close();

    return res;
} /// END of config_reader_main function

/// ================== # 2 # Initial PROCESSING module configuration - reading and output ==================
void config_reader_processing(std::string &source_path, std::vector<string> &sequence_source_paths, std::vector<vector<double>> &max_fractions_vectors, std::vector<vector<double>> &max_cfractions_vectors, double &mu, double &sigma, std::vector<string> &ptype_vector, std::vector<string> &ctype_vector, std::vector<double> &pindex_vector, std::ofstream &Out_logfile_stream) {
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
        if (collection.has("pp_index"))
            pindex_vector.at(3) = stod(processing_ini.get("polyhedrons").get("pp_index"));
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
        if (collection.has("pf_index"))
            pindex_vector.at(2) = stod(processing_ini.get("faces").get("pf_index"));
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
        if (collection.has("pe_index"))
            pindex_vector.at(1) = stod(processing_ini.get("edges").get("pe_index"));
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
        if (collection.has("pn_index"))
            pindex_vector.at(0) = stod(processing_ini.get("nodes").get("pn_index"));
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

    /// sequences
    sequence_source_paths = {nseq_source, eseq_source, fseq_source, pseq_source};

vector<double> max_fractions_output(3, 0); // temporary vector serving as an output template for max fractions
/// Output to the screen/console
    cout<< "______________________________________________________________________________________" << endl;
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
    cout<< "______________________________________________________________________________________" << endl;

/// Output into .log file
    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
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
    Out_logfile_stream<< "______________________________________________________________________________________" << endl;

    return;
} /// END of config_reader_processing function

/// ================== # 6 # Initial WRITER module configuration - reading and output ==================

void config_reader_writer(std::string &source_path, std::vector<int> &writer_specifications, std::ofstream &Out_logfile_stream) {
/// writer_specifications vector ::
int    isSequencesOutput;      // - >     [0]
int    isDesignvectorsOutput;  // - >     [1]
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

/// Output to the screen/console
//    cout << endl;
//    cout<< "______________________________________________________________________________________" << endl;
    cout << "The Writer module specifications:\t\t" << endl;
//    cout << endl;
    cout << "Sequences output \t\t\t"s << writer_specifications.at(0) << endl;
    cout << "Design vectors output \t\t"s << writer_specifications.at(1) << endl;
    cout << "Configuration Edges entropy \t\t\t"s << writer_specifications.at(2) << endl;
    cout << "Special Edge fractions \t\t\t"s << writer_specifications.at(3) << endl;
    cout << "Special Edge degree fractions \t\t\t"s << writer_specifications.at(4) << endl;
    cout << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
    cout << "Analytical Edge degree fractions \t\t\t"s << writer_specifications.at(5) << endl;
    cout << "Analytical configuration Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
    cout << "Laplacians and Betti numbers \t\t\t"s << writer_specifications.at(7) << endl;
    cout << endl;

/// Output into .log file
//    Out_logfile_stream << endl;
//    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
    Out_logfile_stream << "The Writer module specifications:\t\t" << endl;
//    Out_logfile_stream << endl;
    Out_logfile_stream << "Sequences output \t\t\t"s << writer_specifications.at(0) << endl;
    Out_logfile_stream << "Design vectors output \t\t"s << writer_specifications.at(1) << endl;
    Out_logfile_stream << "Configuration Edges entropy \t\t\t"s << writer_specifications.at(2) << endl;
    Out_logfile_stream << "Special Edge fractions \t\t\t"s << writer_specifications.at(3) << endl;
    Out_logfile_stream << "Special Edge degree fractions \t\t\t"s << writer_specifications.at(4) << endl;
    Out_logfile_stream << "Analytical Edge fractions \t\t\t"s << writer_specifications.at(5) << endl;
    Out_logfile_stream << "Analytical Edge degree fractions \t\t\t"s << writer_specifications.at(5) << endl;
    Out_logfile_stream << "Analytical configuration Edges entropy \t\t\t"s << writer_specifications.at(6) << endl;
    Out_logfile_stream << "Laplacians and Betti numbers \t\t\t"s << writer_specifications.at(7) << endl;
    Out_logfile_stream << endl;

    return;
} /// END of config_reader_writer function

/// ================== # 5 # Initial SUBCOMPLEX module configuration - reading and output ==================
void config_reader_subcomplex(std::string &source_path, std::string &sctype, double &cut_length, std::ofstream &Out_logfile_stream) {
/*    std::string line;
    ifstream inConf(config);
    bool isSectionON = 0;

    if (inConf) { // If the file was successfully open, then
        while(getline(inConf, line, '\n'))
// REPAIR            cout << line << endl;
            if (line.compare("DCC_Section SWITCHED ON"s)) {
                isSectionON = 1;
                if (time_step_one == 1) cout << "ON    | DCC_Section"s << endl;
                return isSectionON;
            }
    } else cout << "SubcomplexON() error: The file " << config << " cannot be read" << endl; // If something goes wrong

    if (time_step_one == 1) cout << "OFF   | DCC_Section"s << endl;
    return isSectionON;
    */
} /// end of the bool SubcomplexON() function

/// ================== # 6 # Initial MULTIPFYSICS module configuration - physical dimesions and all ==================
void config_reader_multiphysics(std::string &source_path, std::tuple<double, double, double> &sample_dimensions, std::ofstream &Out_logfile_stream) {

    double lx_size = 0.0, ly_size = 0.0, lz_size = 0.0; // sample dimensions
// ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "multiphysics.ini"s);
    mINI::INIStructure multiphysics_ini;
    file.read(multiphysics_ini);

// I
// sequences and designs output
    if (multiphysics_ini.has("physical_dimensions")) {
        auto& collection = multiphysics_ini["physical_dimensions"];
        if (collection.has("sample_dimension_x"))
            lx_size = stod(multiphysics_ini.get("physical_dimensions").get("sample_dimension_x"));
        if (collection.has("sample_dimension_y"))
            ly_size = stod(multiphysics_ini.get("physical_dimensions").get("sample_dimension_y"));
        if (collection.has("sample_dimension_z"))
            lz_size = stod(multiphysics_ini.get("physical_dimensions").get("sample_dimension_z"));
    } // end of  if (multiphysics_ini.has("physical_dimensions"))

    std::get<0>(sample_dimensions) = lx_size;
    std::get<1>(sample_dimensions) = ly_size;
    std::get<2>(sample_dimensions) = lz_size;

/// Output to the screen/console
//    cout<< "______________________________________________________________________________________" << endl;
    cout << "The Multiphysics module specifications:\t\t" << endl;
    cout << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(sample_dimensions) << ", y: " << std::get<1>(sample_dimensions) << ", z: " << std::get<2>(sample_dimensions) << endl;
    cout << endl;

/// Output into .log file
//    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
    Out_logfile_stream << "The Writer module specifications:\t\t" << endl;
    Out_logfile_stream << "Sample dimensions are \t\t\t"s << " x: " << std::get<0>(sample_dimensions) << ", y: " << std::get<1>(sample_dimensions) << ", z: " << std::get<2>(sample_dimensions) << endl;
    Out_logfile_stream << endl;

    return;
} /// end of the bool MultiphysicsON() function
