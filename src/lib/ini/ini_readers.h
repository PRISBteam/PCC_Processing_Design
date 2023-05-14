/// Author: Dr ELijah Borodin, 2023
/// Manchester
/// Library of specific functions related to the PCC Processing Design code for reading its .ini files

using namespace std; //Standard namespace

/// ================== # 1 # Initial configuration - reading and output ==================
std::vector<int> config_reader_main(std::string &source_path, std::string &source_dir, std::string &output_dir, std::string &main_type, std::ofstream &Out_logfile_stream) {
    std::vector<int> res;
    bool isSectionON = 0, isProcessingON = 0, isCharacterisationON = 0, isKineticON = 0, isMultiphysicsON = 0, isWriterON = 0;
    std::string isSection, isProcessing, isCharacterisation, isKinetic, isMultiphysics, isWriter;

    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "main.ini"s);
    mINI::INIStructure main_ini;
    file.read(main_ini);

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
    res.push_back(stoi(problem_dimension)); // res[0]

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
    if (main_ini.has("modules")) {
        auto& collection = main_ini["modules"];
        if (collection.has("PCC_Writer"))
            isWriter = main_ini.get("modules").get("PCC_Writer");
    }
    // ON/OFF IDs
    if (isSection == "ON") { isSectionON = 1; res.push_back(1); } else res.push_back(0); // res[1] - Section
    if (isProcessing == "ON") { isProcessingON = 1; res.push_back(1); } else res.push_back(0); // res[2] - Processing
    if (isCharacterisation == "ON") { isCharacterisationON = 1; res.push_back(1); } else res.push_back(0); // res[3] - Characterisation
    if (isMultiphysics == "ON") { isMultiphysicsON = 1; res.push_back(1); } else res.push_back(0); // res[4] - Multiphysics
    if (isKinetic == "ON") { isKineticON = 1; res.push_back(1); } else res.push_back(0); // res[5] - Kinetic
    if (isWriter == "ON") { isWriterON = 1; res.push_back(1); } else res.push_back(0); // res[6] - Writer

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
    cout << "Simulation type:\t"s << main_type << endl;
    cout << endl;
    cout << "Source directory:\t"s << source_dir << endl;
    cout << "Output directory:\t"s << output_dir << endl;
    cout << endl;
    if (isSectionON == 1) cout << "ON    | PCC_Section"s << endl;
    else cout << "OFF    | PCC_Section"s << endl;
    if (isProcessingON == 1) cout << "ON    | PCC_Processing"s << endl;
    else cout << "OFF    | PCC_Processing"s << endl;
    if (isCharacterisationON == 1) cout << "ON    | PCC_Characterisation"s << endl;
    else cout << "OFF    | PCC_Characterisation"s << endl;
    if (isMultiphysicsON == 1) cout << "ON    | PCC_Multiphysics"s << endl;
    else cout << "OFF    | PCC_Multiphysics"s << endl;
    if (isKineticON == 1) cout << "ON    | PCC_Kinetic"s << endl;
    else cout << "OFF    | PCC_Kinetic"s << endl;
    if (isWriterON == 1) cout << "ON    | PCC_Writer"s << endl;
    else cout << "OFF    | PCC_Writer"s << endl;
    cout << endl;

/// Output into .log file
    Out_logfile_stream << "The problem dimension that is the maximum value k_max of k-cells in the PCC:\t\t|\t\t"s << "dim = " << res.at(0) << endl;
    Out_logfile_stream << "Simulation type:\t"s << main_type << endl;
    Out_logfile_stream << endl;
    Out_logfile_stream << "Source directory:\t"s << source_dir << endl;
    Out_logfile_stream << "Output directory:\t"s << output_dir << endl;
    Out_logfile_stream << endl;
    if (isSectionON == 1) Out_logfile_stream << "ON    | PCC_Section"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Section"s << endl;
    if (isProcessingON == 1) Out_logfile_stream << "ON    | PCC_Processing"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Processing"s << endl;
    if (isCharacterisationON == 1) Out_logfile_stream << "ON    | PCC_Characterisation"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Characterisation"s << endl;
    if (isMultiphysicsON == 1) Out_logfile_stream << "ON    | PCC_Multiphysics"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Multiphysics"s << endl;
    if (isKineticON == 1) Out_logfile_stream << "ON    | PCC_Kinetic"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Kinetic"s << endl;
    if (isWriterON == 1) Out_logfile_stream << "ON    | PCC_Writer"s << endl;
    else Out_logfile_stream << "OFF    | PCC_Writer"s << endl;
    Out_logfile_stream << endl;

    return res;
} /// END of config_reader_main function



/// ================== # 2 # Initial PROCESSING module configuration - reading and output ==================
void config_reader_processing(std::string &source_path, std::vector<string> &sequence_source_paths, std::vector<vector<double>> &max_fractions_vectors, double &mu, double &sigma, std::vector<string> &ptype_vector, std::vector<double> &pindex_vector, std::ofstream &Out_logfile_stream) {
    // ini files reader - external (MIT license) library
    mINI::INIFile file(source_path + "processing.ini"s);
    mINI::INIStructure processing_ini;
    file.read(processing_ini);

// I: cell types and max fractions and processing modes
/// Polyhedrons
//processing_mode
    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("pp_mode"))
        {
            ptype_vector.at(3) = processing_ini.get("polyhedrons").get("pp_mode");
        } }

    string pseq_source;
    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("source"))
            pseq_source = processing_ini.get("polyhedrons").get("source");
    }

    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("pp_index"))
            pindex_vector.at(3) = stod(processing_ini.get("polyhedrons").get("pp_index"));
        // R(0) - R, S(1) - Smax, S(0) - Smin, I(x.x) - index mode
    }

    string ptypes_number_string;
    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("polyhedron_types_number"))
            ptypes_number_string = processing_ini.get("polyhedrons").get("polyhedron_types_number"); // [2]
    }

// fractions
    string p1_max, p2_max, p3_max;
    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction1"))
            p1_max = processing_ini.get("polyhedrons").get("pmax_fraction1");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p1_max)); // 3 - polyhedra

    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction2"))
            p2_max = processing_ini.get("polyhedrons").get("pmax_fraction2");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p2_max)); // 3 - polyhedra

    if (processing_ini.has("polyhedrons")) {
        auto& collection = processing_ini["polyhedrons"];
        if (collection.has("pmax_fraction3"))
            p3_max = processing_ini.get("polyhedrons").get("pmax_fraction3");
    }
    if (stoi(ptypes_number_string) > 0) max_fractions_vectors.at(3).push_back(stod(p3_max)); // 3 - polyhedra

/// Faces
//processing_mode
    if (processing_ini.has("faces")) {
        auto& collection = processing_ini["faces"];
        if (collection.has("pf_mode"))
        {
            ptype_vector.at(2) = processing_ini.get("faces").get("pf_mode");
        } }

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
            eseq_source = processing_ini.get("nodes").get("source");
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
    if (etypes_number_string != "0") {
        //edges
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

    /// Output to the screen/console
    cout<< "______________________________________________________________________________________" << endl;
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
    cout<< "______________________________________________________________________________________" << endl;

/// Output into .log file
    Out_logfile_stream<< "______________________________________________________________________________________" << endl;
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
    Out_logfile_stream<< "______________________________________________________________________________________" << endl;

    return config_characterisation_vector;
}

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
