/// Attached user defined C++ libraries:
///-------------------------------------
///-------------------------------------
// #include "TJsLab.h"
 #include "LaplaciansLab.h"
#include "analytical_solutions.h"
///-------------------------------------

///Structure characterisation tool
ProcessedComplex PCC_StructureCharacterisation(CellsDesign &new_cells_design) {
    ProcessedComplex PCC_characteristics; // module output
/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen
    std::vector<int> charlabs_polyhedrons, charlabs_faces, charlabs_edges, charlabs_nodes;
    std::vector<double> ConfigVector = config_reader_characterisation(source_path, charlabs_polyhedrons, charlabs_faces, charlabs_edges, charlabs_nodes, Out_logfile_stream); // vector<double> (!) from ini_readers.h library

/// # 1 # Creation of the State Vectors
// Faces
    std::vector<unsigned int> face_sequences_vector; // number of different sequences based on the face design vector
    std::vector<int> face_states_vector;

    for (auto fseq : new_cells_design.Get_f_sequence()) {
        face_sequences_vector.push_back(fseq);
        /// forming a new vector for characterisation
        if ( face_sequences_vector.size() % (int) (2.0*std::log(new_cells_design.Get_f_sequence().size())) == 0) // only some values - each 2.0*std::log(*.size()) - are written in the face_sequences_vector
            PCC_characteristics.face_process_seq.push_back(face_sequences_vector);

    } // end for(auto fseq...)

    for (int i = 0; i < 4; ++i) { // for all types of cells

    /// # 2 # Measures: fractions and entropies
        if(i == 1 && charlabs_edges.at(0) == 1) { // Edges lab
            vector<int> EdgeTypes_char(CellNumbs.at(1 + (dim - 3)), 0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
            std::vector<double> j_edge_fractions(dim+1,0), d_edge_fractions(dim,0); // fractions of (1) edges of different types and (2) edges of different degrees

            /// for each state:
            tuple<double,double> conf_entropy_t;
            for (auto current_seq : PCC_characteristics.face_process_seq) {
                EdgeTypes_char.clear();
                std::fill(j_edge_fractions.begin(), j_edge_fractions.end(),0);
                std::fill(d_edge_fractions.begin(), d_edge_fractions.end(),0);

                EdgeTypes_char = Edge_types_byFaces(CellNumbs, current_seq, j_edge_fractions, d_edge_fractions);

                conf_entropy_t = Configuration_Entropy_tuple(j_edge_fractions); // conf entropy

                PCC_characteristics.e_entropy_mean_vector.push_back(std::get<0>(conf_entropy_t)); // std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions) based on Edges fraction vector
                PCC_characteristics.e_entropy_skrew_vector.push_back(std::get<1>(conf_entropy_t)); // tuple
                //// MUST BE IMPROVED (!) :
                PCC_characteristics.e_entropy_full_vector.push_back( std::get<0>(conf_entropy_t) + std::get<1>(conf_entropy_t));

                PCC_characteristics.je_fractions_vector.push_back(j_edge_fractions);
                PCC_characteristics.de_fractions_vector.push_back(d_edge_fractions);
            }

            if(charlabs_edges.at(4) == 1) { // analytical
            /// Measures: Analytical solutions for edge fractions and configuration entropies
                std::vector<tuple<double, double>> AnalyticalRandEntropies, AnalyticalCrystEntropies;
                PCC_characteristics.j_analytical_rand_vector = TJsAnalytics_random(CellNumbs.at(2 + (dim - 3)), AnalyticalRandEntropies);
                PCC_characteristics.AnRandEntropies_vector = AnalyticalRandEntropies;

                PCC_characteristics.j_analytical_cryst_vector = TJsAnalytics_crystallography(CellNumbs.at(2 + (dim - 3)), AnalyticalCrystEntropies);
                PCC_characteristics.AnCrystEntropies_vector = AnalyticalCrystEntropies;

                // analytical edge degree fractions
                PCC_characteristics.d_analytical_rand_vector = TJDsAnalytics_random(CellNumbs.at(2 + (dim - 3)));
                PCC_characteristics.d_analytical_cryst_vector = TJDsAnalytics_crystallography(CellNumbs.at(2 + (dim - 3)));
            } // enf if(charlabs_edges.at(4) == 1)

            } //         if(i == 1 && charlabs_edges.at(0) == 1) { // Edges lab
    } // end for (int i = 0; i < 4; ++i)

    return PCC_characteristics;
}