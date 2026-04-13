///======================================= PCC Characterisation module ============================================================================///
///==========================================================================================================================================///
///* The interface use functions from PCC_Characterisation/functions C++ libraries to analyse both the design vectors and sequences of     *///
///* State vectors provided by the PCC_Processing modules and written in an object of the ProcessedComplex class from PCC_Objects library *///
///* -----------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design code (CPD code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years *///
/// (https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1)                                                            *///
///==================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include "../../../src/lib/external/Eigen-5.0/SparseCore"

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Measures.h"
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
#include "functions/spectral_analysis.h"
// #include "functions/analytical_solutions.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::string source_path;
extern ofstream Out_logfile_stream;
extern std::vector<unsigned int> CellNumbs;
extern std::vector<std::string> paths_to_PCC_matrices;
extern int PCC_dimension;

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'

#include "PCC_Characterisation.h"
///* ========================================================= PCC CHARACTERISATION FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///

ProcessedComplex PCC_StructureCharacterisation(CellDesign &new_cells_design) {

/// Main output of the module
    ProcessedComplex PCC_characteristics; // module output

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen.
    std::vector<int> charlabs_polyhedrons, charlabs_faces, charlabs_edges, charlabs_nodes, charlabs_laplacians;
    std::vector<double> ConfigVector = config_reader_characterisation(charlabs_polyhedrons, charlabs_faces, charlabs_edges, charlabs_nodes, charlabs_laplacians); // vector<double> (!) from ini_readers.h library

/// # 1 # Creation of the State Vectors
// Faces
    std::vector<unsigned int> face_sequences_vector; // number of different sequences based on the face design vector
    std::vector<int> face_states_vector;

//    for (auto fseq : new_cells_design.Get_f_special_sequence()) {
        for (auto fseq : new_cells_design.Get_f_induced_sequence()) {

        face_sequences_vector.push_back(fseq);
        /// forming a new vector for characterisation
        ///if ( face_sequences_vector.size() % (int) (2.0*std::log(new_cells_design.Get_f_special_sequence().size())) == 0) // only some values - each 2.0*std::log(*.size()) - are written in the face_sequences_vector
        if ( face_sequences_vector.size() % (int) (2.0*std::log(new_cells_design.Get_f_induced_sequence().size())) == 0) // only some values - each 2.0*std::log(*.size()) - are written in the face_sequences_vector

            PCC_characteristics.face_process_seq.push_back(face_sequences_vector);
    } // end for(auto fseq...)
///    cout << "new_cells_design.Get_f_special_sequence() size" << new_cells_design.Get_f_special_sequence().size() << endl;
    cout << "new_cells_design.Get_f_induced_sequence() size" << new_cells_design.Get_f_induced_sequence().size() << endl;
    Out_logfile_stream << "new_cells_design.Get_f_induced_sequence() size" << new_cells_design.Get_f_induced_sequence().size() << endl;
/**
    for (int i = 0; i < 4; ++i) { /// for all types of cells

        /// # 2 # Measures: fractions and entropies
        if(i == 1 && charlabs_edges.size() > 0) { // Edges lab
            if (charlabs_edges.at(0) == 1) {
                vector<int> EdgeTypes_char(CellNumbs.at(1 + (dim - 3)),
                                           0); // vector<int> in the form [ 0 2 3 3 2 1 ...] with the TJs type ID as its values
                std::vector<double> j_edge_fractions(dim + 1, 0), d_edge_fractions(dim,0); // fractions of (1) edges of different types and (2) edges of different degrees

                /// for each state:
                tuple<double, double> conf_entropy_t;
                for (auto current_seq: PCC_characteristics.face_process_seq) {
                    EdgeTypes_char.clear();
                    std::fill(j_edge_fractions.begin(), j_edge_fractions.end(), 0);
                    std::fill(d_edge_fractions.begin(), d_edge_fractions.end(), 0);

                    EdgeTypes_char = Edge_types_byFaces(CellNumbs, current_seq, j_edge_fractions, d_edge_fractions);

                    conf_entropy_t = Configuration_Entropy_tuple(j_edge_fractions); // conf entropy

                    PCC_characteristics.e_entropy_mean_vector.push_back(std::get<0>(
                            conf_entropy_t)); // std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions) based on Edges fraction vector
                    PCC_characteristics.e_entropy_skrew_vector.push_back(std::get<1>(conf_entropy_t)); // tuple
                    //// MUST BE IMPROVED (!) :
                    PCC_characteristics.e_entropy_full_vector.push_back(
                            std::get<0>(conf_entropy_t) + std::get<1>(conf_entropy_t));

                    PCC_characteristics.je_fractions_vector.push_back(j_edge_fractions);
                    PCC_characteristics.de_fractions_vector.push_back(d_edge_fractions);
                } // for (auto current_seq : PCC_characteristics.face_process_seq)

                if (charlabs_edges.at(4) == 1) { // analytical
                    /// Measures: Analytical solutions for edge fractions and configuration entropies
                    std::vector<tuple<double, double>> AnalyticalRandEntropies, AnalyticalCrystEntropies;
                    PCC_characteristics.j_analytical_rand_vector = TJsAnalytics_random(CellNumbs.at(2 + (dim - 3)),
                                                                                       AnalyticalRandEntropies);
                    PCC_characteristics.AnRandEntropies_vector = AnalyticalRandEntropies;

                    PCC_characteristics.j_analytical_cryst_vector = TJsAnalytics_crystallography(
                            CellNumbs.at(2 + (dim - 3)), AnalyticalCrystEntropies);
                    PCC_characteristics.AnCrystEntropies_vector = AnalyticalCrystEntropies;

                    // analytical edge degree fractions
                    PCC_characteristics.d_analytical_rand_vector = TJDsAnalytics_random(CellNumbs.at(2 + (dim - 3)));
                    PCC_characteristics.d_analytical_cryst_vector = TJDsAnalytics_crystallography(
                            CellNumbs.at(2 + (dim - 3)));
                } // enf if(charlabs_edges.at(4) == 1)

            } //  end of if(charlabs_edges.at(0) == 1)
        } //  end of if(i == 1 && charlabs_edges.at(0) == 1) { // Edges lab

    } // end for (int i = 0; i < 4; ++i)
**/

/// # 3 # Laplacians of specail cells
    if(charlabs_laplacians.at(1) == 1) {
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values // odir - source path

/// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
//    SpMat ANS(CellNumbs.at(0), CellNumbs.at(0)), AES(CellNumbs.at(1 + (dim - 3)), CellNumbs.at(1 + (dim - 3))), AFS(CellNumbs.at(2 + (dim - 3)), CellNumbs.at(2 + (dim - 3))), AGS(CellNumbs.at(3 + (dim - 3)), CellNumbs.at(3 + (dim - 3)));
//    ANS = SMatrixReader(paths.at(0), (CellNumbs.at(0)), (CellNumbs.at(0))); //all Nodes
//    ANS = 0.5 * (ANS + SparseMatrix<double>(ANS.transpose())); // Full matrix instead of triagonal
//    AES = SMatrixReader(paths.at(1 + (dim - 3)), (CellNumbs.at(1 + (dim - 3))), (CellNumbs.at(1 + (dim - 3)))); //all Edges
//    AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose())); // Full matrix instead of triagonal
//    AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2 + (dim - 3))), (CellNumbs.at(2 + (dim - 3)))); //all Faces
//    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal
//    AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3 + (dim - 3))), (CellNumbs.at(3 + (dim - 3)))); //all Volumes
//    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
        SpMat ENS(CellNumbs.at(0), CellNumbs.at(1)), FES(CellNumbs.at(1 + (PCC_dimension - 3)), CellNumbs.at(2 + (PCC_dimension - 3))), GFS(
                CellNumbs.at(2 + (PCC_dimension - 3)), CellNumbs.at(3 + (PCC_dimension - 3)));

        ENS = SMatrixReader(paths_to_PCC_matrices.at(4 + (PCC_dimension - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
        FES = SMatrixReader(paths_to_PCC_matrices.at(5 + (PCC_dimension - 3)), (CellNumbs.at(1 + (PCC_dimension - 3))),
                            (CellNumbs.at(2 + (PCC_dimension - 3)))); //all Edges-Faces
        GFS = SMatrixReader(paths_to_PCC_matrices.at(6 + (PCC_dimension - 3)), (CellNumbs.at(2 + (PCC_dimension - 3))),
                            (CellNumbs.at(3 + (PCC_dimension - 3)))); //all Faces-Grains

        double number_of_steps = (double) charlabs_laplacians.at(0); // reading from Characterisation.ini file (!)
        std::vector<unsigned int> new_current_seq;
        std::vector<double> new_Betti_numbers;

        std::vector<unsigned int> Betti_calc_time;


//        cout << endl;
//        cout << "we are here!!" << endl;
//        cout << endl;

//#pragma omp parallel for // parallel execution by OpenMP
        unsigned int d_seq = std::floor((double) PCC_characteristics.face_process_seq.size() / number_of_steps);
        for (unsigned int i = 0; i <= number_of_steps; ++i) { // ? number_of_steps ? 10
            new_current_seq = PCC_characteristics.face_process_seq.at(i*d_seq);
            cout << "---------------------------------------------------------------------------" << endl;
            cout << "[CHAR] Current special faces fraction: " << new_current_seq.size()/ (double) CellNumbs.at(2 + (PCC_dimension - 3)) << endl;
            Out_logfile_stream << "---------------------------------------------------------------------------" << endl;
            Out_logfile_stream << "[CHAR] Current special faces fraction: " << new_current_seq.size()/ (double) CellNumbs.at(2 + (PCC_dimension - 3)) << endl;
            // Elapsing time for 3 Betti numbers
            clock_t t;
            Betti_calc_time.push_back(clock());

            new_Betti_numbers = OperatorsBetti(new_current_seq, ENS, FES, GFS);
            cout << "[CHAR] p & Betti numbers:  " << new_current_seq.size()/ (double) CellNumbs.at(2) << "  " << new_Betti_numbers.at(0) << "  " << new_Betti_numbers.at(1) << "  " << new_Betti_numbers.at(2) << endl;
            Out_logfile_stream << "[CHAR] p & Betti numbers:  " << new_current_seq.size()/ (double) CellNumbs.at(2) << "  " << new_Betti_numbers.at(0) << "  " << new_Betti_numbers.at(1) << "  " << new_Betti_numbers.at(2) << endl;
            //
            Betti_calc_time.back() = clock() - Betti_calc_time.back();
            cout << "[CHAR] Calculation time [s] of the Betti numbers is equal to   " << Betti_calc_time.back()/ pow(10.0,6.0)<< endl;
            Out_logfile_stream << "[CHAR] Calculation time [s] of the Betti numbers is equal to   " << Betti_calc_time.back()/ pow(10.0,6.0)<< endl;

            PCC_characteristics.Betti_vector.push_back({new_current_seq.size()/(double) CellNumbs.at(2), new_Betti_numbers.at(0), new_Betti_numbers.at(1), new_Betti_numbers.at(2)});
        }
    } // end of  for (int i = 0; i < number_of_steps, ++i)

   return PCC_characteristics;
 } /// END of the ProcessedComplex PCC_StructureCharacterisation(CellDesign &new_cells_design) function

