///======================================= PCC Processing module ================================================================================ ///
///============================================================================================================================================= ///
///* The interface use functions from Processing_<***>_functions.h C++ libraries to generate quasi-random or non-random processes of           *///
///* labelling of the k-cells ( k={0,1,2,3} ) of the pre-constructed polytopal cell complex (PCC) 'M' using the whole set of incidence        *///
///* and adjacency matrices                                                                                                                   *///
///* ---------------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2024 years as a module of the PCC Processing Design code (CPD code) *///
///* A part of the MATERiA codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years   *///
/// https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1                                                                  *///
///===================================================================================================================================== ///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include <Eigen/SparseCore>

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../ini/ini_readers.h"
#include "../PCC_Objects.h"
#include "../PCC_Measures.h"

// Local
///---------------------------------------------------------
#include "functions/processing_assigned_labelling.h"
#include "functions/processing_induced_labelling.h"
#include "functions/processing_indexing.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs;
extern ofstream Out_logfile_stream;
extern string source_path, output_dir;
extern std::vector<std::string> PCCpaths;
extern int dim;

#include "PCC_Processing.h"
///* ========================================================= PCC PROCESSING FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Employ functions from Processing_<***>_functions.h C++ libraries to generate quasi-random or non-random processes of labelling of the k-cells ( k={0,1,2,3} )
 * of the pre-constructed polytopal cell complex (PCC) 'M' using the whole set of incidence and adjacency matrices.
 * All the initial settings are written in 'processing.ini' file, including PCCpaths to the corresponding directories and execution types.
 * The key theoretical concepts are special (SCS), ordinary (OCS) k-cell sequences, state vectors (SVs), design vectors (DVs) and configurations.
 * @param configuration
 * @return CellDesign object
 */
CellDesign PCC_Processing(Config &configuration) {
/// Main output of the module 'special_cells_design' (CD) - class. In particular, it contains (1) special_nodes_sequence, (2) special_edges_sequence, (3) special_faces_sequence (in the 2D and 3D cases), and (4) special_polyhedrons_sequence (in the 3D case)
    CellDesign CD;

    std::vector<std::vector<unsigned int>> Configuration_sState = configuration.Get_Configuration_sState(); // definition of the local State Vectors of special cells with values from the object of Config class
    std::vector<std::vector<unsigned int>> Configuration_cState = configuration.Get_Configuration_iState(); // definition of the local State Vectors of induced cells with values from the object of Config class

/// Processing vector variables
    std::vector<unsigned int> cell_strip_distribution; // vector of positive integers containing "discrete" length distribution of special chains/strips of k-cells
    std::vector <unsigned int> special_x_sequence; // variable sequence of 'special' cell structure (SCS) (different from 'ordinary' cell  structure (OCS)) k-cells in order of their generation.
    std::vector <unsigned int> induced_x_sequence; // variable sequence of 'induced' cell  structure (ICS) k-cells in order of their generation as the result of a kinetic process (historically they are also called "fractured" cells).
    std::vector<std::vector<unsigned int>> special_x_series; // vector of series of special k-cells (strips/chains)
    std::vector<std::vector<unsigned int>> induced_x_series; // vector of series of induced k-cells (like crack paths)
    std::vector<Agglomeration> agglomeration_x_sequence; // vector of agglomerations

    double mu_L = 1.0, sigma_L = 0.0; // mean and dispersion for a lengthy defect sequences distribution - used only in the case of the log-normal distribution for 'L' PCC_Processing execution mode ('pp_mode' in the config/processing.ini file).
    unsigned int bins_number_L = 10;

    Out_logfile_stream.open(output_dir + "Processing_Design.log"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    cout << "=========================================================================" << endl; Out_logfile_stream << "==============================================================================================================================================================" << endl;

/// Read configuration from processing.ini file :: the number of special cell types and calculation parameters.
    std::vector<vector<double>> max_sfractions_vectors(4), max_ifractions_vectors(4); // double vectors containing maximal values of fractions for Processing module execution :: [0][..] - nodes, [1][..] - edges, [2][..] - faces, [3][..] - polyhedrons
    std::vector<string> stype_vector(4), itype_vector(4); // special_types and induced_types vector of strings corresponding to the Processing execution type ON/OFF in 'config/main.ini file' - {0,1} for all the possible 'k' values of k-cell read from the file 'config/processing.ini'; In the both vectors :: [0] - nodes, [1] - edges, [2] - faces, [3] - polyhedrons
    std::vector<string> sequence_source_paths(4); // k-sequence paths for the reading them from file(s) - used only in the 'S' PCC_Processing execution mode ('pp_mode' in the config/processing.ini file)
    bool multiplexity = 0; // Bool variable indicating one (multiplexity = 0) or many (multiplexity = 1) labels can be assigned for each k-cell in a PCC
    std::vector<double> pindex_vector(4); // supplementary index used for multiplexity bool parameter - read from the 'config/processing.ini' file

// Reading of the configuration from the 'config/processing.ini' file
    config_reader_processing(source_path, sequence_source_paths, max_sfractions_vectors, max_ifractions_vectors, mu_L, sigma_L, bins_number_L, stype_vector, itype_vector, pindex_vector, Out_logfile_stream); // void function

///* Cases for Processing types // cell type k :: 0 - nodes, 1 - edges, 2 - faces, 3 -polyhedrons - must coincide with the indexing of the CellNumbs.at(cell_type) vector (!) *///
///* ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- *///
/// Here '(dim - 3)' term is important for 1D and 2D cases (!) as in these cases there are no polyhedrons or even faces (2D) in the PCC                                      ///
    for (int cell_type = (3 + (dim - 3)); cell_type >= 0; --cell_type) { // Loop over all types of k-cells in the PCC

        special_x_sequence.clear(); special_x_series.clear(); // clearence of vectors for each new cell type
        agglomeration_x_sequence.clear(); // clearence of agglomeration file

    /// I. Beginning of the processing of 'special' k-cells
    ///=======================================================
        if (stype_vector.at(cell_type) == "R" && max_sfractions_vectors.at(cell_type).size() > 0) { //  Random separate cells generation processing
            cout << "Random (R) mode processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "Random (R) mode processing in operation: cell_type : "s << cell_type << endl;

            multiplexity = (bool) pindex_vector.at(cell_type); // convert 'pindex' read from the 'config/processing.ini' file for the specific 'cell_type' to a bool variable 'multiplexity'.
            special_x_series = Processing_Random(cell_type, Configuration_sState, max_sfractions_vectors, multiplexity); // defined in the assigned labelling library; "\functions" subfolder

        /// Writing 'special_x_sequence' - each k-cell number appeared only once in the sequence. Configuration_sState and State Vectors show possible 'agglomeration' - several similar label per k-cell
            special_x_sequence.clear();
            for (auto it = Configuration_sState[cell_type].begin(); it != Configuration_sState[cell_type].end(); ++it)
                if(*it > 0) {
                    special_x_sequence.push_back(distance(Configuration_sState[cell_type].begin(), it)); // add new element to the special_x_cells_sequence
                } // end if(it)

            if(multiplexity != 0) {  /// Agglomerations counter
                std::vector<unsigned int> probe_state_vector(CellNumbs.at(cell_type), 0); // probe vector of the total x_series size filled with 0s
                std::fill(probe_state_vector.begin(), probe_state_vector.end(),0);
                for (auto kcell_seq: special_x_series) {// each strip/chain in the series of strips
                    for (unsigned int kcell: kcell_seq) { // in each strip/chain
                        probe_state_vector.at(kcell) += 1; // change element of the State Vector
                    }
                }

                for (auto it = probe_state_vector.begin(); it != probe_state_vector.end(); ++it) {
                    if (*it > 1) {
                        agglomeration_x_sequence.push_back( {(unsigned int) distance(probe_state_vector.begin(), it),*it} ); // Agglomeration(unsigned int AFace, unsigned int AglPower);
//REPAIR     cout << " number2: " << agglomeration_x_sequence.back().Get_agglomeration_kcell_number()<< " power2: " << agglomeration_x_sequence.back().Get_agglomeration_power() << endl;
                    } // end if()
                } // end for()
            } // end of if(multiplexity != 0) - agglomerations counter

        } // End of 'R' type simulations - // end  if (stype_vector.at(cell_type) == "R" && ..)

        else if (stype_vector.at(cell_type) == "L" && max_sfractions_vectors.at(cell_type).size() > 0) { //  Random lengthy strips (chains) of cells generation processing
            cout << "Random strips/chains (L) mode processing in operation: cell_type : "s << cell_type << endl << endl; Out_logfile_stream << "Random strips/chains (L) mode processing in operation: cell_type : "s << cell_type << endl << endl;
            cout << "Average (mu) and dispersion (sigma): " << endl << mu_L << "  " << sigma_L << endl << "Number of bins: " << bins_number_L << endl; Out_logfile_stream << "Average (mu) and dispersion (sigma): " << endl << mu_L << "  " << sigma_L << endl << "Number of bins: " << bins_number_L << endl;
            /// Obtaining the distribution of strip/chain lengths
            std::vector<double> strip_lenghts_distribution = Log_normal_distribution(mu_L, sigma_L, bins_number_L); // double valued "continuous" distribution, where Log_normal_distribution() function is for obtaining strip_lenghts_distribution

// REPAIR             ofstream Distributions_stream; Distributions_stream.open(output_dir + "Strip_distributions.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
            cell_strip_distribution.clear(); // clearing strip/chain length distribution vector for new 'cell_type' iterations
            for (auto itr = strip_lenghts_distribution.begin(); itr != strip_lenghts_distribution.end(); ++itr) {
// REPAIR                cout << " strip_lenghts_distribution " << max_sfractions_vectors[cell_type][0] * CellNumbs.at(cell_type) * (*itr) << endl;
///                    cell_strip_distribution.push_back(vop.at(0) * CellNumbs.at(2) * (*itr) / (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
                cell_strip_distribution.push_back( (max_sfractions_vectors.at(cell_type).at(0) * CellNumbs.at(2) * (*itr) * 5.0/ (double) bins_number_L) / (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
//                    cout << cell_strip_distribution.back() << "  ";
            }

            /// Output of the strip/chain lengths distribution
            double itr_sum = 0, num = 1;
// REPAIR               for (auto idtr = strip_lenghts_distribution.begin(); idtr != strip_lenghts_distribution.end(); ++idtr) {
            for (auto idtr = cell_strip_distribution.begin(); idtr != cell_strip_distribution.end(); ++idtr) {
                cout << *idtr << "  "; Out_logfile_stream << *idtr << "  ";

// REPAIR       Distributions_stream << *idtr << " \t";
                itr_sum += *idtr*num; //*(std::distance(cell_strip_distribution.begin(), idtr) + 1.0);
                ++ num;
            }
// REPAIR             Distributions_stream << endl;
            cout << endl << "Inclusion strip number:\t\t" << itr_sum << "\t\tSpecial faces Number:\t\t" << CellNumbs.at(2)*max_sfractions_vectors.at(cell_type).at(0) << endl << endl; //<< "\t\tDifference:\t\t" << abs(itr_sum - CellNumbs.at(2)*vop.at(0))*100.0/ double(CellNumbs.at(2)*vop.at(0))
            Out_logfile_stream << endl << "Inclusion strip number:\t\t" << itr_sum << "\t\tSpecial faces Number:\t\t" << CellNumbs.at(2)*max_sfractions_vectors.at(cell_type).at(0) << endl << endl;
            itr_sum = 0;
// REPAIR     }    Distributions_stream.close(); //exit(0);

            /// Random_Strips_Distribution() function call
            special_x_series = Processing_Random_Strips(cell_type, cell_strip_distribution, Configuration_sState, max_sfractions_vectors); // series of k-cells for each strip/chain

            /// Writing 'special_x_sequence' - each k-cell number appeared only once in the sequence. Configuration_sState and State Vectors show possible 'agglomeration' - several similar label per k-cell
            for (auto it = Configuration_sState[cell_type].begin(); it != Configuration_sState[cell_type].end(); ++it)
                if(*it > 0) {
                    special_x_sequence.push_back(distance(Configuration_sState[cell_type].begin(), it)); // add new element to the s_cells_sequence
                } // end if()

            /// Agglomeration counter
            std::vector<unsigned int> probe_state_vector(CellNumbs.at(cell_type)); // probe vector of the total x_series size filled with 0s
            std::fill(probe_state_vector.begin(), probe_state_vector.end(),0);
            for (auto kcell_seq : special_x_series) {// each strip/chain in the series of strips
                for (unsigned int kcell : kcell_seq) { // in each strip/chain
                    probe_state_vector.at(kcell) += 1; // change element of the State Vector
                }
            }

            for (unsigned int kcell : probe_state_vector) {
                for (auto it = probe_state_vector.begin(); it != probe_state_vector.end(); ++it) {
                    if (*it > 1) {
                        agglomeration_x_sequence.push_back( {(unsigned int) distance(probe_state_vector.begin(),it), *it} ); // Agglomeration(unsigned int AFace, unsigned int AglPower);
//REPAIR cout << " number2: " << agglomeration_x_sequence.back().Get_agglomeration_kcell_number()<< " power2: " << agglomeration_x_sequence.back().Get_agglomeration_power() << endl;
                    } // end if()
                } // end for()
            } // end for (unsigned int kcell : probe_state_vector)

        } // End of 'L' type simulations - else if (stype_vector.at(cell_type) == "L" && .. )
  /*
        else if (stype_vector.at(cell_type) == "F" && max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
            // processing index :: 0 - direct special faces assignment;  1 - crystallographic ; 2 - configurational TJs-based entropy (deviatoric); //        if (pindex_vector.at(cell_type) == 0) { //        } else if (pindex_vector.at(cell_type) == 1) {
            cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
          // if(cell_type == 2 + (dim - 3))             // cell type = 2 -> faces
//            double Configuration_Entropy(std::vector<int> const &TJsTypes);
            double (*conf_entropy_jfractions) (std::vector<double> const&j_fractions);
            conf_entropy_jfractions = Configuration_Entropy;

//            double Configuration_Entropy(std::vector<int> const &TJsTypes);
            double (*conf_entropy_jtypes) (std::vector<double> const&TJsTypes);
            conf_entropy_jtypes = Configuration_Entropy;

            //std::vector<unsigned int> Processing_maxFunctional(int cell_type, std::vector<std::vector<int>> &Configuration_State, std::vector<std::vector<double>> const &max_fractions_vectors, bool multiplexity);
            special_x_sequence = Processing_maxFunctional(cell_type, Configuration_sState, max_sfractions_vectors, multiplexity, conf_entropy_jtypes);
        } // End of 'F' type simulations (elseif)

        else if (stype_vector.at(cell_type) == "D" && max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
            cout << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "Min (MAX-deviator) Functional processing in operation: cell_type : "s << cell_type << endl;
///            if (max_fractions_vectors.at(cell_type).size() > 0)
///            special_x_sequence = Processing_minConfEntropy(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));

        } // End of 'D' [S min] type simulations (elseif)
*/
        else if (stype_vector.at(cell_type) == "S") { /// Reading structure from file
            vector<unsigned int> special_x_design;
            char* kseq_sourcepath = const_cast<char*>(sequence_source_paths.at(cell_type).c_str());
            special_x_design = VectorIReader(kseq_sourcepath);

            special_x_sequence.clear();
            for (auto it = special_x_design.begin(); it != special_x_design.end(); ++it) {
                if(*it > 0) {
                    special_x_sequence.push_back(distance(special_x_design.begin(), it));
                    // cout << special_x_sequence.back() << " ";
                }
            }
        // (!) Output +1 like in Neper, so he numbers should be modified back as -1
//            for (auto it = special_x_sequence.begin(); it != special_x_sequence.end(); ++it)
  //              special_x_sequence.at(distance(special_x_sequence.begin(), it)) = *it - 1;

            cout << endl;
            cout << "S processing mode! Special_x_sequence size: " << special_x_sequence.size() << endl << endl;
            cout << " Fraction " << cell_type << "-cells: " << (double) special_x_sequence.size()/ CellNumbs.at(cell_type) << endl;

/**
            // Cut up to max_fraction (!!)
            std::vector<unsigned int> temp_x_sequence = special_x_sequence; // temporarily new vector
            double total_max_sCell_fraction_processing = 0;
            for (int j = 0; j < max_sfractions_vectors[cell_type].size(); ++j)
                if(max_sfractions_vectors[cell_type][j] > 0)
                    total_max_sCell_fraction_processing += max_sfractions_vectors[cell_type][j];

            temp_x_sequence.clear();
            for (unsigned int k = 0; k < total_max_sCell_fraction_processing*CellNumbs.at(cell_type); ++k)
                temp_x_sequence.push_back(special_x_sequence.at(k));

            special_x_sequence.clear(); // putting everything back
            for(unsigned int p : temp_x_sequence)
                special_x_sequence.push_back(p);
            temp_x_sequence.clear();
**/
        // Update of the corresponding Configuration State vector
            Configuration_sState[cell_type].clear();
            std::vector<int> State_vector(CellNumbs.at(cell_type), 0);
            for (unsigned int k_cell : special_x_sequence) // fill state vector from cells_sequence
                State_vector.at(k_cell) = 1;

            for (int var : State_vector)
                Configuration_sState[cell_type].push_back(var);
        } // End of 'S' [reading from file] type simulations (elseif)

        else if (stype_vector.at(cell_type) == "Pi") { /// Reading structure from file
            cout << " Start TopDown_cell_indexing() based on the cell type k+1: " << cell_type << endl;

            /// II. Beginning of the processing of 'indexing' k-cells (TopDown and BottomUp)
            ///=======================================================
            std::vector<unsigned int> StateVector_indexing;
            StateVector_indexing = TopDown_cell_indexing(cell_type, Configuration_sState);

            /// 'special_x_sequence' - each k-cell number appeared only once in the sequence obtained by Configuration_sState
            for (auto it = StateVector_indexing.begin(); it != StateVector_indexing.end(); ++it)
                if(*it > 0) {
                    special_x_sequence.push_back(distance(StateVector_indexing.begin(), it)); // add new element to the s_cells_sequence
                } // end if(it)
            cout << " TopDown_cell_indexing() special_x_sequence size: " << special_x_sequence.size() << endl << endl;
        } // End of 'Pind' [TopDown indexing] type simulations (elseif)

        else if(max_sfractions_vectors[cell_type].size() > 0) cout << "ERROR [Processing] : unknown simulation type - please replace with 'R', 'L', 'F', 'D' or 'S'..!" << endl;

        ///* ONLY for 2-cells or 'grain boundaries' *///
        if (cell_type == (2 + (dim - 3))) { // '..+ (dim - 3)' because in the 2D case, boundaries become edges (!) and grains become faces of the corresponding 2D tessellation

            if (stype_vector.at(cell_type) == "Cm" &&
                max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
///                if (max_sfractions_vectors.at(cell_type).size() > 0)
///                    special_x_sequence = Processing_maxF_crystallographic(2, Configuration_sState, max_sfractions_vectors, pindex_vector.at(2));
            } // end of 'Cm' type simulations (elseif)

            else if (stype_vector.at(cell_type) == "Cr" &&
                     max_sfractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "MaxFunctional processing in operation: cell_type : "s << cell_type << endl;
///                if (max_sfractions_vectors.at(cell_type).size() > 0)
///                    special_x_sequence = Processing_maxP_crystallographic(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));
                    //special_x_sequence = Processing_Random_crystallographic(2, Configuration_sState, max_fractions_vectors, pindex_vector.at(2));
            } // end of 'Cr' type simulations (elseif)

            /// III. Beginning of the processing of 'induced' (of 'fractured') k-cells
            ///=======================================================

            if (itype_vector.at(cell_type) == "Km" && max_ifractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "Induced processing in operation: cell_type : "s << cell_type << endl;
                Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type << endl;

                if (max_ifractions_vectors.at(cell_type).size() > 0)
                     induced_x_sequence = PCC_Kinematic_cracking(cell_type, special_x_sequence, Configuration_cState, max_ifractions_vectors);
            } // End of 'Km' type simulations (elseif)

            else if (itype_vector.at(cell_type) == "Kn" && max_ifractions_vectors[cell_type].size() > 0) { // Maximum <functional> production
                cout << "Induced processing in operation: cell_type : "s << cell_type << endl; Out_logfile_stream << "Induced processing in operation: cell_type : "s << cell_type << endl;

///                if (max_ifractions_vectors.at(cell_type).size() > 0)
///                    induced_x_sequence = PCC_Kinetic_cracking(Configuration_sState, face_elastic_energies, large_crack);
            } // End of 'Km' type simulations (elseif)

            else if(max_ifractions_vectors[cell_type].size() > 0) cout << "ERROR [Processing] : unknown induced simulation type - please replace with 'Km' or 'Kn'..!" << endl;

        } // End of if (cell_type == 2)
// REPAIR    cout << "ctype_vector " << ctype_vector.at(cell_type + (3 - 3)) << "  " << max_cfractions_vectors[cell_type + (3 - 3)].size() << endl;

    /// Assigned sequences and series:
    cout << " special_x_sequence size " << special_x_sequence.size() << " cell type " << cell_type << endl << endl;
        CD.Set_special_sequence(special_x_sequence, cell_type); // (sequence, id)
        CD.Set_special_series(special_x_series, cell_type); // (series, id)
    // agglomerations
        CD.Set_agglomeration_sequence(agglomeration_x_sequence, cell_type); // (sequence, k-cell ID)
        CD.Set_special_configuration(Configuration_sState.at(cell_type), cell_type); // (configuration/design, id) - design vector with types

    /// Induced sequences and series:
        CD.Set_induced_sequence(induced_x_sequence, cell_type); // (sequence, ctype)

    } // END of for (int cell_type = 3; cell_type >= 0; --cell_type)

    if (configuration.Get_main_type() == "LIST"s) {
        cout << endl;
        Out_logfile_stream << endl;
        cout << "n-sequence size: " << CD.Get_n_special_sequence().size() << endl; Out_logfile_stream << "n-sequence size: " << CD.Get_n_special_sequence().size() << endl;
        cout << "e-sequence size: " << CD.Get_e_special_sequence().size() << endl; Out_logfile_stream << "e-sequence size: " << CD.Get_e_special_sequence().size() << endl;
        cout << "f-sequence size: " << CD.Get_f_special_sequence().size() << endl; Out_logfile_stream << "f-sequence size: " << CD.Get_f_special_sequence().size() << endl;
        cout << "p-sequence size: " << CD.Get_p_special_sequence().size() << endl; Out_logfile_stream << "p-sequence size: " << CD.Get_p_special_sequence().size() << endl;
        cout << endl; Out_logfile_stream << endl;

        cout << "n-induced-sequence size: " << CD.Get_n_induced_sequence().size() << endl; Out_logfile_stream << "n-induced-sequence size: " << CD.Get_n_induced_sequence().size() << endl;
        cout << "e-induced-sequence size: " << CD.Get_e_induced_sequence().size() << endl; Out_logfile_stream << "e-induced-sequence size: " << CD.Get_e_induced_sequence().size() << endl;
        cout << "f-induced-sequence size: " << CD.Get_f_induced_sequence().size() << endl; Out_logfile_stream << "f-induced-sequence size: " << CD.Get_f_induced_sequence().size() << endl;
        cout << "p-induced-sequence size: " << CD.Get_p_induced_sequence().size() << endl; Out_logfile_stream << "p-induced-sequence size: " << CD.Get_p_induced_sequence().size() << endl;
        cout << endl; Out_logfile_stream << endl;

        cout << "n-design vector size: " << CD.Get_n_design().size() << endl; Out_logfile_stream << "n-design vector size: " << CD.Get_n_design().size() << endl;
        cout << "e-design vector size: " << CD.Get_e_design().size() << endl; Out_logfile_stream << "e-design vector size: " << CD.Get_e_design().size() << endl;
        cout << "f-design vector size: " << CD.Get_f_design().size() << endl; Out_logfile_stream << "f-design vector size: " << CD.Get_f_design().size() << endl;
        cout << "p-design vector size: " << CD.Get_p_design().size() << endl; Out_logfile_stream << "p-design vector size: " << CD.Get_p_design().size() << endl;
        cout << endl; Out_logfile_stream << endl;
    } // End of if (configuration.Get_main_type() == "LIST"s)
    Out_logfile_stream.close();

    return CD;
} /// The END of PCC_Processing() function

/*
 *
        cell_strip_distribution.clear(); // clearing strip/chain length distribution vector for new 'cell_type' iterations
            for (auto  itr = strip_lenghts_distribution.begin(); itr != strip_lenghts_distribution.end(); ++itr) {
//REPAIR                cout << " strip_lenghts_distribution " << max_sfractions_vectors[cell_type][0] * CellNumbs.at(cell_type) * (*itr) << endl;
                cell_strip_distribution.push_back(max_sfractions_vectors[cell_type][0] * CellNumbs.at(cell_type) * (*itr) / (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
            }

            /// Output of the strip/chain lengths distribution
            for (auto  itr = cell_strip_distribution.begin(); itr != cell_strip_distribution.end(); ++itr) {
                cout << *itr << "  "; Out_logfile_stream << *itr << "  ";
            }
            cout << endl << endl; Out_logfile_stream << endl << endl;

 */