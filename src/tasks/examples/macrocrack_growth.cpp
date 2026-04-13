/// Simulates Macrocrack growth as a sequence of cracks
/// Dr Elijah Borodin, 2024

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <random>

using namespace std;

#include "../lib/PCC_Support_Functions.h" // It must be here - first in this list of libraries (!)

#include "../lib/PCC_Measures.h"

#include "../lib/PCC_Objects.h"

#include "../lib/PCC_Subcomplex/PCC_Subcomplex.h"

#include "../lib/PCC_Multiphysics/PCC_Multiphysics.h"

#include "../lib/PCC_Processing/PCC_Processing.h"
#include "../../src/lib/PCC_Processing/functions/processing_assigned_labelling.h"
#include "../../src/lib/PCC_Processing/functions/processing_induced_labelling.h"

#include "../lib/PCC_Characterisation/PCC_Characterisation.h"
#include "../lib/PCC_Characterisation/functions/spectral_analysis.h"

extern std::ofstream Out_macrocrack_growth_stream;
extern std::string output_dir;
extern std::vector<unsigned int> CellNumbs;
extern std::vector<std::string> paths_to_PCC_matrices;
/*!
 * @details
 */
ProcessedComplex Macrocrack_growth(Config &configuration) {
    ProcessedComplex pcc_processed;  // function output

    double Main_time = 0.0, S_time = 0.0, M_time = 0.0, P_time = 0.0, C_time = 0.0, W_time = 0.0;

    /// KEYS =================================================================
    bool key_subcomplex = 1, key_processing = 1, key_multiphysics = 1, key_induced = 0, key_fBetti = 1, key_iBetti = 0;
    /// =================================================================
    /// =================================================================

    std::vector<std::vector<unsigned int>> special_face_series;
    std::vector <unsigned int> special_face_sequence, induced_face_sequence;
    std::vector<Agglomeration> agglomeration_face_sequence;
    std::vector<vector<double>> max_sfractions_vectors(4);
    max_sfractions_vectors[2].push_back({0});

    /// Defects
    std::vector<Macrocrack> crack_growth_series;
    std::vector<CellEnergies> new_cells_energies; // a class described in PCC_Objects.h contained (1) all the k-cell elastic energies and (2) all the k-cell thermal energies in the PCC

//    std::vector<std::vector<unsigned int>> Configuration_sState(4); // = configuration.Get_Configuration_sState(); // definition of the local State Vectors of special cells with values from the object of Config class
//    Configuration_sState = {State_p_vector, State_f_vector, State_e_vector, State_n_vector },
    std::vector<std::vector<unsigned int>> Configuration_sState = configuration.Get_Configuration_sState(); // definition of the local State Vectors of special cells with values from the object of Config class
    std::vector<std::vector<unsigned int>> Configuration_cState = configuration.Get_Configuration_iState(); // definition of the local State Vectors of induced cells with values from the object of Config class

    /// ====================== II. PCC Processing module ======================
    CellDesign new_cells_design; // a class described in PCC_Objects.h contained (1) all special k-cell sequences and (2) all the design_<*>_vectors for all k-cells in the PCC
    std::vector<CellDesign> new_cells_design_vector;
    cout << "-------------------------------------------------------------------------" << endl;
    Out_macrocrack_growth_stream << "-------------------------------------------------------------------------" << endl;
    cout << "START of the PCC Processing module " << endl;
    Out_macrocrack_growth_stream << "START of the PCC Processing module " << endl;

    SpMat ENS(CellNumbs.at(0), CellNumbs.at(1));
    ENS = SMatrixReader(paths_to_PCC_matrices.at(4), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges

    SpMat FES(CellNumbs.at(1), CellNumbs.at(2));
    FES = SMatrixReader(paths_to_PCC_matrices.at(5), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces

    SpMat GFS(CellNumbs.at(2), CellNumbs.at(3));
    GFS = SMatrixReader(paths_to_PCC_matrices.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains


    /// INCUSION POWDERS
    std::vector<Config> powder_confs(12);
// %
///            powder_confs.at(powder).Set_config();
//           double p1 = 0.01, p2 = 0.1, p3 = 0.3, mu_min = 0.00001*CellNumbs.at(2), mu_max = 0.0005*CellNumbs.at(2), sig_min = 10.0, sig_max = 100.0, bins = 20;
///    double p1 = 0.99, p2 = 0.99, p3 = 0.99, mu_min = 0.0, mu_max = 1.0, sig_min = 0.2, sig_max = 0.7, bins = 30;
    double p1 = 0.02, p2 = 0.1, p3 = 0.3, mu_min = 0.0, mu_max = 1.0, sig_min = 0.2, sig_max = 0.7, bins = 30;

    std::vector<std::vector<double>> vector_of_powders = {
            {p1, mu_min, sig_min, bins}, // 1
            {p2, mu_min, sig_min, bins}, // 2
            {p3, mu_min, sig_min, bins}, // 3
            {p1, mu_max, sig_min, bins}, // 4
            {p2, mu_max, sig_min, bins}, // 5
            {p3, mu_max, sig_min, bins}, // 6
            {p1, mu_max, sig_max, bins}, // 7
            {p2, mu_max, sig_max, bins}, // 8
            {p3, mu_max, sig_max, bins}, // 9
            {p1, mu_min, sig_max, bins}, // 10
            {p2, mu_min, sig_max, bins}, // 11
            {p3, mu_min, sig_max, bins}  // 12
    };
/*    std::vector<std::vector<double>> vector_of_powders = {
            {0.05, mu_max, sig_min, bins}, // 1
            {0.1, mu_max, sig_min, bins}, // 2
            {0.15, mu_max, sig_min, bins}, // 3
            {0.2, mu_max, sig_min, bins}, // 4
            {0.25, mu_max, sig_min, bins}, // 5
            {0.3, mu_max, sig_min, bins}, // 6
            {0.35, mu_max, sig_min, bins}, // 7
            {0.4, mu_max, sig_min, bins}, // 8
            {0.45, mu_max, sig_min, bins}, // 9
            {0.5, mu_max, sig_min, bins}, // 10
            {0.55, mu_max, sig_min, bins}, // 11
            {0.6, mu_max, sig_min, bins}  // 12
    };
*/
    std::vector<unsigned int> cell_strip_distribution; // vector of positive integers containing "discrete" length distribution of special chains/strips of k-cells
    std::vector<double> strip_lenghts_distribution;

    std::ofstream Face_Edge_indices, Cracked_pcc_out, Cracked_stat_out, Microcracks_out, Cracked_sfaces_numbers_in_plane_out, Yj_stat_out, Dspace_ifaces_out;
    Face_Edge_indices.open(output_dir + "Inclusion_Face_Edge_indices.txt"s, ios::trunc);
    Face_Edge_indices.close();
    Microcracks_out.open(output_dir + "Microcrack_faces_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Microcracks_out.close();
    Cracked_pcc_out.open(output_dir + "Macrocrack_faces_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Cracked_pcc_out.close();
    //Cracked_stat_out.open(output_dir + "Macrocrack_sfaces_numbers_in_planes_qeq10.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    //Cracked_stat_out.close();
    Cracked_sfaces_numbers_in_plane_out.open(output_dir + "Macrocrack_sfaces_numbers_in_planes.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Cracked_sfaces_numbers_in_plane_out.close();
    Yj_stat_out.open(output_dir + "Macrocrack_Y_junction_numbers_100k.txt"s, ios::trunc);
    Yj_stat_out.close();

    Face_Edge_indices.open(output_dir + "Inclusion_Face_Edge_indices.txt"s, ios::app);
    Face_Edge_indices << "\t\tBL\t\t" << "\t\tQL\t\t" << endl;
    Microcracks_out.open(output_dir + "Microcrack_faces_coordinates.txt"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    Cracked_pcc_out.open(output_dir + "Macrocrack_faces_coordinates.txt"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    //Cracked_stat_out.open(output_dir + "Macrocrack_sfaces_numbers_in_planes_qeq10.txt"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    Yj_stat_out.open(output_dir + "Macrocrack_Y_junction_numbers_100k.txt"s, ios::app);
    Cracked_sfaces_numbers_in_plane_out.open(output_dir + "Macrocrack_sfaces_numbers_in_planes.txt"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function
    Dspace_ifaces_out.open(output_dir + "Macrocrack_D_iFace_fractions.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    int powder_iterator = 0;

    /// ITERATION over the 'vector_of_powders'
    /// ==================================================================================
    /// ==================================================================================
    int vop_counter = 0;
    for (auto vop : vector_of_powders) {
        vop_counter++;
        Cracked_sfaces_numbers_in_plane_out << "\tPowder number\t" << vop_counter << endl;

        /// p according to the powder
        max_sfractions_vectors[2][0] = vop.at(0);

        if( key_processing == 1) {
            strip_lenghts_distribution.clear();

        strip_lenghts_distribution = Log_normal_distribution(vop.at(1), vop.at(2), vop.at(3)); // double valued "continuous" distribution, where Log_normal_distribution() function is for obtaining strip_lenghts_distribution
//REPAIR                 double itr_sum1 = 0.0; for (auto idtr = strip_lenghts_distribution.begin(); idtr != strip_lenghts_distribution.end(); ++idtr)  cout << *idtr << "  "; itr_sum1 += *idtr*5.0/ (double) bins; } cout << endl << "Sum:\t\t" << itr_sum1  << endl << endl;  exit(0);

        cell_strip_distribution.clear(); // clearing strip/chain length distribution vector for new 'cell_type' iterations
        for (auto itr = strip_lenghts_distribution.begin(); itr != strip_lenghts_distribution.end(); ++itr) {
// REPAIR                cout << " strip_lenghts_distribution " << max_sfractions_vectors[cell_type][0] * CellNumbs.at(cell_type) * (*itr) << endl;
///                    cell_strip_distribution.push_back(vop.at(0) * CellNumbs.at(2) * (*itr) / (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
            cell_strip_distribution.push_back((vop.at(0) * CellNumbs.at(2) * (*itr) * 5.0 / (double) bins) /
                                              (std::distance(strip_lenghts_distribution.begin(), itr) + 1.0)); // only for the first 'Xmax_fraction1' in the 'config/processing.ini' file values of max k-cell fractions
//                    cout << cell_strip_distribution.back() << "  ";
        }

        /// Output of the strip/chain lengths distribution
        double itr_sum = 0, num = 1;
// REPAIR               for (auto idtr = strip_lenghts_distribution.begin(); idtr != strip_lenghts_distribution.end(); ++idtr) {
cout << " Cell Strip Distribution:\t\t";
        for (auto idtr = cell_strip_distribution.begin(); idtr != cell_strip_distribution.end(); ++idtr) {
            cout << *idtr << "  ";
            Out_macrocrack_growth_stream << *idtr << " \t ";

// REPAIR       Distributions_stream << *idtr << " \t";
            itr_sum += *idtr * num; //*(std::distance(cell_strip_distribution.begin(), idtr) + 1.0);
            ++num;
        }
        Out_macrocrack_growth_stream << endl;
// REPAIR             Distributions_stream << endl;
        cout << endl << "Inclusion strip number:\t\t" << itr_sum << "\t\tSpecial faces Number:\t\t"
             << CellNumbs.at(2) * vop.at(0) << endl
             << endl; //<< "\t\tDifference:\t\t" << abs(itr_sum - CellNumbs.at(2)*vop.at(0))*100.0/ double(CellNumbs.at(2)*vop.at(0))
        Out_macrocrack_growth_stream << endl << "Inclusion strip number:\t\t" << itr_sum << "\t\tSpecial faces Number:\t\t"
                           << CellNumbs.at(2) * vop.at(0) << endl << endl;
        itr_sum = 0;
// REPAIR     }    Distributions_stream.close(); //exit(0);

/// A LOOP OVER AVAILABLE INCUSION POWDERS
///====================================================================================

        /// Random_Strips_Distribution() function call
        ///==================================================================
        special_face_series.clear();
///            cout << "special_in_pcc\t\t" << ; cout << endl;
//        for(unsigned int it = 0; it < Configuration_sState.size(); ++it)
//            Configuration_sState[2].at(it) = 0; //

            //// Zero initial state only here!!!
            std::fill(Configuration_sState[2].begin(), Configuration_sState[2].end(), 0);
           special_face_series = Processing_Random_Strips(2, cell_strip_distribution, Configuration_sState, max_sfractions_vectors); // series of k-cells for each strip/chain

        //        cout << "cell_type\t" << 2 << "\tConfiguration_State.at(cell_type)\t" << Configuration_sState.at(2).size() << endl; //"S_Vector " << S_Vector.size() << endl;
///            bool multiplexity = 1; // convert 'pindex' read from the 'config/processing.ini' file for the specific 'cell_type' to a bool variable 'multiplexity'.
//////// (alternative!code) special_face_series =  Processing_Random(2, Configuration_sState, max_sfractions_vectors, 1); // random flakes

        /// Writing 'special_x_sequence' - each k-cell number appeared only once in the sequence. Configuration_sState and State Vectors show possible 'agglomeration' - several similar label per k-cell
        special_face_sequence.clear();
            for (auto it = Configuration_sState[2].begin(); it != Configuration_sState[2].end(); ++it)
            if (*it > 0) {
                special_face_sequence.push_back(distance(Configuration_sState[2].begin(), it)); // add new element to the s_cells_sequence
            } // end if()

        // ================ Elapsing time for the Processing module ================
        unsigned int Processing_time = clock();
        P_time = (double) Processing_time - S_time - M_time - Main_time;
        cout << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl
             << endl; //cout << "-------------------------------------------------------------------------" << endl;
        Out_macrocrack_growth_stream << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl
                           << endl;

        cout << "-------------------------------------------------------------------------" << endl;
        Out_macrocrack_growth_stream << "-------------------------------------------------------------------------" << endl;
        cout << "START of the PCC Characterisation module" << endl;
        Out_macrocrack_growth_stream << "START of the PCC Characterisation module" << endl;
        cout << "=========================================================================" << endl;
        Out_macrocrack_growth_stream << "==============================================================================================================================================================" << endl;

        double FE_index = Face_edge_index(special_face_sequence, FES, 14.0);
        double NE_index = Node_edge_index(special_face_sequence, ENS, 3300.0);
            Face_Edge_indices << "Powder_number\t\t" << vop_counter << endl;
            Face_Edge_indices << FE_index << "\t" << NE_index << endl;

        /// Agglomeration counter
        ///===============================
        std::vector<unsigned int> probe_state_vector(CellNumbs.at(2)); // probe vector of the total x_series size filled with 0s
        std::fill(probe_state_vector.begin(), probe_state_vector.end(), 0);

            //                 cout << distance(Configuration_sState[2].begin(), it) << "\t";
            for (auto kcell_seq : special_face_series) {// each strip/chain in the series of strips
                for (unsigned int kcell : kcell_seq) { // in each strip/chain
                probe_state_vector.at(kcell) += 1; // change element of the State Vector
                 }
            }

            std::vector<Agglomeration> agglomerations_in_pcc;
///        for (unsigned int kcell: probe_state_vector) {
        for (auto it = probe_state_vector.begin(); it != probe_state_vector.end(); ++it) {
            if (*it > 1) {
                agglomerations_in_pcc.push_back({(unsigned int) distance(probe_state_vector.begin(), it),*it}); // Agglomeration(unsigned int AFace, unsigned int AglPower);
//REPAIR cout << " number2: " << agglomeration_x_sequence.back().Get_agglomeration_kcell_number()<< " power2: " << agglomeration_x_sequence.back().Get_agglomeration_power() << endl;
            } // end if()
        } // end for()
///        } // end for (unsigned int kcell : probe_state_vector)
        pcc_processed.agglomerations_in_powders.push_back(agglomerations_in_pcc); // Agglomeration(unsigned int AFace, unsigned int AglPower);

        /// de_fractions_vector
        pcc_processed.de_fractions_sface_vector.push_back(d_fractions_vector(EdgesTypesCalc(CellNumbs, special_face_sequence, FES)));

        /// Betti_Li_vector
            if (key_fBetti == 1) {
                std::vector<double> Betti_Li_vector = OperatorsBetti(special_face_sequence, ENS, FES, GFS);
                pcc_processed.Betti_0_sface.push_back(Betti_Li_vector.at(0));
                pcc_processed.Betti_1_sface.push_back(Betti_Li_vector.at(1));
                pcc_processed.Betti_2_sface.push_back(Betti_Li_vector.at(2));
                pcc_processed.inverse_connectivity_sface.push_back(std::log(Betti_Li_vector.at(0) / Betti_Li_vector.at(2)));
            }
        // ===== Elapsing time for the Characterisation module ================
        unsigned int Characterisation_time = clock();
        C_time = (double) Characterisation_time - S_time - M_time - P_time - Main_time;
        cout << "Characterisation time is equal to  " << C_time / pow(10.0, 6.0) << "  seconds" << endl
             << endl; //cout << "-------------------------------------------------------------------------" << endl;
        Out_macrocrack_growth_stream << "Characterisation time is equal to  " << C_time / pow(10.0, 6.0) << "  seconds"
                           << endl
                           << endl;

    } /// end if (key) processing

    cout << "-------------------------------------------------------------------------" << endl;
    Out_macrocrack_growth_stream << "-------------------------------------------------------------------------" << endl;
    cout << " START of the PCC Subcomplex module " << endl;
    Out_macrocrack_growth_stream << " START of the PCC Subcomplex module " << endl;
    std::vector<Subcomplex> pcc_subcomplexes; // vector containing all the PCC subcomplexes (cuts, k-order grain neighbours, etc)

    if (key_subcomplex == 1) {
            pcc_subcomplexes = PCC_Subcomplex(configuration, special_face_sequence);

    cout << " PCC_Subcomplexes size\t=\t\t" << pcc_subcomplexes.size() << endl;

    // ================ Elapsing time for the Subcomplex module ================
    unsigned int Subcomplex_time = clock();
    S_time = (double) Subcomplex_time - Main_time;
    cout << "Section time is equal to  " << S_time / pow(10.0, 6.0) << "  seconds" << endl;
    cout << "-------------------------------------------------------" << endl;
    Out_macrocrack_growth_stream << "Section time is equal to  " << S_time / pow(10.0, 6.0) << "  seconds" << endl;
    Out_macrocrack_growth_stream << "-------------------------------------------------------" << endl;

} /// end if(key) subcomplex

    cout << "Get_sfaces_sequence\t" << pcc_subcomplexes.at(0).Get_sub_sfaces_sequence().size() << endl;

if (key_multiphysics == 1) {
    cout << "-------------------------------------------------------------------------" << endl;
    Out_macrocrack_growth_stream << "-------------------------------------------------------------------------" << endl;
    cout << "START of the PCC Multiphysics module " << endl << endl;
    Out_macrocrack_growth_stream << "START of the PCC Multiphysics module " << endl << endl;

    crack_growth_series.clear();
    new_cells_energies = PCC_Multiphysics(configuration, pcc_subcomplexes, crack_growth_series);

    cout << "Get_crack_length" << crack_growth_series.back().Get_crack_length() << endl;
    cout << "Get_sfaces_sequence SIZE:\t\t" << crack_growth_series.back().Get_sfaces_sequence().size() << endl;
    cout << "Get_sfaces_coordinates SIZE:\t\t" << crack_growth_series.back().Get_sfaces_coordinates().size() << endl;
//    exit(0);
    cout <<  endl; Cracked_stat_out << endl;
    for(auto mcc : crack_growth_series) {
        cout << mcc.Get_crack_length() << "\t";
        Cracked_stat_out << mcc.Get_crack_length() << "\t";
    }
        cout << "\tp\t" <<  vop.at(0) <<  "\tmu\t" <<  vop.at(1) << "\tsigm\t" <<  vop.at(2) <<  "\tbins\t" <<  vop.at(3) << endl << endl;
        Cracked_stat_out  << "\tp\t" <<  vop.at(0) <<  "\tmu\t" <<  vop.at(1) << "\tsigm\t" <<  vop.at(2) <<  "\tbins\t" <<  vop.at(3) << endl << endl;

    std::vector <unsigned int> large_inclusions_seq;
    unsigned int full_strips_lengths = 0;
    Cracked_sfaces_numbers_in_plane_out << "inclusions_seq.size()" << endl;
    for(auto mcc : crack_growth_series) {

        large_inclusions_seq.clear();
        for(unsigned int sface : mcc.Get_sfaces_sequence()) {
// REPAIR            cout << "special_face_series " <<  special_face_series.size() << endl;

            for (auto strip : special_face_series) {
                if (strip.size() > 2 && strip.size() < 100 && std::find(strip.begin(), strip.end(), sface) != strip.end()) {
                    large_inclusions_seq.push_back(sface);
                }

            } // end for (auto strip : special_face_series)
        } // end for(unsigned int sface : mcc.Get_sfaces_sequence())

        //      for (auto mcsq : mcc.Get_sfaces_sequence()) {
        cout << "\tlarge_inclusions_seq\t" << large_inclusions_seq.size() << "\t";
        Cracked_sfaces_numbers_in_plane_out << large_inclusions_seq.size() << "\t";

    } // end for(auto mcc : crack_growth_series)

    cout <<  endl; Cracked_sfaces_numbers_in_plane_out << endl;

    unsigned int Yjunctions_number = 0;
    std::vector<Agglomeration> agglomerations_powder_vector;
    std::vector<std::vector<unsigned int>> agglomeration_faces_powers;
    Cracked_sfaces_numbers_in_plane_out << "\tfull_strips_lengths\t" << endl;
    for(auto mcc : crack_growth_series) {
        full_strips_lengths = 0;
        Yjunctions_number = 0;
        for(unsigned int sface : mcc.Get_sfaces_sequence()) {
// REPAIR            cout << "special_face_series " <<  special_face_series.size() << endl;

            for (auto strip : special_face_series) {
                if (std::find(strip.begin(), strip.end(), sface) != strip.end()) {
                    full_strips_lengths += strip.size();
///
//                        cout << "\tvop_counter\t" << vop_counter << endl;
                        agglomerations_powder_vector = pcc_processed.agglomerations_in_powders.at(vop_counter-1);
//                        cout << "\tagglomerations_powder_vector SIZE\t" << agglomerations_powder_vector.size() << endl;
                    agglomeration_faces_powers.clear();
                    agglomeration_faces_powers.resize(2);
                    for (auto apv : agglomerations_powder_vector) {
                        agglomeration_faces_powers[0].push_back(apv.Get_agglomeration_kcell_number());
                        agglomeration_faces_powers[1].push_back(apv.Get_agglomeration_power());
                    }
///                    cout << "\tagglomeration_faces_powers[0] SIZE\t" << agglomeration_faces_powers[0].size() << endl;
///                    cout << "\tagglomeration_faces_powers[1] SIZE\t" << agglomeration_faces_powers[1].size() << endl;

                    for (auto sf : strip) {
                        for (auto ifn = agglomeration_faces_powers[0].begin(); ifn != agglomeration_faces_powers[0].end(); ++ifn) { //                            cout << "\tagglomeration_faces_powers[0].at(distance(agglomeration_faces_powers[0].begin(),ifn))\t" << agglomeration_faces_powers[0].at(distance(agglomeration_faces_powers[0].begin(),ifn)) << endl; //                            cout << "\tdistance(agglomeration_faces_powers[1].begin(),ifn)\t" << distance(agglomeration_faces_powers[0].begin(),ifn) << endl; //                            cout << "\tagglomeration_faces_powers[1].at(distance(agglomeration_faces_powers[0].begin(),ifn))\t" << agglomeration_faces_powers[1].at(distance(agglomeration_faces_powers[0].begin(),ifn)) << endl;
                            if (agglomeration_faces_powers[0].at( distance(agglomeration_faces_powers[0].begin(), ifn)) == sf)
                                Yjunctions_number += agglomeration_faces_powers[1].at(distance(agglomeration_faces_powers[0].begin(), ifn));
//                            cout << "\tYjunctions_number\t" << Yjunctions_number << endl;
                        }
                    }// end  for (auto sf : strip)
                } // end if
            } // end for (auto strip : special_face_series)
        }// end for(unsigned int sface : mcc.Get_sfaces_sequence())
        Cracked_sfaces_numbers_in_plane_out << full_strips_lengths << "\t";
        Yj_stat_out << Yjunctions_number << "\t"; //      for (auto mcsq : mcc.Get_sfaces_sequence()) {
    } // end for(auto mcc : crack_growth_series)
    Yj_stat_out << endl;  Cracked_sfaces_numbers_in_plane_out << endl;
    cout <<  endl;

    for(auto mcc : crack_growth_series) {
        large_inclusions_seq.clear();
        full_strips_lengths = 0;
        for(unsigned int sface : mcc.Get_sfaces_sequence()) {
// REPAIR            cout << "special_face_series " <<  special_face_series.size() << endl;

            for (auto strip : special_face_series) {
                if (strip.size() >= 2 && std::find(strip.begin(), strip.end(), sface) != strip.end()) {
                    large_inclusions_seq.push_back(sface);
                }

                if (std::find(strip.begin(), strip.end(), sface) != strip.end()) {
                    full_strips_lengths += strip.size();
                }
            } // end for (auto strip : special_face_series)
        } // end for(unsigned int sface : mcc.Get_sfaces_sequence())

        cout << "\tratio\t" << large_inclusions_seq.size()/(double) full_strips_lengths << "\t";
        Cracked_stat_out << large_inclusions_seq.size()/(double) full_strips_lengths << "\t";

///            cout << mcc.Get_sfaces_sequence().size() << "\t";
///            Cracked_stat_out << mcc.Get_sfaces_sequence().size() << "\t";
    } // end for(auto mcc : crack_growth_series)

    cout <<  endl; Cracked_stat_out << endl;

//        std::vector<unsigned int> sfaces_sub_seq =  // pcc_subcomplexes.at(0). //Get_sfaces_sequence();
    std::vector<std::tuple<double, double, double>> internal_face_coordinates = crack_growth_series.back().plane_subcomplex.Get_sub_polytope_coordinates(); //kCell_barycentre_coordinates(2, sfaces_sub_seq);

//    cout << "Internal_face_coordinates SIZE:\t\t" << internal_face_coordinates.size() << endl;
//    cout << "Internal_sface_sequence SIZE:\t\t" << crack_growth_series.back().plane_subcomplex.Get_sfaces_sequence().size() << endl;

    for (auto acc: internal_face_coordinates) {
        cout << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << endl;
        Cracked_pcc_out << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t"
                        << endl;
    }
    cout << endl;
    Cracked_pcc_out << endl;

// REPAIR     Vector_of_vectors_ui_cout(macrocrack_sfaces_vector, "macrocrack_sfaces"s);
    // ================ Elapsing time for the Multiphysics module ================
    unsigned int Multiphysics_time = clock();
    M_time = (double) Multiphysics_time - S_time - Main_time;
    cout << endl << "Multiphysics time is equal to  " << M_time / pow(10.0, 6.0) << "  seconds" << endl
         << endl; //cout << "-------------------------------------------------------------------------" << endl;
    Out_macrocrack_growth_stream << endl << "Multiphysics time is equal to  " << M_time / pow(10.0, 6.0) << "  seconds"
                       << endl
                       << endl;

} /// end if(key) multyphysics

    } /// END for (auto vop: vector_of_powders) {
///    exit(0);
///=======================================================================================
///=======================================================================================

/////////////////////////////////////////////
        /// ====================== II. PCC Processing module ======================
        cout << "-------------------------------------------------------------------------" << endl;
        Out_macrocrack_growth_stream << "-------------------------------------------------------------------------" << endl;
        cout << "START of the PCC Processing INDUCED module " << endl;
        Out_macrocrack_growth_stream << "START of the PCC Processing INDUCED module " << endl;

    /// Studied material with inclusions
    Material Zr_NC_rGO("Zr_NC"s, "rGO"s);

    /// Loop over the all considered defect configurations:
        for (int powder_id = 0; powder_id < vector_of_powders.size(); ++powder_id) {
    cout << "pcc_processed.agglomerations_in_powders.at(powder_id) SIZE\t" << pcc_processed.agglomerations_in_powders.at(powder_id).size() << endl;
            std::vector<Agglomeration> agglomerations_vector = pcc_processed.agglomerations_in_powders.at(powder_id);
            std::vector<std::vector<double>> max_cfractions_vectors(4);
            max_cfractions_vectors[2].push_back({0});
            max_cfractions_vectors[2][0] = 0.1; ////// ARBITRARY NOW ///////

            induced_face_sequence.clear();
            /// For ALL including the first pass (!)
            std::vector<unsigned int> cState_n_vector(CellNumbs.at(0), 0), cState_e_vector(CellNumbs.at(1), 0), cState_f_vector(CellNumbs.at(2), 0), cState_p_vector(CellNumbs.at(3), 0);
            Configuration_cState = {cState_n_vector, cState_e_vector, cState_f_vector, cState_p_vector};

            /// MICROCRACKING ==================================
            /////// !!! 0.0
            double inclusion_r_intensity_factor = 1.0, microcrack_c_intensity_factor = 0.0;
            std::vector<vector<unsigned int>> ifs_vector;
    if(key_induced == 1) {
        ifs_vector.clear();
        pcc_processed.de_fractions_iface_vector.clear();
        for (int crack_id = 0; crack_id < new_cells_energies.size(); ++crack_id) {
            CellEnergies new_face_energy = new_cells_energies.at(crack_id);
            induced_face_sequence = PCC_Kinematic2_cracking(2, special_face_sequence, Configuration_cState,
                                                            max_cfractions_vectors, inclusion_r_intensity_factor,
                                                            microcrack_c_intensity_factor, new_face_energy,
                                                            agglomerations_vector, Zr_NC_rGO);
            //        std::vector <unsigned int> PCC_Kinematic2_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<std::vector<unsigned int>> &Configuration_cState, std::vector<std::vector<double>> const &max_cfractions_vectors, double inclusion_sintensity_factor, double crack_sintensity_factor, std::vector<CellEnergies> &new_cells_energies, Material &Mat_id) {
            ifs_vector.push_back(induced_face_sequence);

            /// de_fractions_vector
            pcc_processed.de_fractions_iface_vector.push_back(d_fractions_vector(EdgesTypesCalc(CellNumbs, induced_face_sequence, FES)));
        }

        /// Output to file
        Dspace_ifaces_out << "Powder #\t" << powder_iterator++ << endl;
        for (auto difrac : pcc_processed.de_fractions_iface_vector) {
            for (auto itdf : difrac) {
// REPAIR                cout << itdf << "\t";
                Dspace_ifaces_out << itdf << "\t";
            }
            Dspace_ifaces_out << "\t";
        }
        Dspace_ifaces_out << endl;

        std::vector<std::vector<std::tuple<double, double, double>>> microcrack_coordinates_tensor;
        microcrack_coordinates_tensor.clear();
        for (int crack_id = 0; crack_id < new_cells_energies.size(); ++crack_id)
///            microcrack_coordinates_tensor.push_back(kSequence_barycentre_coordinates(2, ifs_vector.at(crack_id))); // coordinates for each crack_id
            microcrack_coordinates_tensor.push_back(kFaceSeq_barycentre_coordinates(2, ifs_vector.at(crack_id))); // coordinates for each crack_id

        /// for (auto mc : microcrack_coordinates_vector) {
        std::tuple<double, double, double> mbeg, mback;
//        double x = get<0>(face_coordinates.at(itr)) * std::get<0>(sample_dimensions);
//        double y = get<1>(face_coordinates.at(itr)) * std::get<1>(sample_dimensions);
//        double z = get<2>(face_coordinates.at(itr)) * std::get<2>(sample_dimensions);

        for (unsigned int i = 0; i < microcrack_coordinates_tensor.begin()->size(); ++i) {
            mbeg = microcrack_coordinates_tensor[0][i];
            mback = microcrack_coordinates_tensor.back().at(i);
            Microcracks_out << get<0>(mbeg) * 10.0 << "\t" << get<1>(mbeg) * 10.0 << "\t" << get<2>(mbeg) * 10.0
                            << "\t\t" << get<0>(mback) * 10.0 << "\t" << get<1>(mback) * 10.0 << "\t"
                            << get<2>(mback) * 10.0 << endl;
//                    Microcracks_out << get<0>(mc) * 10.0 << "\t" << get<1>(mc) * 10.0 << "\t" << get<2>(mc) * 10.0 << endl;
        }
        Microcracks_out << endl << endl;

        /// Betti_Li_vector
        if (key_iBetti == 1) {
            std::vector<double> Betti_Li_microcrack_vector = OperatorsBetti(induced_face_sequence, ENS, FES, GFS);
            pcc_processed.Betti_0_iface.push_back(Betti_Li_microcrack_vector.at(0));
            pcc_processed.Betti_1_iface.push_back(Betti_Li_microcrack_vector.at(1));
            pcc_processed.Betti_2_iface.push_back(Betti_Li_microcrack_vector.at(2));
            pcc_processed.inverse_connectivity_iface.push_back(
                    std::log(Betti_Li_microcrack_vector.at(0) / Betti_Li_microcrack_vector.at(2)));
        }
        cout << endl << "\tCracked/induced_face_sequence SIZE\t" << induced_face_sequence.size() << endl << endl;
    }// if(key_induced == 1)

    } // end of for (int powder_id = 0; powder_id < vector_of_powders.size(); ++powder_id)

        // ================ Elapsing time for the Processing module ================
        unsigned int Processing_time = clock();
        P_time = (double) Processing_time - S_time - M_time - Main_time;
        cout << "Processing INDUCED time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl
             << endl; //cout << "-------------------------------------------------------------------------" << endl;
        Out_macrocrack_growth_stream << "Processing INDUCED time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl
                           << endl;

    Face_Edge_indices.close();
    Microcracks_out.close();
    Cracked_stat_out.close();
    Cracked_pcc_out.close();
    Dspace_ifaces_out.close();
    Cracked_sfaces_numbers_in_plane_out.close();
    Yj_stat_out.close();
////    exit(0);
///////////////////////////////////
/**
    std::vector<std::vector <unsigned int>> macrocrack_sfaces_vector;
    macrocrack_sfaces_vector.clear();
    int c_num = 0;
    for (Macrocrack cg_series : crack_growth_series) {
        macrocrack_sfaces_vector.push_back(cg_series.Get_sfaces_sequence()); // Macrocrack cg_series
    }
    pcc_processed.Set_macrocrack_sfaces(macrocrack_sfaces_vector);

    PCC current_PCC;
    current_PCC.Set_face_barycentre_coordinates();
    std::vector<std::tuple<double, double, double>> FBC = current_PCC.Get_face_barycentre_coordinates();

    std::ofstream Cracked_pcc_out, Cracked_stat_out;
    Cracked_pcc_out.open(output_dir + "Macrocrack_sfaces_coordinates.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Cracked_stat_out.open(output_dir + "Macrocrack_sfaces_fractions.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function

    std::vector <unsigned int>  sf_seq;
    std::vector<std::tuple<double, double, double>> sf_seq_coord;
    std::vector<double> sfaces_macrocrack_fractions;
    for(Macrocrack cgs : crack_growth_series) {
///        for(auto imc = crack_growth_series.begin(); imc != crack_growth_series.end(); ++imc) {
///            crack_growth_series.at(imc).Set_sfaces_sequence(special_face_sequence);
        cgs.Set_sfaces_sequence(special_face_sequence);
        sf_seq = cgs.sfaces_sequence;
// REPAIR  string vid = "macrocrack_sfaces"s; Vector_ui_cout(sf_seq, vid);
    sf_seq_coord.clear();
//            cout << "sf_seq SIZE\t" << sf_seq.size() << endl; cout << "FBC SIZE\t" << FBC.size() << endl;

        for (auto fsbc : face_sequence_barycentre_coordinates(sf_seq, FBC)) {
           if (get<2>(fsbc) <= cgs.Get_crack_length())
               sf_seq_coord.push_back(fsbc);
       }
//            cout << "sf_seq_coord SIZE\t" << sf_seq_coord.size() << endl;
        cgs.sfaces_coord = sf_seq_coord;

//           cout << "sf_seq_coord SIZE\t" << sf_seq_coord.size() << endl; cout << "Get_sfaces_coordinates SIZE\t" << cgs.sfaces_coord.size() << endl;
        cout << "\tCrack length ID:\t"<< c_num << endl; Cracked_pcc_out << "\tCrack length ID:\t"<< c_num++ << endl;
        for (auto acc : cgs.sfaces_coord) {
            cout << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << endl;
            Cracked_pcc_out << get<0>(acc) * 10.0 << "\t" << get<1>(acc) * 10.0 << "\t" << get<2>(acc) * 10.0 << "\t" << endl;
        }
        cout << endl; Cracked_pcc_out << endl;

        sfaces_macrocrack_fractions.push_back(cgs.sfaces_sequence.size()/cgs.plane_subcomplex.Get_sub_faces_set().size());

        cout << sfaces_macrocrack_fractions.back() << endl;
        Cracked_stat_out << sfaces_macrocrack_fractions.back() << endl;

    } // end of for(Macrocrack cgs : crack_growth_series)
    Cracked_pcc_out.close();
    Cracked_stat_out.close();
///        for(Macrocrack cms : crack_growth_series) { cout << "cms.Get_sfaces_coordinates() SIZE\t" << cms.sfaces_coord.size() << endl; }



//        void Macrocrack::Set_sfaces_sequence(std::vector <unsigned int> const &special_face_sequence) {
// REPAIR     Vector_of_vectors_ui_cout(macrocrack_sfaces_vector, "macrocrack_sfaces"s);

    std::vector<std::set <unsigned int>> macrocrack_sfaces_series;
    std::vector <unsigned int> macrocrack_sface;
    std::vector<std::vector <unsigned int>> macrocrack_sfaces; // = pcc_processed.Get_macrocrack_sfaces();
    /// Finding inclusions in the current powder

    for (Macrocrack mc : crack_growth_series) {// over the growth series of a macrocrack
        macrocrack_sface = mc.Get_sfaces_sequence();
        macrocrack_sfaces.push_back(macrocrack_sface);
    }

    pcc_processed.Set_macrocrack_sfaces(macrocrack_sfaces);
    **/

///    } /// END of the iteration loop over all POWDERS:     for (auto vop: vector_of_powders)

return pcc_processed;
} /// END of Macrocrack_growth(Config &configuration)
