///================================ A part of the PCC Writer module ======================================================================///
///======================================================================================================================================///
/** The library contains functions for the formatted output of subcomplexes associated with a PCC                                      **/
///===================================================================================================================================///

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "../../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../../PCC_Objects.h"
#include "../../ini/ini_readers.h"

using namespace std; // standard namespace

extern std::ofstream writer_logfile_stream;
extern std::string output_dir, source_path;
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::vector<std::string> paths_to_PCC_matrices; // The vector containing the paths to all the PCC's matrices, measures and other supplementary data files

#include "writer_subcomplexes.h"

void PCC_Subcomplex_Writer(std::vector<Subcomplex> &pcc_subcomplexes, int output_counter) {

    // Off-streams
    ofstream Out_subcomplex_max_cells; // List of vectors of a PCC's maximal cells formed a subcomplex each
    // File names and output directories
    string seq_subcomplex_odir = output_dir + "pcc_subcomplex_max_cells.txt"s; // output file

    // Output to file
    if (Out_subcomplex_max_cells) {
        /// Creation of the empty files
        Out_subcomplex_max_cells.open(seq_subcomplex_odir, ios::trunc);

        // special ASSIGNED cell sequences
            for (auto subcomplex : pcc_subcomplexes) {
                if (subcomplex.Get_sub_polytope_set().size() > 0) {
                    for (auto grain_id : subcomplex.Get_sub_polytope_set())
                        Out_subcomplex_max_cells << grain_id + 1 << " "; /// vit + 1 !!! for compatibility with the Neper output
                    Out_subcomplex_max_cells << endl;
                } // enf if (subcomplex.size() > 0)
            } // end for (auto subcomplex : .. )

    } else cout << "Error: No such a directory for\t" << seq_subcomplex_odir << endl;

    return;
} // END of PCC_Subcomplex_Writer()
