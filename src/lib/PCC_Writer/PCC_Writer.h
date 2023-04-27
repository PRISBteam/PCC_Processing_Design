///================================ DCC Writer module ===================================================================================///
///=====================================================================================================================================///
/** The function in this library perform a structured output of the calculation results **/
///===================================================================================================================================///

/// Attached user defined C++ libraries:
///-------------------------------------
#include "Writer_Functions.h"
//#include "DCC_Writer/Writer_Functions.h"
///-------------------------------------
void PCC_Writer(CellsDesign &new_cells_design, ProcessedComplex &pcc_processed) {
// Read PCC Writer specifications from the writer.ini file and output of the current configuration to the screen and .log file
std::vector<int> writer_specifications; // vector<int> containing writer specifications and formats
config_reader_writer(source_path, writer_specifications, Out_logfile_stream); // Read and output the initial configuration from the writer.ini file

int output_counter = 0; // special counter for output numeration

    /// Output special and ordinary face sequences to the output directory specified in config.txt
    PCC_CellSequences_Writer(new_cells_design, output_counter);
    PCC_Entropic_Writer(pcc_processed, output_counter);
    PCC_AnalyticalFractions_Writer(pcc_processed, output_counter); // analytical solutions

/*
       char* stype = const_cast<char*>(P_type.c_str()); // Processing type variable
        char P_simulation_type = *stype;

  if (P_simulation_type == 'L') {
        Strips_distribution_Writer(output_counter); // strips distribution
        RWseries_Writer(RW_series_vector, output_counter); // vector of all strips

        /// agglomerations  output
        vector<agglomeration> agglomerations_vector; // vector containing aglomarations as its objects - only for L type simulations (!)
        Agglomerations_Writer(output_counter, RW_series_vector, mu_f, sigm_f);

    } // end of if(P_simulation_type == 'L')
*/
     return ;
} // END of PCC Writer module

/// * Related functions * ///