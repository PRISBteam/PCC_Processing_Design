///================================ DCC Writer module ===================================================================================///
///=====================================================================================================================================///
/** The function in this library perform a structured output of the calculation results **/
///===================================================================================================================================///

/// Attached user defined C++ libraries:
///-------------------------------------
#include "Writer_Functions.h"
///-------------------------------------

void DCC_Writer(std::vector<unsigned int> const  &sfaces_sequence, std::vector<unsigned int> const  &kface_sequence) {

    /// Output special and ordinary face sequences to the output directory specified in config.txt
    DCC_sequences_Writer(sfaces_sequence, kface_sequence);

}