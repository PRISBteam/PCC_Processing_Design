#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>
#include <tuple>

using namespace std; //Standard namespace

/* Various useful functions (it must be here - first in this list! ) */
#include "../lib/PCC_SupportFunctions.h"

/* Various set measures */
#include "../lib/measures.h"

/* Objects library contains classes of various objects related to PCC substructures and elements (it must be here - second in this list! ) */
#include "../lib/PCC_Objects.h"

#include "../lib/PCC_Processing/PCC_Processing.h"

void Grain_design(CellsDesign &new_cells_design, std::vector<vector<int>> &Configuration_State, std::vector<vector<int>> &Configuration_cState) {
    cout << " START of the PCC Processing module " << endl;
    new_cells_design = PCC_Processing(Configuration_State, Configuration_cState);
}