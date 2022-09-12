* =================================================================================================================================== *
* =========================	 DCC Structural Design tool	 ==================================================================== *
* =================================================================================================================================== *

Version: 1.0
Author: Dr Elijah Borodin (Research Fellow in Materials Physics at the University of Manchester; Mechanics and Physics of Solids research group)
Contacts: Elijah.Borodin@icloud.com
Release date: 17/09/2022

DCC_Structural_Design_tool ("SDT code") is the Materials Design tool intended to design defect microstructure evolution during material processing. The code provides an effective tool for the studies and design of both advanced material microstructure and processing routes.

Release notes:
The several sets of combinatorial and topological parameters can be set directly in the file design_request.txt. Using the Monte-Carlo type of stochastic simulations (Ising-like model), the code will create the structure of special elements arranged to fit all the parameters as close as possible (in their average) to the requested ones. 
… 
The code is written and tested in C++17. It works well with CMake 3.23 (cmake.org), g++ compiler (gcc.gnu.org) and CLion IDE (jetbrains.com/clion).

Acknowledgements:
This code has been written as a part of the EPSRC-funded project EP/V022687/1 "Patterns recognition inside shear bands: tailoring microstructure against localisation" (PRISB).

*  =========================================================== < < * > > =========================================================== *

The principal files, modules and concepts of the project:

The code consists (1) of the main_StrDesign.cpp where all functions launching and the data reading from files occur; three libraries (or modules) with related functions :
(2) DCC_Processing, 
(3) DCC_Characterisation,
(4) DCC_Writer.
These libraries contain all the corresponding functions, and only part of them (such as DCCprocessing.h, DCCcharacterisation.h, DCCwriter.h) was used in the SDT code.

*config.txt* file
The .txt file contains
1. *input* directory address, 
2. *output* directory address,
3. the number of calculation points for different values of *p*,
4. the precision of the Monte-Carlo algorithm,
and other parameters necessary for the code.

*design_requests.txt* file
The .txt file containing only the matrix of indices in the form
/==============================
N i(S) i(Sd) i(chi) i(sigma) ...
1 X X X X ...
2 X X X X ...
3 X X X X ...
......
/==============================
where each *X* can may contain any value between 0 and 1, and # are the numbers numerating different designs. 
(!) It is important that each design is applied to the whole process, which is the changing of the fraction of special faces *p* from 0 to 1.

Each of the *index* for any variable *A* is equal to

i(A) = A - A_min / A_max - A_min,

where *A* is the current value of the variable, A_min is its minimal possible value, and A_max is its maximum possible value. So, for any variable it’s index evenly distributed between 0 < i(A) < 1. The requirement of a minimum value corresponds to the i(A) = 0, while the maximum - i(A) = 1 instead of *X* in the request.txt file. 

*s2c_sequences*
Means “special 2-cells sequences”. It contains the list of numbers (the numeration is consistent with the numeration of faces (2-cells) elements) in the discrete cell complex (DCC) possessing the “special” IDs. 

Processing:
For the whole process of changing the fraction of special faces *p* from 0 to 1
1. Analyses the *request* matrix,
2. Using the Monte-Carlo-like random optimisation algorithm calculates *s2c_sequences* and saves them in the corresponding matrices,
3. Calculates the purely “random” *s2c_sequence* as the reference.

Characterisation:
1. In this code, this module just calculates the precision matrix of all the requested design parameters: the difference in per cent between the requested (in request.txt) and actual fitted value. If only a single parameter was optimised, the difference is typically minimal, but for several parameters, it can be impossible to fit well simultaneously all of them. 

Writer:
1. Writes the sets of the designed s2c_sequences as the matrix (each sequence as a separate *line*) in *s2c_design.txt* file alongside the random case sequence in a separate *s2c_random.txt* file.
2. Writes the requested design precision file *design_precision.txt* containing the matrix similar to the one in the *request.txt* file, but with the difference in per cent between the requested (in request.txt) and actual fitted value, instead of *X*. 

Practical tutorial:
*** to be continued..

* =================================================================================================================================== *
* ================================	End of file	============================================================================= *
* =================================================================================================================================== *

	
