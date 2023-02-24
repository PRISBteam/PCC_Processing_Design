<h1>Discrete Processing Design (DPD code)</h1>

Manual version: 0.3.0 <br>
First manual release date: 21/09/2022 <br>
Current manual release date: 15 January 2023 <br>

<p> DPD code is a <i> Materials Design tool</i> intended to design defect microstructure evolution during material processing. The code provides an effective tool for the studies and design of both material microstructures and the effect of the specific material processing routes. Its key feature is the usage of polyhedral cell complexes, which provide a discrete space for designing realistic material defect structures of different dimensions and types. </p>

<p> This is a C++ based software project consisting of several modules working with pre-created Polyhedral Cell Complex (PCC) as the set of its incidence and adjacency matrices represented in a sparse matrix form. The code contains 4 main modules (Processing, Characterisation and Writer), the library of Simulation tasks and the Main module merging everything together, reading configurations and data files,  and launching all the other modules.
The code works both with 3-complexes (3D) and 2-complexes (2D) BUT - to make results consistent with 2D/3D EBSD scans - it assumes that the grains are 3-cells in the 3D, 2-cells in the 2D case, and so on for grain boundaries and other element types. So it actually replaces definitions of k-cells with ( k + (dim - 3) )-cells, where dim = {2, 3} is the dimension of the problem. In such a way 2-cells are ALWAYS associated with grain boundaries on EBSD maps and are edges for 2D case! </p>
  
<h2>General specifications</h2>
  
<p>The code is written and tested in C++ 17 with the parallelised verson used some features of C++ 20. It is used explicitly Eigen and Spectra libraries which must be <a href="https://spectralib.org/download.html"> downloaded</a> and copied to the directory containing all the STL C++ libraries on the local PC.

The computational costs of different calculation types, functions and tasks are hugely different: for instance, the component analysis or Monte-Carlo simulations are a very time consuming procedures, while the random generation of special chains are fast.<\p>

<h2>Basic definitions</h2>
 
<p> 
<ol>
  <li>Polyhedral cell complex (PCC): </li>
  <li>k-Cells: </li>
  <li>Nodes: </li>
  <li>Edges: </li>
  <li>Faces: </li>
  <li>Volumes: </li>
  <li>cell number: </li>
  <li>cell fraction: </li>
  <li>Ordinary cells: including OCellsNumb or  OrdinaryCellNumbs; related variables New2CellNumb </li>
  <li>Special cells: including SCellsNumb or  SpecialCellNumbs </li>
  <li>sfaces_sequence</li>
  <li>cfaces_sequence</li>
  <li>induced topology of the defect structure: </li>
  <li>State_sVector: If there are several types of special k-cells, only pair of s_faces_sequence and State_sVector describing types for each of them </li>
  <li>element types: including NewFaceType</li>
  <li>Normal and reduced Incidence and Adjacency matrices</li>
  <li>Subcomplex: including Plane cut (a,b,c,D), (reduced (k-1)-complex)</li>
  <li>k-Skeleton: </li>
  <li>k-Chain: </li>
  <li>k-sChain: </li>
  <li>k-cChain: </li> 
</ol>
  
  </p>
  
<h2> More sources on the mathematics related to the Polyhedral cell complexes</h2>

An excellent simple introduction to the DCC with their various applications is given in the <a href="https://link.springer.com/book/10.1007/978-1-84996-290-2" target="_blank"> book </a> of Leo Grady and Jonathan Polimeni _“Discrete Calculus. Applied Analysis on Graphs for Computational Science. (2010)_

<h2>Modules</h2>
All the modules except the Main consist of the interface and the core parts.  The interfaces contains pre- and post-processing of data for this particular module, adapting input from the Main, and manage the function implementations according to the calculation types specified in the configuration file. 

<ol>
  <li>PCC_Main: include libraries, global variables, reading from files, launching the other modules and simulation tasks; </li>
  <li>PCC_Processing: assign chains of special k-cells. Essentially, the ultimate goal of the module is to create a list of k-chains (or k-sequences) in the order of appearance of new elements during the consideration process; </li>
  <li>PCC_Subcomplex: including Plane cut (a,b,c,D), (reduced (k-1)-complex)</li>
  <li>PCC_Characterisation: performs characterisation of special structures; </li>
  <li>PCC_Writer: output data in files; </li>
  <li>Objects.h </li>
  <li>Support Functions </li> 
</ol>

Main

I. SUBCOMPLEX

II. PROCESSING

Aim: Set new structures (chains) on a PCC elements.

Several <i>generation types</i> deigned for various simulation tasks are available:

One-element generation functions
<ol>
<li> Random choice of a single element </li>
	- Use advanced algorithm  (C++ …. Library) for random number generation;
	- Used in all the other more complicated structure generation processes.
<li> Random Walker with leaps </li>
	- The random walk algorithm choosing travelling across the any k-skeleton (for any k) with some possibility of larger-scale leaps;
	-   Every new regular (without leap) step it choses between the neighbouring k-Cells using the Random choice function( );  
 </ol>

Elements’ chain generation functions
<ol>
<li> Random one by one choice of special elements </li>
	- Currently allow to set several types of special faces that encoded in the couple of sfaces_sequence (showing the generation order and allowing then easily change the fraction of special faces) and  State_Vector  (showing the type for each of the faces);
	- As its output generates set of vectors of  s<element>_sequences and  corresponding State_<element>Vector 
	- In a binary case (0 and 1 types only) sfaces_sequence is enough for the whole representation of special chains;
<li> Strip distributions </li>
<li> Max/Min configuration entropy production principle </li>
<li> Metropolis-like algorithm with the "energy" minimisation (design module) </li>
<li> Ising-like model </li>
</ol>

III. CHARACTERISATION

IV. WRITER
  
Objects and Support functions.h

<h2>External files</h2>
<ul>	
<li> Configuration file </li>
<li> PCC sparse matrices and their reading </li>  
<li> Input and output folders </li>
</ul>
	
<h2> Where to take a complex? </h2>
The discrete cell complex is a pretty well-known object that originated from the field of algebraic topology, so it can be obtained in many various ways Below is just a concise review of a couple of flexible tools developed in the Mechanics and Physics of Solids research group in the University of Manchester providing DCCs based on Voronoi and a few others tessellations of space by convex polygons. 

<h3> Tessellations of space provided by Neper software </h3>

The Voronoi tesselation provided by Neper supposed to be a <i>dual</i> complex and so all the other tessellations provided by the Neper output with the <a href="https://neper.info/doc/neper_t.html#morphology-options" target="_blank"> morphology </a> option <i> -morpho <morphology> </i> like <i> cube, square, tocta, lamellar, etc. </i> different from <i>voronoi</i>.

Please, see more <a href="https://neper.info/doc/neper_t.html#examples" target="_blank"> examples </a> on the Neper webpage.


<h3> DCC Generator Tool </h3>

Based on the Poisson-Voronoi tessellation of space provided by the <a href="https://neper.info" target=”_blank”> Neper </a> software the code creates discrete (combinatorial) cell complex (DCC) as the set of sparse matrices. The  DCC Generator Tool generates a sparse representation of matrices: for any matrix element _a_(_i_, _j_) = _c_, the files of the matrices contain the list of triplets in the form (_i_, _j_, _c_). The enumeration of indices starts from 0, and, for instance, the line "5, 7, 1" in the adjacency matrix A<sub>k</sub> means that the _k_-cell #6 is the neighbour of the _k_-cell #8. For any incidence matrices B<sub>k</sub>,  the same triplet "5, 7, 1" means that the (_k_-1)-cell #6 is on the boundary of the _k_-cell #8, and their orientations coincide (_c_ = -1 for the opposite orientations). 

All the other information on the GitHub page of the <a href="https://github.com/PRISBteam/Voronoi_DCC_Analyser/blob/main/README.md" target=”_blank”> project </a>
The latest release of the code can be downloaded from the <a href="[https://github.com/PRISBteam/Voronoi_DCC_Analyser/blob/main/README.md](https://github.com/PRISBteam/Voronoi_DCC_Analyser/tags)" target=”_blank”> DCGT </a> project page.


<h3> FCC and BCC primal slip planes </h3>

The package DCC_Structure contains Python modules to build a discrete cell complex (DCC) based on the slip planes of crystal nanostructures (simple cubic, FCC, BCC; HCP not yet available). The script execute.py is a summarised execution of the whole package. It takes as input:
…

<h2> DSD code Tutorial </h2>


<h3> 1. Input files and the requested design </h3>

All the input files must be in a single folder specified as the ‘input’ directory in the ‘config.txt’ file.

The several sets of combinatorial and topological parameters can be set directly in the file ‘design_request.txt’. Using the Monte-Carlo type of stochastic simulations (Ising-like model), the code creates the chains of special cells arranged to fit all the parameters as close to the requested ones as it is possible in their average. 

DCC itself as the set of sparse matrices..

‘config.txt’ file contains
<ol>
<li> ‘input’ directory address </li>
<li>  ‘output’ directory address </li>
<li>   the number of calculation points for different values of <i>p</i></li>
<li>   the precision of the Monte-Carlo algorithm,
and other parameters necessary for the code</li>

‘design_requests.txt’ file containing only the matrix of indices in the form {Nd i(S) i(Sd) i(chi) i(sigma)} where each term can contain any value between 0 and 1, and Nd are the numbers numerating different designs. 
Each of the *index* for any variable *A* is equal to
i(A) = A - A_min / A_max - A_min,
where *A* is the current value of the variable, A_min is its minimal possible value, and A_max is its maximum possible value. So, for any variable it’s index evenly distributed between 0 < i(A) < 1. The requirement of a minimum value corresponds to the i(A) = 0, while the maximum - i(A) = 1 instead of *X* in the request.txt file. 
  
It is important to mention that each design is applied to the whole process, which is the changing of the fraction of special faces *p* from 0 to 1.

Sparse matrices of the DCC and their names

<h3> 2. Compilation </h3>

What needs to be installed..
Command line commands for compilation of the code..

The code is written and tested in C++17. It works well with CMake 3.23 (cmake.org), g++ compiler (gcc.gnu.org) and CLion IDE (jetbrains.com/clion).


<h3> 3. Execution </h3>

Command line commands to execute the code..


<h3> 4. Calculation process </h3>

What is happening during the calculation process.. 

<h3> 5. Output files </h3>

The only output of the code is ‘s2c_sequences.txt’ file, which means <i> “special 2-cells sequences”</I>. It contains the list of numbers (the numeration is consistent with the numeration of faces (2-cells) elements) in the discrete cell complex (DCC) possessing the “special” IDs. 

In the order of appearance during the considering process.. 

<h2> Project internal structure </h2>

The principal files, modules and concepts of the project are
<ol>
<li> the ‘main.cpp’ where all functions launching and the data reading from files occur </li>
and three libraries (or modules) with related functions:
<li> DCC_Processing </li>
<li> DCC_Characterisation </li>
<li>  DCC_Writer </li>
</ol>

These libraries are part of the larger software development and contain a large set of functions, only a small part of which (such as DCCprocessing.h, DCCcharacterisation.h, DCCwriter.h) is used in the DSD code.


<h3>1. Processing module </h3>

For the whole process of changing the fraction of special faces *p* from 0 to 1
1. Analyses the *request* matrix,
2. Using the Monte-Carlo-like random optimisation algorithm, calculates *s2c_sequences* and saves them in the corresponding matrices,
3. Calculates the purely “random” *s2c_sequence* as the reference.


<h3>2. Characterisation module </h3>

1. In this code, this module just calculates the precision matrix of all the requested design parameters: the difference in per cent between the requested (in request.txt) and actual fitted value. If only a single parameter was optimised, the difference is typically minimal, but for several parameters, it can be impossible to fit well simultaneously all of them. 


<h3>2. Writer module </h3>

1. Writes the sets of the designed s2c_sequences as the matrix (each sequence as a separate *line*) in *s2c_design.txt* file alongside the random case sequence in a separate *s2c_random.txt* file.
2. Writes the requested design precision file *design_precision.txt* containing the matrix similar to the one in the *request.txt* file, but with the difference in per cent between the requested (in request.txt) and actual fitted value, instead of *X*. 


<h2> After the Structural Design tool </h2>


The DSD tool is typically just a basic tool for obtaining specific defect configurations on a given DCC. The next two obvious steps are _characterisation_ (which means combinatorial and topological analysis) of the defect structure and simulations of various physical and mechanical processes in the materials possessed with this microstructure.

As it was already mentioned, the design of the processing routes leading to a given microstructure development is of particular interest here. Such a tool provides a unique opportunity for detailed studies of the physical processes supplementing the microstructure evolution.

<h2> Applications of DSD tool </h2>

<ol>
<li> E.N. Borodin, A.P. Jivkov, A.G. Sheinerman, M.Yu. Gutkin, 2021. Optimisation of rGO-enriched nanoceramics by combinatorial analysis. Materials & Design 212, 110191. [doi: 10.1016/j.matdes.2021.110191.](https://doi.org/10.1016/j.matdes.2021.110191) </li>
</ol>

<h2> Acknowledgements </h2>

This code has been created as a part of the EPSRC funded projects EP/V022687/1 _“Patterns recognition inside shear bands: tailoring microstructure against localisation”_ (PRISB) and EP/N026136/1 _"Geometric Mechanics of Solids: a new analysis of modern engineering materials"_ (GEMS).


<h2> License </h2>

Distributed under the GNU General Public License v3.0. See LICENSE.txt for more information.


<h2> Contacts </h2>
Dr Elijah Borodin (Research Fellow in Materials Physics at the University of Manchester; Mechanics and Physics of Solids research group)
<a href=“ Elijah.Borodin@icloud.com” Send e-mail> to Elijah Borodin (any queries regarding the code) 

<!-- FUTURE DIRECTIONS:
v. +1 More types of "special" faces (like the fractured ones)
v. +1 More types of "special" not 2-elements (like grains or nodes) -->
