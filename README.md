<h1>Polyhedral Cell Complex (PCC) Processing Design  (DPD code)</h1>

Manual version: 0.3.0 <br>
Current manual release date: 17 August 2023 <br>

<p>  <i> PCC Processing Design </i>  or   <i> Discrete Processing Design  (DPD code)</i> is a software aimed to provide a wholly discrete representation of defect microstructure evolution during the materials processing. The code is an effective tool both for studies microstructure network evolution during the materials processing and for design of advanced microstructures in composite polycrystalline materials. Its key feature is the usage of polyhedral cell complexes (PCCs) as objects of algebraic topology, which provide a convenient discrete space for large-scale design of realistic material defect microstructures of different dimensions (line defects, surface defects, different volumetric phases). Such PCCs can be created based on the 2D/3D space tessellation with the convex polyhedron’s such as the Voronoi tessellation of space realistically imitating a polycrystalline material.

An excellent simple introduction to the cell complexes with their various applications is given in the <a href="https://link.springer.com/book/10.1007/978-1-84996-290-2" target="_blank"> book </a> of Leo Grady and Jonathan Polimeni <i>"Discrete Calculus. Applied Analysis on Graphs for Computational Science" (2010).</i> 
More rigorous mathematically introduction to the main concepts and methods of combinatorial algebraic topology is given in the <a href="https://doi.org/10.1007/978-3-540-71962-5" target="_blank"> book </a> of Dmitry Kozlov <i> "Combinatorial Algebraic Topology" (2008) </i>. An excellent classical introduction to the graph theory is provided in the classical <a href="https://doi.org/10.1007/978-1-4612-0619-4" target="_blank"> book </a> of Béla Bollobás <i> "Modern Graph Theory" (1998)</i>.

For effective use of the DPD code, a clear understanding of several basic concepts, such as Polyhedral Cell Complex, state vectors and special sequences of cells is necessary. All these concepts (and much more) are described in the Theoretical Manual. Even more comprehensive discussions of some particular topics with examples can be found in the following publications:
<ol>
<li> S. Zhu, E.N. Borodin, A.P. Jivkov, 2021. Triple junctions network as the key structure for characterisation of SPD processed copper alloys. Materials & Design 198(24), 109352. doi: 10.1016/j.matdes.2020.109352. </li>
<li> E. N. Borodin, A. P. Jivkov, 2019. Evolution of triple junctions’ network during severe plastic deformation of copper alloys – a discrete stochastic modelling. Philosophical Magazine 100(4), 467-485. doi: 10.1080/14786435.2019.1695071. </li>
<li> E.N. Borodin, A.P. Jivkov, A.G. Sheinerman, M.Yu. Gutkin, 2021. Optimisation of rGO-enriched nanoceramics by combinatorial analysis. Materials & Design 212, 110191. doi: 10.1016/j.matdes.2021.110191. </li>
</ol>
In this Technical Guidance, only the architecture of the code and practical guidance for its effective use with examples are provided. Much other related stuff, including already created PCCs with the corresponding discrete space tessellations are presented on the  <a href="http://materia.team" target="_blank"> Materia</a> software project. 
</p>

<h2>A few necessary definitions</h2>
<p> All the k-cells with a specific value of k form a PCC <i>skeleton</i>, so that the whole PCC can be thought as an agglomeration of these several skeletons. </br>
Additionally to the scalar and vector parameters defined on the PCC cells, for each cell can be assigned a label characterising its own type. Hence all k-cells can be divided in two classes of <i>ordinary</i> (without labels or with the labels 0) and <i>special</i> (labelled) cells. </br>
Inside the class of special cells several types can be defined. The whole set of these labels assigned on the each skeleton form a configuration <i>state</i> which can be expressed by a single <i>state vector</i> encoding all the cell types like the well-known DNA sequences. </br>
The set of all the state-vectors in a PCC form its <i>structure</i>.  </br>
The notion of a <i>process</i> can be referred to as the discrete sequence of all the states expressed by their state vectors. 
Three different types of the structures can be distinguished: directly <i>assigned</i>, <i>imposed</i> by the cells of other dimensions and <i>induced</i> by a pre-existing assigned structure as a result of a kinetic process. </br>
The last component of the proposed approach is a function or <i>map</i> relating structural (such as fractions of special grain boundaries) and physical (time, plastic strain, mass fraction, etc. ) characteristics which can be  obtained based on the particular experimental data. </br>
Such a theoretical framework allows to describe and analyse space ordering inside a real materials by the well-developed tools of topology, statistics and graph theory.   
</p>

<h2> How to execute the code </h2>
<p> This is a C++ based software project consisting of several modules-libraries and the main.cpp file which is considered as a separate module. The code works with pre-created PCC as the set of its incidence and adjacency matrices represented in a sparse matrix form. It is intended to be launched as a CMake project using CMakeLists.txt file.
For successful code execution, both the C++ compiler and CMake software must be installed. 
 
The code is written and tested in C++17. It works well with <a href="http://cmake.org" target="_blank"> CMake 3.23 </a>, <a href="http://gcc.gnu.org" target="_blank"> g++ compiler </a> and <a href="http://jetbrains.com/clion" target="_blank"> CLion IDE </a>. It is partly parallelised with the <a href="https://www.openmp.org" target="_blank"> Open MP </a> libraries for its effective execution using simultaneously several cores within one CPU. 
It is used explicitly the Eigen and Spectra libraries which must be <a href="http://spectralib.org/download.html"> downloaded</a> and copied to the directory containing all the STL C++ libraries on the local PC.
As an example, to compile the project in a command line tool one needs to change the working directory (cd command in Linux) to the one containing the project CMakeLists.txt file and then launch CMake as:
```
cmake -B buildtree
cmake --build buildtree
```
With the CLion and other IDEs, everything is straightforward: a new C++ project must be created (if it contains its own main.cpp “Hello world!” file by default, it must be deleted or ignored) and then executed. 
The computational costs of different calculation types, functions and tasks are hugely different: for instance, the component analysis or Monte-Carlo simulations are a very time consuming procedures, while the random generation of special chains are fast. <br>
The code works equally good with 3D and 2D tessellations of space. In the 2D case, there are no 3D polyhedrons (volumes) and 2D polytopes are associated with faces or 2-cells of the corresponding PCC. All the project functions works similarly in these two cases.
</p>
  
<h2>Project files</h2>
<p>
The project directory contains several folders: <br>
<ul>
<li> <i>\src</i> — contains all the source files, *.h libraries of the project and the main.cpp file.</li>
<li> <i>\config</i> — contains all the *.ini files using for the initial definition of parameters governing the execution of the corresponding modules.</li>
<li> <i>\PCC_examples</i> — contains a few examples of the PCCs created based on the 2D and 3D Voronoi tessellations. A large PCC library is published on the webpage of the <a href="http://materia.team" target="_blank"> Materia</a> software project.</li>
<li> <i>\data</i> — contains supplementary data files for each of the project modules.</li>
</ul>

The <i>\src</i> "source" directory contains several subfolders with libraries:
<ul>
<li> <i>\lib</i> — contains all the *.h project libraries Each of the modules is placed in the folder with the same name as its own. These subdirectories contains, in their turn, *.h files with the same name as the corresponding module, and a subdirectory named \functions containing *.h libraries with all the functions used in this particular module. Besides, it contains the library <i>SupportFunctions.h</i> with the additional functions used in several modules, <i>measures.h</i> with functions for calculations of various structural measures, and another library <i>Objects.h</i> as the only place in the project containing all the definitions of the code-specific classes. The subdirectory <i>\ini</i> contains the code-specific readers of the project's *.ini files which nedds to be amended at their any significant changes.</li>
<li> <i>\task</i> — contains all the additional user-defined tasks written as separate functions and included in the main.cpp file. These functions became active in the mode "TASK" instead of the "LIST" of the [simulation_mode] option in the <i> main.ini</i> file. </li>
<li> <i>\other</i> — contains additional external libraries such as the simple <a href="https://github.com/pulzed/mINI"> mINI </a> reader of the *.ini files.</li>
</ul>

<h2> Project modules by their *.ini files </h2>
The main part of the “user interface” contains a few *.ini files governing the behaviour of each of the modules. Among them:
<ol>
<li> Main </li>
Provides a common “environment” for all the other project parts. Currently, it is the only .cpp file compiling in the project (see CMakeLists.txt), including all the other project libraries. All the variables, names and paths are defined here and then it calls all the other modules marked “ON” in the main.ini file.
<li> Processing </li>
The central part of the code - generates labelling of a PCC k-cells of different dimensions according to some governing principles. 
<li> Characterisation </li>
Takes as its input the state_vectors (for assigned, imposed and induced types of faces) as its input and performs different characterisation tasks to output an object of the C++ class Processed_Complex containing all structural characteristics of the special cell structures requested in the characterisation.ini file. The syntax of the characterisation.ini file is especially simple: it contains only boolean type parameters with a value equal to 1 means calculation of this particular characteristic, and 0 means that the code will not calculate it. 
<li> Writer </li>
Takes as its input the object Processed_Complex and performs only the output of various characteristics contained in this object to the “output” directory specified in the main.ini file with the pre-defined names. The writer.ini file also contains only boolean type of parameters with a value equal to 1 means writing these characteristics to the corresponding file, and 0 means that the parameter will not be written. 
</ol>

Some other modules like Multiphysics (stress and energies), Subcomplex (a subdivision of a PCC into subcomplexes), and Kinetic (processes related to state variables and do not determined by changes in the cell types of a PCC) still unfinished and are under active development. Each of the modules is a *.h library of the corresponding functions written in C++ and, in their turn, consists of two sub-modules - one governing the code implementation based on the requests written in the corresponding *.ini file, and another one (in the corresponding \functions subdirectory) contains the functions themselves. A couple of other *.h function libraries in the <i>\lib</i> directory like <i> SupportFunctions.h</i> contain many useful functions using in different modules for supplementary simple tasks like reading from files to matrices. For reading *.ini files the code exploits external <a href="https://github.com/pulzed/mINI"> mINI </a> library. 
</p>

<h3> main.ini </h3>
[general] <br>
<ul>
<li> “dim = 3” or 2 - is the problem dimension for 3D or 2D space tessellations, respectively. </li>
<li> “source = … \“ set the path to the directory containing PCC in its algebraic representation as a set of all adjacency and incidence matrices with some additional data about the corresponding space tessellation such as polyhedra volumes, face areas, face normals, etc. </li>
It is important to use “\” symbol at the end of the source path!
<li> “output=…\” set the output directory for the Writer module - where all the calculation results will be written. </li>
It is important to use “\” symbol at the end of the output path!
</ul>

[modules] <br>
All the rest in the main.ini file is only the list of all MODULES with the two variants: “ON” - for switching on the module execution, and “OFF” - for switching off the module execution. <br>
There are two different [simulation_mode] of the Main module: <br>
”LIST” - by default, launch all the modules one after another strictly according to the data from *.ini files. <br>
“TASK” - assumes tailored execution of the code using the functions, *.cpp and *.h files included explicitly inside the else if(task) {..} statement in the main.cpp module INSTEAD of the “LIST” mode sequence of modules. The ”task” mode is supposed to provide scientific freedom of the code execution and can ignore any instructions listed in the *.ini files. <br>

Each of the following *.ini files contains an almost similar list of settings for every type of cell in the tessellation and the corresponding PCC: polyhedrons (3-cells), faces (2-cells), edges (1-cells), nodes (0-cells). In the Processing module, any algorithm calculates as its output the lists (vectors) of “special” cells of different types described in the corresponding “state vectors”. The are three distinct sub-modules: (1) assigned structures: the algorithm picks cells and assigns them some type ID (label), writing the cell number in the corresponding special_sequence (s_sequence) vector (example: the random assignment of “special” type for some number of faces); (2) imposed structures: assigned types for low-dimensional (k-1)-cells or higher- dimensional (k+1)-cells according to some specific rule based on the already created assigned structures for k-cells (example: classification of face junctions according to the number of special faces incident to each junction); (3) induced structures: assigned types for the k-cells of the same dimension based on the already created assigned structures for k-cells (example: introducing fractured or cracked faces based on the initially assigned structure of faces containing inclusions). 

<h3> processing.ini </h3>
The file is divided in several parts reflecting the dimensions of the cells: <br>
[polyhedrons] <br>
Containing only instructions for the assignment of the polyhedrons (3-cells) types.
The set of parameters for polyhedra is in full analogy (possibly less in their number) with one for [faces]. Please read the detailed description below. <br>
[faces] <br>
Containing only instructions for the assignment of the faces (2-cells) types. Several parameters below define the settings for the Assignment type of Processing module.<br>
<ul>
<li> “face_types_number =..” — the number of distinct face types (normally from 0 to 3) where 0 means that there are no special faces and the module does nothing here. </li>
<li> “pf_mode = ..” — (if face_types_number >0) choose the specific processing type from the list of functions in the /src/lib/PCC_Processing/functions directory. </li>
S — reading from source *.txt file (s_sequence.txt ) the list of special faces (for this particular PCC) created before by some of the processing modes listed below;
<li> “source = /…/s_cells_sequence.txt” — the path to the *.txt file containing a list of numbers of faces of special types. This “source” affects only S processing mode and does not affect any other parts of the code. </li>
R — simple random choice of new special faces during the assignment process;
F — choice of new special faces governing by the maximum configuration entropy production principle (MEPP);
D — choice of new special faces governing by the minimum configuration entropy production principle;
Cr — determination of new special faces by effective random rotations of grains (applicable only for crystallography-related problems);  
Cm — determination of new special faces by effective rotations of grains governing by the minimum configuration entropy production principle (MEPP) (applicable only for crystallography-related problems);  
L — the random choice of new special faces with some restrictions that allow to the creation of elongated chains of special cells, whose lengths are normally distributed with the average “mu” and dispersion “sigma”. 
<li> “pf_index = 0” — supplementary index for more flexibility in the code execution, it does not affect anything in the default mode; </li>
The fractions from 0 to 1 for three possible face types specified above in the face_types_number parameter in their order. It is the fractions of special faces which will be assigned by the PCC_Processing module. If there is only one special type, only fmax_fraction1 should be above 0,  face_types_number =1, and all the rest fractions will be ignored. 
<li> fmax_fraction1 = 0.9 </li>
<li> fmax_fraction2 = 0.0 </li>
<li> fmax_fraction3 = 0.0 </li>
The following statements set parameters for the Induced type of Processing module and for “historical” reasons called cracked faces. By default, there is a possibility to set only one type of such induced faces. This part of the Processing module is always following (in execution time) after assignment one, and, by definition, uses the assignment face types for calculation of the corresponding list of induced faces type. 
<li> “crack_types_number = ..” — similarly to the face_types_number set the number of types and currently only two options are allowed: 0 - there are no induced faces, and 1 - means that the induced part of the Processing module is “on”. </li>
<li> “cf_mode = Km” — similarly to the pf_mode set the specific mode of the choice of induced faces. </li>
Km — currently only one mode of the “kinematic fracture” is allowed. 
<li> “cfmax_fraction = …” — similarly to the fmax_fraction parameters for assigned face types it sets the fraction of induced faces (in the range from 0 to 1). </li>
</ul>

[edges]<br>
Containing only instructions for the assignment of the edges (1-cells) types.
The set of parameters for edges is in full analogy (possibly less in their number) with one for [faces]. Please read the detailed description above. <br>
[nodes]<br>
Containing only instructions for the assignment of the nodes (0-cells) types.
The set of parameters for nodes is in full analogy (possibly less in their number) with one for [faces]. Please read the detailed description above. 
Finally, the “distribution” is a very special category relevant for the only case of the elongated chains of special cells, whose lengths are normally distributed with the average “mu” and dispersion “sigma”. It does not affect any other processing modes. <br>
[distribution]<br>
mu = 1.0<br>
sigma = 0.0

<h3>characterisation.ini</h3> 
The characterisation module is divided in several "labs" [polyhedrons_lab], [faces_lab], [edges_lab], [nodes_lab], corresponding to the each type of k-cells in a 3-complex (PCC) with the similar set of structural characteristics: <br>
"*l_active = .." - switch on/off the calculations of the characteristics related to the corresponding k-cells (polyhedra, feaces, edges, nodes); <br>
"config_entropy=.." - calculation of the configuration entropy with its mean (if "S_mean = 1") and deviatoric (if "S_skew = 1") parts.<br>
For the [edges_lab] in addition the calculation of the imposed by faces characteristics of "j_fractions = .." (special edge fractions), "d_fractions = .." (special edge degree fractions), and "analytical = .." (analytical solutions) is posdible. <br>
Finally, [spectra_lab] contains parameters of the corresponding Laplacian spectra. Here<br>
"calc_steps_numb = .." - set the number of points where the spectrum will be calculated;  "laplacians = .." - switch on/off the calculation of the calculation of the corresponding matrix of the combinatorial Laplacian; "laplacians_spectra = .." - switch on/off the calculation of the Laplacian's spectrum (the list of all its eigenvalues); "laplacians_betti = .." - switch on/off the Laplacian's Betti numbers as the dimensions of its null-space (the number of zero eigenvalues).

<h3>writer.ini</h3> 

[sequences]
isSequencesOutput = 1
isDesignvectorsOutput = 1

[entropic_polyhedrons]
isPolyhedronFractions = 0

[entropic_faces]
isFaceFractions = 0
isConfEntropy = 0

[entropic_edges]
isConfEntropy = 1
isFractions = 1
isDegreeFractions = 1

[entropic_nodes]
isNodeFractions = 0

[entropic_analytical]
isEdgeFractions = 0
isEdgeConfEntropies = 0

[component_analysis]
isBetti = 0

<h2> Where to take a complex? </h2>
The discrete cell complex is a pretty well-known object that originated from the field of algebraic topology, so it can be obtained in many various ways Below is just a concise review of a couple of flexible tools developed in the Mechanics and Physics of Solids research group in the University of Manchester providing DCCs based on Voronoi and a few others tessellations of space by convex polygons. 

<h3> Tessellations of space provided by Neper software </h3>

The Voronoi tesselation provided by Neper supposed to be a <i>dual</i> complex and so all the other tessellations provided by the Neper output with the <a href="https://neper.info/doc/neper_t.html#morphology-options" target="_blank"> morphology </a> option <i> -morpho <morphology> </i> like <i> cube, square, tocta, lamellar, etc. </i> different from <i>voronoi</i>.

Please, see more <a href="https://neper.info/doc/neper_t.html#examples" target="_blank"> examples </a> on the Neper webpage.


<h3> PCC Generator Tool </h3>

Based on the Poisson-Voronoi tessellation of space provided by the <a href="https://neper.info" target=”_blank”> Neper </a> software the code creates discrete (combinatorial) cell complex (DCC) as the set of sparse matrices. The  DCC Generator Tool generates a sparse representation of matrices: for any matrix element _a_(_i_, _j_) = _c_, the files of the matrices contain the list of triplets in the form (_i_, _j_, _c_). The enumeration of indices starts from 0, and, for instance, the line "5, 7, 1" in the adjacency matrix A<sub>k</sub> means that the _k_-cell #6 is the neighbour of the _k_-cell #8. For any incidence matrices B<sub>k</sub>,  the same triplet "5, 7, 1" means that the (_k_-1)-cell #6 is on the boundary of the _k_-cell #8, and their orientations coincide (_c_ = -1 for the opposite orientations). 

All the other information on the GitHub page of the <a href="https://github.com/PRISBteam/Voronoi_DCC_Analyser/blob/main/README.md" target=”_blank”> project </a>
The latest release of the code can be downloaded from the <a href="[https://github.com/PRISBteam/Voronoi_DCC_Analyser/blob/main/README.md](https://github.com/PRISBteam/Voronoi_DCC_Analyser/tags)" target=”_blank”> DCGT </a> project page.


<h3> FCC and BCC primal slip planes </h3>

The package DCC_Structure contains Python modules to build a discrete cell complex (DCC) based on the slip planes of crystal nanostructures (simple cubic, FCC, BCC; HCP not yet available). The script execute.py is a summarised execution of the whole package. It takes as input:
…

<h2> Applications of DPD code </h2>
<ol>
<li> S. Zhu, E. Borodin, A. P. Jivkov, Topological phase transitions of grain boundary networks during severe plastic deformations of copper alloys. Acta Materialia (Under review), 2023.</li>
<li> E.N. Borodin, A.P. Jivkov, A.G. Sheinerman, M.Yu. Gutkin, 2021. Optimisation of rGO-enriched nanoceramics by combinatorial analysis. Materials & Design 212, 110191. [doi: 10.1016/j.matdes.2021.110191.](https://doi.org/10.1016/j.matdes.2021.110191) </li>
</ol>

<h2> Acknowledgements </h2>

This code has been created as a part of the EPSRC funded projects EP/V022687/1 _“Patterns recognition inside shear bands: tailoring microstructure against localisation”_ (PRISB) and EP/N026136/1 _"Geometric Mechanics of Solids: a new analysis of modern engineering materials"_ (GEMS).


<h2> License </h2>

Distributed under the GNU General Public License v3.0. See LICENSE.txt for more information.


<h2> Contacts </h2>
Please feel free <a href = "mailto: Elijah.Borodin@icloud.com"> e-mail </a> to Dr Elijah Borodin (Research Fellow in Materials Physics at the University of Manchester, Department of Solids and Structures) any queries relating with the code.
