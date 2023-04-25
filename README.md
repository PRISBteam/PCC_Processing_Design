<h1>Polyhedral Cell Complex (PCC) Processing Design  (DPD code)</h1>

Manual version: 0.1.0 <br>
Current manual release date: 25 April 2023 <br>

<p>  <i> PCC Processing Design  (DPD code)</i>  or   <i> Discrete Processing Design  (DPD code)</i> is a software aimed to provide a wholly discrete representation of defect microstructure evolution during the materials processing. The code is an effective tool both for studies microstructure network evolution during the materials processing and for design of advanced microstructures in composite polycrystalline materials. Its key feature is the usage of polyhedral cell complexes (PCCs) as objects of algebraic topology [1], which provide a convenient discrete space for large-scale design of realistic material defect structures of different dimensions (line defects, surface defects, different volumetric phases). Such PCCs can be constructed based on the 2D/3D space tessellation with the convex polyhedron’s such as the Voronoi tessellation of space realistically imitating a polycrystalline material.

The two following introductory sections discuss the basics of the PCCs theory needed for understanding the code input and output and the practical ways to obtain different PCCs as the sets of sparse matrices of their operators. </p>

<ol>
<li> K. Berbatov, P.D. Boom, A.L. Hazel, A.P. Jivkov, Applied Mathematical Modelling 110 (2022) 172-192. </li>
</ol>

<h2>Some necessary definitions</h2>

<p> All the k-cells with a specific value of k form a PCC <i>skeleton</i>, so that the whole PCC can be thought as an agglomeration of these several skeletons. </br>
Additionally to the scalar and vector parameters defined on the PCC cells, for each cell can be assigned a label characterising its own type. Hence all k-cells can be divided in two classes of <i>ordinary</i> (without labels or with the labels 0) and <i>special</i> (labelled) cells. </br>
Inside the class of special cells several types can be defined. The whole set of these labels assigned on the each skeleton form a configuration <i>state</i> which can be expressed by a single <i>state vector</i> encoding all the cell types like the well-known DNA sequences. </br>
The set of all the state-vectors in a PCC form its <i>structure</i>.  </br>
The notion of a <i>process</i> can be referred to as the discrete sequence of all the states expressed by their state vectors. 
Three different types of the structures can be distinguished: directly <i>assigned</i>, <i>imposed</i> by the cells of other dimensions and <i>induced</i> by a pre-existing assigned structure as a result of a kinetic process. </br>
The last component of the proposed approach is a function or <i>map</i> relating structural (such as fractions of special grain boundaries) and physical (time, plastic strain, mass fraction, etc. ) characteristics which can be  obtained based on the particular experimental data. </br>
Such a framework allows to describe and analyse space ordering inside a real materials by the well-developed tools of topology, statistics and graph theory.   </p>

<p> This is a C++ based software project consisting of several modules working with pre-created Polyhedral Cell Complex (PCC) as the set of its incidence and adjacency matrices represented in a sparse matrix form. The code contains 4 main modules (Processing, Characterisation and Writer), the library of Simulation tasks and the Main module merging everything together, reading configurations and data files,  and launching all the other modules.
The code works both with 3-complexes (3D) and 2-complexes (2D) BUT - to make results consistent with 2D/3D EBSD scans - it assumes that the grains are 3-cells in the 3D, 2-cells in the 2D case, and so on for grain boundaries and other element types. So it actually replaces definitions of k-cells with ( k + (dim - 3) )-cells, where dim = {2, 3} is the dimension of the problem. In such a way 2-cells are ALWAYS associated with grain boundaries on EBSD maps and are edges for 2D case! </p>
  
<h2>General specifications</h2>
  
<p>The code is written and tested in C++ 17 with the parallelised verson used some features of C++ 20. It is used explicitly Eigen and Spectra libraries which must be <a href="https://spectralib.org/download.html"> downloaded</a> and copied to the directory containing all the STL C++ libraries on the local PC.

The computational costs of different calculation types, functions and tasks are hugely different: for instance, the component analysis or Monte-Carlo simulations are a very time consuming procedures, while the random generation of special chains are fast.<\p>

<h2> More sources on the mathematics related to the Polyhedral cell complexes</h2>

An excellent simple introduction to the DCC with their various applications is given in the <a href="https://link.springer.com/book/10.1007/978-1-84996-290-2" target="_blank"> book </a> of Leo Grady and Jonathan Polimeni <i>“Discrete Calculus. Applied Analysis on Graphs for Computational Science. (2010)</i>

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
Please feel free <a href=“Elijah.Borodin@icloud.com” e-mail > to Dr Elijah Borodin (Research Fellow in Materials Physics at the University of Manchester, Department of Solids and Structures) any queries relating with the code.
