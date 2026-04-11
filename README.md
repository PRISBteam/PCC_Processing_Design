<h1>Polytopal Cell Complex (PCC) Processing Design  (CPD code)</h1>

Manual version: 0.5.0 <br>
Current manual release date: 3 November 2025 <br>

<p> <i> PCC Processing Design </i> or <i> CPD code </i> is a software designed to deliver an entirely discrete, combinatorial representation of complex multidimensional and multiscale defect microstructures, along with detailed characterisation of their evolution during various processing routes. Additionally, it can be utilized as a highly efficient design tool for microstructure optimisation in composite materials. 
The key methodological feature is the use of polytopal (polygons in the 2D case or polyhedra in 3D) cell complexes (PCCs) as objects of algebraic topology, which provide a convenient discrete space for large-scale design of realistic material defect microstructures of different dimensions, such as point and line defects, interfaces, and various volumetric phases. Such PCCs can be created based on the 2D/3D tessellation of space with the convex polytopes, such as Voronoi polyhedra.

An excellent, simple introduction to the cell complexes with their various applications is given in the <a href="https://link.springer.com/book/10.1007/978-1-84996-290-2" target="_blank"> book </a> of Leo Grady and Jonathan Polimeni <i>"Discrete Calculus. Applied Analysis on Graphs for Computational Science" (2010).</i>
More rigorous mathematical introduction to the main concepts and methods of combinatorial algebraic topology is given in the <a href="https://doi.org/10.1007/978-3-540-71962-5" target="_blank"> book </a> of Dmitry Kozlov <i> "Combinatorial Algebraic Topology" (2008) </i>. A comprehensive classical introduction to graph theory is provided in the classical <a href="https://doi.org/10.1007/978-1-4612-0619-4" target="_blank"> book </a> of Béla Bollobás <i> "Modern Graph Theory" (1998)</i>.

For the effective use of the CPD code, a clear understanding of several basic concepts, such as Polytopal Cell Complexes, state vectors, and special sequences of cells, is necessary. All these concepts (and much more) are described in the <b>Theoretical Manual</b> included in the project repository. Even more comprehensive discussions of some particular topics with examples can be found in the following publications of the Mechanics of Solids research group from the University of Manchester (see also  <a href="https://materia.team/publications.html" target="_blank"> MATERiA publication list </a>):
<ol>
<li> E.N. Borodin, A.G. Sheinerman, O.Yu. Bushuev, M.Yu. Gutkin, A.P. Jivkov (2024) Defect-induced fracture topologies in Al2O3 ceramic-graphene nanocomposites. Materials & Design, 239, 112783. <a href="https://doi.org/10.1016/j.matdes.2024.112783" target="_blank"> 10.1016/j.matdes.2024.112783 </a> </li>
<li> Siying Zhu, Elijah Borodin, Andrey P. Jivkov (2023) Topological characteristics of grain boundary networks during severe plastic deformations of copper alloys. Acta Materialia, 259, 119290. <a href="https://doi.org/10.1016/j.actamat.2023.119290" target="_blank"> 10.1016/j.actamat.2023.119290 </a> </li>
<li> E.N. Borodin, A.P. Jivkov, A.G. Sheinerman, M.Yu. Gutkin (2021) Optimisation of rGO-enriched nanoceramics by combinatorial analysis. Materials & Design, 212, 110191. <a href="https://doi.org/10.1016/j.matdes.2021.110191" target="_blank"> 10.1016/j.matdes.2021.110191 </a> </li>
<li> S. Zhu, E.N. Borodin, A.P. Jivkov (2021) Triple junctions network as the key structure for characterisation of SPD processed copper alloys. Materials & Design, 198(24), 109352. <a href="https://doi.org/10.1016/j.matdes.2020.109352" target="_blank"> 10.1016/j.matdes.2020.109352 </a> </li>
<li> E. N. Borodin, A. P. Jivkov (2019) Evolution of triple junctions’ network during severe plastic deformation of copper alloys – a discrete stochastic modelling. Philosophical Magazine, 100(4), 467-485.<a href="https://doi.org/10.1080/14786435.2019.1695071" target="_blank"> 10.1080/14786435.2019.1695071 </a> </li>
</ol>

All the necessary technical details about the CPD code internal architecture and capabilities can be read in the <b>Technical Manual</b> included in the project repository.

<h2> How to execute the code </h2>
<p> This is a C++ based software project consisting of several modules, libraries, and the main.cpp file, which is considered a separate module. The code works with a pre-created PCC, represented as a set of its incidence and adjacency matrices in a sparse matrix form. It is intended to be launched as a CMake project using the CMakeLists.txt file. For successful code execution, both the C++ compiler and CMake software must be preinstalled. 

The code is written and tested in C++17. It works well with <a href="http://cmake.org" target="_blank"> CMake 3.23 </a>, <a href="http://gcc.gnu.org" target="_blank"> g++ compiler </a> and <a href="http://jetbrains.com/clion" target="_blank"> CLion IDE </a>.
<-- It is partly parallelized with the <a href="https://www.openmp.org" target="_blank">OpenMP</a> libraries for its effective execution using multiple cores simultaneously within a single CPU. -->
It explicitly uses the <a href="http://spectralib.org/download.html">Eigen and Spectra</a> C++ libraries included in the project as its external libraries.

All the “user interface” allowing interaction with the code contains in a few 'config/*.ini' files governing the behaviour of each of the modules. For writing tailored C++ scripts that implement project modules and libraries in loops or in an arbitrary order, the 'TASK' execution mode can be used, with user code written in 'tasks/*.cpp' files. Please refer to the <b>Technical Manual</b> included in the project repository for more detailed guidance.

As an example, to compile the project in a command line tool, one needs to change the working directory (cd command in Linux) to the one containing the project CMakeLists.txt file and then launch CMake as:
```
cmake -B buildtree
cmake --build buildtree
```
After that, one needs to change the working directory again to the newly created \buildtree and execute the built file like
```
cd buildtree
 ./PCC_Processing_Design
```
With CLion and other IDEs, everything is even more straightforward: a new C++ project must be created (if it contains its own main.cpp file with the “Hello world!” code by default, it must be deleted or ignored) and then executed.

The computational costs of different calculation types, functions, and tasks vary significantly: for instance, spectral analysis of matrices or genetic algorithms are very time-consuming procedures, while the random walk algorithm and generation of chains of special cells are fast. <br>

The code works equally well with 3D and 2D tessellations of space. In the 2D case, there are no 3D polyhedra (volumes), and 2D polytopes are associated with faces or 2-cells of the corresponding PCC. All the project functions work similarly in these two cases.
</p>

<h2> Where to take a PCC? </h2>

A large PCC library is published on the <a href="https://materia.team/Voronoi_library.html" target="_blank"> MATERiA</a> project web page.

<h3> Tessellations of space provided by Neper software </h3>

A variety of polytopal tessellations of 2D and 3D domains can be generated by <a href="https://neper.info/" target="_blank">Neper</a> software. It
allows output with with several different  the <a href="https://neper.info/doc/neper_t.html#morphology-options" target="_blank"> morphologies </a> and tesselation options <i> -morpho <morphology> </i>, such as <i> cube, square, tocta, lamellar, etc. </i> different from <i>Voronoi</i>. Please, see more <a href="https://neper.info/doc/neper_t.html#examples" target="_blank"> examples </a> on the Neper webpage.


<h3> PCC Generator Tool </h3>

The code Voronoi PCC Analyser in the same MATERiA software project GitHub organisation allows to create a PCC as a set on incidence and adjacency matrices directly from the *.tess file as an output of the Neper tesselation module. See all the other related information on the <a href="https://github.com/PRISBteam/Voronoi_PCC_Analyser" target=”_blank”> GitHub repository </a>. Its latest release can be downloaded from the <a href="https://github.com/PRISBteam/Voronoi_PCC_Analyser/releases" target=”_blank”> page of releases </a> or installed directly by pip (by 'pip install PCCanalyser==0.1.0' terminal command) or from the <a href="https://pypi.org/project/PCCanalyser/0.1.0/" target=”_blank”> PCCanalyser </a> web page.


<h2> Applications of the CPD code </h2>

Please refer to the <a href="https://materia.team" target="_blank"> MATERiA</a> project web page, and, in particular, to its <a href="https://materia.team/simulations.html" target="_blank">  Simulation Examples </a> section for more practical examples and the applications of the CPD code.


<h2> Acknowledgements </h2>

This code has been created as a part of the EPSRC-funded projects EP/V022687/1 _“Patterns recognition inside shear bands: tailoring microstructure against localisation”_ (PRISB) and EP/N026136/1 _"Geometric Mechanics of Solids: a new analysis of modern engineering materials"_ (GEMS).


<h2> License </h2>

Distributed under the standard MIT License. Please see the 'LICENSE.txt' file in the code repository for more information.


<h2> Contacts </h2>

Feel free to <a href = "mailto: Elijah.Borodin@icloud.com"> e-mail </a> to Dr Elijah Borodin (Lecturer in Mechanics of Solids at the University of Manchester, School of Engineering) with any code-related queries.