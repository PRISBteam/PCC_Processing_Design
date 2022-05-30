///================================ DCC Processing module =============================================================///
///=======================================================================================================================///
/** The function in this library generate quasi-random process of changes in the elements of the pre-constructed          ///
*   discrete sell complex (DCC) with the whole set of incidence and adjacency martices                                   **/
///=======================================================================================================================///

/// Standard (STL) C++ libraries:
///------------------------------
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
///------------------------------

/// Attached user defined C++ libraries:
///-------------------------------------
#include "Eigen/Core"
#include "Spectra/GenEigsSolver.h"
#include "Spectra/SymEigsSolver.h"
#include "Spectra/MatOp/SparseGenMatProd.h"

#include "TJsLab.h"
///-------------------------------------

using namespace std;
using namespace Eigen;
using namespace Spectra;

void HAGBsProbability3D(char* AN, char* AE, char* AF, char* AG, char* MEN, char* MFE, char* MGF, char* CNumbs, char* config, char* odir, char stype) {
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// odir - source path
// stype - simulation type ('R', 'S', 'I',....)

/// Declaration of FUNCTIONS, see the function bodies at the end of file.
    void confCout(vector<bool> configuration, int number_of_types); // Configuration's of the problem output read from file "config"
    Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols); // The function read any matrices from lists of triplets and create corresponding sparse matrices
    std::vector<unsigned int> VectorReader(char* FilePath); // The function read any matrices from lists of triplets and create corresponding sparse matrices
    std::vector<int> confCout(char* config, vector<int> const& configuration); // Read and Output configuration

////////////////////// Matrices initialisation part //////////////////////
// Triplets in the form T = T(i,j,value), where i and j element's indices in the corresponding dense matrix
    typedef Triplet<double> Tr; // <Eigen library class> Declares a triplet's type name - Tr
    typedef SparseMatrix<double> SpMat; // <Eigen library class> Declares a column-major sparse matrix type of doubles name - SpMat
    std::vector<Tr> SFaces_Triplet_list; // Probe vector of triplets

    long unsigned int i = 0, t_length = 0;
    std::vector<unsigned int> CellNumbs;

////////==================== Reading from files of all the adjacency matrices of the DCC ===================////

/// Reading vector from the file "number_of_cells" the numbers of cells od different types ///
/// ::      vector components: [0] - Nodes, [1] - Edges, [2] - Faces, [3] - Grains ::     ///
    CellNumbs = VectorReader(CNumbs);

    //Screen output for the numbers of cells in the DCC
    cout << "The number of different k-cells in the DCC:" << endl;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;

/// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
    SpMat ANS(CellNumbs.at(0)+1,CellNumbs.at(0)+1), AES(CellNumbs.at(1)+1,CellNumbs.at(1)+1),
            AFS(CellNumbs.at(2)+1,CellNumbs.at(2)+1), AGS(CellNumbs.at(3)+1,CellNumbs.at(3)+1);
    ANS = SMatrixReader(AN, (CellNumbs.at(0)+1), (CellNumbs.at(0)+1)); // Nodes
    AES = SMatrixReader(AE, (CellNumbs.at(1)+1), (CellNumbs.at(1)+1)); //Edges
    AFS = SMatrixReader(AF, (CellNumbs.at(2)+1), (CellNumbs.at(2)+1)); //Faces
    AGS = SMatrixReader(AG, (CellNumbs.at(3)+1), (CellNumbs.at(3)+1)); // Volumes
/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
    SpMat ENS(CellNumbs.at(1)+1,CellNumbs.at(0)+1), FES(CellNumbs.at(2)+1,CellNumbs.at(1)+1), GFS(CellNumbs.at(3)+1,CellNumbs.at(2)+1);
    ENS = SMatrixReader(MEN, (CellNumbs.at(1)+1), (CellNumbs.at(0)+1)); // Nodes-Edges
    FES = SMatrixReader(MFE, (CellNumbs.at(2)+1), (CellNumbs.at(1)+1)); // Edges-Faces
    GFS = SMatrixReader(MGF, (CellNumbs.at(3)+1), (CellNumbs.at(2)+1)); // Faces-Grains
    cout << "The number of different k-cells in the DCC:" << endl;

    /// Generation of the output streams
    ofstream OutFaceFile; // Output stream for variables related with special 2-Cells in DCC
    ofstream OutFLfile; // Special Face Laplacian
    ofstream OutRWFLfile; // Random Walker Special Face Laplacian

    OutFaceFile.open(odir + "FaceAdjacency.txt"s, ios::trunc);
    OutFaceFile<< "Adjacency Matrix of all Special Faces" << endl;
    OutFaceFile.close();
    OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
    OutFLfile<< "Laplacian Matrix of all Special Faces" << endl;
    OutFLfile.close();
    OutRWFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
    OutRWFLfile<< "Random Walker Laplacian Matrix of all Special Faces" << endl;
    OutRWFLfile.close();

    ofstream OutTJsFile; // Output stream for variables related with special 1-Cells in DCC
    OutTJsFile.open(odir + "TJsLab_TJsTypes.txt"s, ios::trunc);
    OutTJsFile << "Junctions [1-Cells] of different types" << endl;
    OutTJsFile << "Strain\t\t" <<"J0:\t" << "\t J1:\t" << "\t J2:\t" << "\t J3:\t" << endl;

    ofstream OutSFile;// Output stream for configurational entropy related with 1-Cells
    OutSFile.open(odir + "TJsLab_ConTJsEntropy.txt"s, ios::trunc);
    OutSFile << "Strain\t" <<"\t Configurational [1-Cells] Entropy \t" << endl;

////////////////////////////////////// KINETIC PROCESS ///////////////////////////////////////////////
    if (stype == 'R') { //  Random generation case

/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen ///
        std:vector<int> configuration = confCout(config, configuration);
        int number_of_types = configuration.at(0);

///=============================================================================================================================================////
/// ====================================================>  Random generation process   <========================================================////
///=============================================================================================================================================////
        std::vector<bool> SpecialCells(CellNumbs.at(2), 0), OrdinaryCells(CellNumbs.at(2),1); // New vectors initialised with 0 for SpecialCells and 1 for OrdinaryCells

        map<unsigned int,unsigned int> SpecialCellMap; // Mapping [k]->[l] from the set of all 2-Cells in the initial DCC to the set of newly generated special cells
        map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
        double ordinary_faces_fraction = 1.0, special_faces_fraction = 0.0;
        std::vector<unsigned int> SpecialCellNumbs(CellNumbs.at(2), 0); // Vector of the size equal to the total number of faces in DCC initialised with '0's
        for (unsigned int lit = 0; lit < SpecialCellNumbs.size(); lit++) SpecialCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

        /// Numerates newly created Faces during the random generation process
        long unsigned int numerator = 0;
        /// Maximal fraction for simulation loop max_sFaces_fraction = [0,1]
        double max_sFaces_fraction = 0.8;
        /// Step for output (structural analysis and data output will be performed after each output_step calculation steps (= number of newly converted elements))
        int output_step = 700;
///=============================================================================================================================================////
/// ================= Loop over all possible fractions of special cells [0,1] =======================>

/// The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
        srand (time(NULL));

        do { // do{ ... }while(output_step) loop starting point
            int New2CellNumb = 0, NewFaceType = 1; // Only one possible Face type (binary model)

            New2CellNumb = rand() % (SpecialCellNumbs.size()-1); // Random generation of the boundary number in the range from 0 to SpecialCellNumbs.size()-1
            if (number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types
            OrdinaryCells.at(SpecialCellNumbs.at(New2CellNumb)) = 0; // Replace the chosen element with 0 instead of 1 in the Original Faces vector

            SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)] = numerator++; // Assign new elements to the map and increase numerator
            long unsigned int sit = SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)]; // Just a useful variable for the number of newly converted Face (= numerator--)

            /// It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
            SpecialCellNumbs.erase(SpecialCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size BUT

            /// Another key point: creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
            for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
                if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
                    SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

            // Special and Ordinary Faces fraction calculation
            unsigned int OCellAmount = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1);
            ordinary_faces_fraction = OCellAmount / (double) CellNumbs.at(2);
            special_faces_fraction = 1.0 - ordinary_faces_fraction;

            ///=============================================================================================================================================////
            ///=================================== Characterisation module =============================///

            if( OCellAmount % output_step == 0 && ordinary_faces_fraction > 0.05) { // Periodicity of characterisation output

                /// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces /// Sparse adjacency matrix from the corresponding triplet list of special faces
                SpMat Id(numerator, numerator), SAM_FacesGraph(numerator, numerator), Face_Laplacian(numerator, numerator), Sym_Face_Laplacian(numerator, numerator), RW_Face_Laplacian(numerator, numerator);
                SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

                /// Degree vector and matrix
                vector<unsigned int> SFace_degrees;
                for (unsigned int i = 0; i < numerator; i++) { //columns
                    unsigned int degree_Fcounter = 0;
                    for (unsigned int j = 0; j < numerator; j++) { // rows
                        if (SAM_FacesGraph.coeff(j, i) == 1 && i != j)
                            degree_Fcounter++;
                    }
                    SFace_degrees.push_back(degree_Fcounter);
                }
                //Creation of the S-Face degree matrix
                SpMat SFDegree(SFace_degrees.size(), SFace_degrees.size());
                SFDegree.setIdentity();
                unsigned int numerator_degree = 0;
                for (double num: SFace_degrees) {
                    SFDegree.coeffRef(numerator_degree,numerator_degree) = SFace_degrees.at(numerator_degree);
                    numerator_degree++;
                }
                /// Identity matrix
                Id.setIdentity();
                /// Special Faces Laplacian matrix
                Face_Laplacian = SFDegree - SAM_FacesGraph;
                // Symmetric S-Faces Laplacian [D^-1/2 L ]
                /// Left Random-Walk normalized Laplacian matrix
                RW_Face_Laplacian = Id - SFDegree.cwiseInverse() * SAM_FacesGraph;

                /// Output to file Special Faces Laplacian matrix
                OutFaceFile.open(odir + "FaceAdjacency.txt"s, ios::trunc);
                OutFaceFile<< "Adjacency Matrix of all Special Faces" << endl;
                OutFaceFile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;                 // Output of special_faces_fraction to files
                OutFaceFile << SAM_FacesGraph << endl;
                //  Dense matrix output:      OutFaceFile << MatrixXd(SAM_FacesGraph) << endl;
                OutFaceFile.close();

                OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
                OutFLfile<< "Laplacian Matrix of all Special Faces" << endl;
                OutFLfile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;
                OutFLfile << Face_Laplacian << endl;
                // Dense matrix output:        OutFLfile << MatrixXd(Face_Laplacian) << endl;
                OutFLfile.close();

                OutRWFLfile.open(odir + "RWFaceLaplacian"s, ios::trunc);
                OutRWFLfile<< "Random Walker Laplacian Matrix of all Special Faces" << endl;
                OutRWFLfile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;
                OutRWFLfile <<RW_Face_Laplacian << endl;
                // Dense matrix output:         OutRWFLfile << MatrixXd(RW_Face_Laplacian) << endl;
                OutRWFLfile.close();

                /// Creation of the sparse incidence matrix (FES_Graph) for special Faces sparse FES matrices
                //  SpMat FES_Graph(CellNumbs.at(2) + 1,CellNumbs.at(1) + 1);
                // SpMat FES_Graph;
                // FES_Graph*FES_Graph.transpose()

                if (configuration.at(1)) { // Nodes types statistics, indices and configuration entropy
                    //#include "QPsLab.h"
                }
                if (configuration.at(2)) { /// Edges types statistics, indices and configuration entropy
                    EdgesStat(CellNumbs, numerator, SpecialCellMap, FES, odir, special_faces_fraction);
                }
                if (configuration.at(3)) { // Faces types statistics and structural indices
                    //#include "FacesLab.h"
                }
                if (configuration.at(4)) { // Grains types statistics, indices and configuration entropy
                    //#include "GrainsLab.h"
                }
                if (configuration.at(5)) { // Nodes Laplacian with its spectrum for the nodes network
                    //#include"NodeLaplacians.h"
                }
                if (configuration.at(6)) { // Edges Laplacian with its spectrum for the edges network
                    //#include"EdgeLaplacians.h"
                }
                if (configuration.at(7)) { // Face Laplacian with its spectrum for the special faces network
                }
                if (configuration.at(8)) { // Betti numbers (FB0, FB1 and FB2) of the special faces network
                    //#include"BettiNumbers.h"
                }
                if (configuration.at(9)) { // Grain Laplacian with its spectrum for the grain network
                    //#include"GrainLaplacians.h"
                }
                if (configuration.at(10)) { // Tutte polynomial for the special faces network
                    //#include"FaceTutte.h"
                }
                if (configuration.at(11)) { // Statistical physics module
                    //#include"DCCStatistical.h"
                }

                cout << "Fraction of special faces:" << 1.0 - ordinary_faces_fraction << ";\t Size of S-Face Adjacency matrix \t" << SAM_FacesGraph.nonZeros() << endl;
            } // End of analysis and output iterator ( IF: iterator % X == 0 )
        }while(ordinary_faces_fraction > (1.0 - max_sFaces_fraction)); /// End of the Random generation process

/// Closing and deleting
        OutFaceFile.close(); // 2-Cells related matrices
        OutTJsFile.close(); // 1-Cells statistics output
        OutSFile.close(); // 1-Cells Configurational Entropy output
        //Remove all elements anf free the memory from the probe SpecialCells and vector
        SpecialCells.clear();
        OrdinaryCells.clear();
        SpecialCells.shrink_to_fit();
        OrdinaryCells.shrink_to_fit();

    } ///End of 'R' type simulations

///=============================================================================================================================================////
///=============================================================================================================================================////
/// ====================================================>  Maximum entropy production process   <========================================================////
///=============================================================================================================================================////
    else if (stype == 'S') { // Maximum entropy production
/// Read simulation configuration from file :: the number of special face types and calculating parameters. Then Output of the current configuration to the screen ///
        vector<int> configuration = confCout(config, configuration);
        int number_of_types = configuration.at(0);

        std::vector<bool> SpecialCells(CellNumbs.at(2), 0), OrdinaryCells(CellNumbs.at(2),1); // New vectors initialised with 0 for SpecialCells and 1 for OrdinaryCells
        map<unsigned int,unsigned int> SpecialCellMap; // Mapping [k]->[l] from the set of all 2-Cells in the initial DCC to the set of newly generated special cells
        map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
        double ordinary_faces_fraction = 1.0, special_faces_fraction = 0.0;
        std::vector<unsigned int> SpecialCellNumbs(CellNumbs.at(2), 0); // Vector of the size equal to the total number of faces in DCC initialised with '0's
        for (unsigned int lit = 0; lit < SpecialCellNumbs.size(); lit++) SpecialCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

        /// Numerates newly created Faces during the random generation process
        long unsigned int numerator = 0;

/// The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
        srand (time(NULL));

/// At first, random generation of small amounts about 0.05 of special seeds

/////////////////////////////////  Maximum Entropy Production ///////////////////////////////////
        do { // do{ ... }while(output_step) loop for seed points
            int New2CellNumb = 0, NewFaceType = 1; // Only one possible Face type (binary model)

            New2CellNumb = rand() % (SpecialCellNumbs.size()-1); // Random generation of the boundary number in the range from 0 to SpecialCellNumbs.size()-1
            if (number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types
            OrdinaryCells.at(SpecialCellNumbs.at(New2CellNumb)) = 0; // Replace the chosen element with 0 instead of 1 in the Original Faces vector

            SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)] = numerator++; // Assign new elements to the map and increase numerator
            long unsigned int sit = SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)]; // Just a useful variable for the number of newly converted Face (= numerator--)

            /// It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
            SpecialCellNumbs.erase(SpecialCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size BUT

            // Special and Ordinary Faces fraction calculation
            unsigned int OCellAmount = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1);
            ordinary_faces_fraction = OCellAmount / (double) CellNumbs.at(2);
            special_faces_fraction = 1.0 - ordinary_faces_fraction;

        }while(ordinary_faces_fraction > 0.95); /// End of the Random generation process
/////////////////////////////////////////////////////////////////////////////////////////////
        /// Maximal fraction for simulation loop max_sFaces_fraction = [0,1]
        double max_sFaces_fraction = 0.3;
        /// Step for output (structural analysis and data output will be performed after each output_step calculation steps (= number of newly converted elements))
        int output_step = 70;
///=============================================================================================================================================////
/// ================= Loop over all possible fractions of special cells [0,1] =======================>
        /// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces /// Sparse adjacency matrix from the corresponding triplet list of special faces
        SpMat Id(numerator, numerator), SAM_FacesGraph(numerator, numerator), Face_Laplacian(numerator, numerator), Sym_Face_Laplacian(numerator, numerator), RW_Face_Laplacian(numerator, numerator);
        SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

        do { // do{ ... }while(output_step) loop starting point
            int New2CellNumb = 0, NewFaceType = 1; // Only one possible Face type (binary model)

            vector<int> TJsTypes(CellNumbs.at(1)+1,0);
            vector<double> EntropyIncreaseList(CellNumbs.at(2));

            double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
            double Configurational_Face_Entropy = 0;
            for(int i = 1; i < numerator; i++)
                for(int j = 1; j < numerator; j++)
                    if( i != j && SAM_FacesGraph.coeff(i,j) == 1) for(int k = 1; k < CellNumbs.at(1); k++) if( FES.coeff(j,k) == 1) TJsTypes.at(k)++;

            J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
            J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
            J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
            J0 = CellNumbs.at(1) - J1 - J2 - J3;
            Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
            j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;

            /// Configurational Entropy related with Faces
            Configurational_Face_Entropy = - (j0* log2(j0) + j1* log2(j1) + j2* log2(j2) + j3* log2(j3));

            double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
            for (unsigned int k = 0; k < CellNumbs.at(2); k++) { //loop over all the Faces in DCC

                // Loop over each element neighbours
                TJsTypes.clear(); J00 = 0; J10 = 0; J20 = 0; J30 = 0;  J0N = 0; J1N = 0; J2N = 0; J3N = 0;
                for(int j = 1; j < numerator; j++) // For each element loop over all his neighbours
                        if( i != j && SAM_FacesGraph.coeff(i,j) == 1) for(int k = 1; k < CellNumbs.at(1); k++) if( FES.coeff(j,k) == 1) TJsTypes.at(k)++;
                J00 = std::count(TJsTypes.begin(), TJsTypes.end(), 0);
                J10 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
                J20 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
                J30 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
                J0N = J00 + 1; J1N =  J10 + 1; J2N = J20 + 1; J3N = J30 + 1;
                // The entropy increase calculation for a given Face
                CFace_EntropyIncrease = (J0-J00+J0N) * log2(J0-J00+J0N) + (J1-J10+J1N) * log2(J1-J10+J1N) + (J2-J20+J2N) * log2(J2-J20+J2N) + (J3-J30+J3N) * log2(J3-J30+J3N);
                EntropyIncreaseList.push_back(CFace_EntropyIncrease);
            }
            // Number of element giving the maximum increase in configurational entropy at its conversion
            double CFace_max =
            New2CellNumb = (int) *max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)); //'*' because this fuction return a pointer to the element, not the element itself

            // Then all the corresponding maps chan
            OrdinaryCells.at(SpecialCellNumbs.at(New2CellNumb)) = 0; // Replace the chosen element with 0 instead of 1 in the Original Faces vector
            SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)] = numerator++; // Assign new elements to the map and increase numerator
            long unsigned int sit = SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)]; // Just a useful variable for the number of newly converted Face (= numerator--)

            /// Creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
            for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
                if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
                    SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

            // Special and Ordinary Faces fraction calculation
            unsigned int OCellAmount = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1);
            ordinary_faces_fraction = OCellAmount / (double) CellNumbs.at(2);
            special_faces_fraction = 1.0 - ordinary_faces_fraction;

            ///=============================================================================================================================================////
            ///=================================== Characterisation module =============================///
            if( OCellAmount % output_step == 0 && ordinary_faces_fraction > 0.05) { // Periodicity of characterisation output

                /// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces /// Sparse adjacency matrix from the corresponding triplet list of special faces
                SpMat Id(numerator, numerator), SAM_FacesGraph(numerator, numerator), Face_Laplacian(numerator, numerator), Sym_Face_Laplacian(numerator, numerator), RW_Face_Laplacian(numerator, numerator);
                SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

                /// Degree vector and matrix
                vector<unsigned int> SFace_degrees;
                for (unsigned int i = 0; i < numerator; i++) { //columns
                    unsigned int degree_Fcounter = 0;
                    for (unsigned int j = 0; j < numerator; j++) { // rows
                        if (SAM_FacesGraph.coeff(j, i) == 1 && i != j)
                            degree_Fcounter++;
                    }
                    SFace_degrees.push_back(degree_Fcounter);
                }
                //Creation of the S-Face degree matrix
                SpMat SFDegree(SFace_degrees.size(), SFace_degrees.size());
                SFDegree.setIdentity();
                unsigned int numerator_degree = 0;
                for (double num: SFace_degrees) {
                    SFDegree.coeffRef(numerator_degree,numerator_degree) = SFace_degrees.at(numerator_degree);
                    numerator_degree++;
                }
                /// Identity matrix
                Id.setIdentity();
                /// Special Faces Laplacian matrix
                Face_Laplacian = SFDegree - SAM_FacesGraph;
                // Symmetric S-Faces Laplacian [D^-1/2 L ]
                /// Left Random-Walk normalized Laplacian matrix
                RW_Face_Laplacian = Id - SFDegree.cwiseInverse() * SAM_FacesGraph;

                /// Output to file Special Faces Laplacian matrix
                OutFaceFile.open(odir + "FaceAdjacency.txt"s, ios::trunc);
                OutFaceFile<< "Adjacency Matrix of all Special Faces" << endl;
                OutFaceFile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;                 // Output of special_faces_fraction to files
                OutFaceFile << SAM_FacesGraph << endl;
                //  Dense matrix output:      OutFaceFile << MatrixXd(SAM_FacesGraph) << endl;
                OutFaceFile.close();

                OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
                OutFLfile<< "Laplacian Matrix of all Special Faces" << endl;
                OutFLfile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;
                OutFLfile << Face_Laplacian << endl;
                // Dense matrix output:        OutFLfile << MatrixXd(Face_Laplacian) << endl;
                OutFLfile.close();

                OutRWFLfile.open(odir + "RWFaceLaplacian"s, ios::trunc);
                OutRWFLfile<< "Random Walker Laplacian Matrix of all Special Faces" << endl;
                OutRWFLfile << endl << "Special Faces fraction is equal to\t " << special_faces_fraction  << endl << endl;
                OutRWFLfile <<RW_Face_Laplacian << endl;
                // Dense matrix output:         OutRWFLfile << MatrixXd(RW_Face_Laplacian) << endl;
                OutRWFLfile.close();

                /// Creation of the sparse incidence matrix (FES_Graph) for special Faces sparse FES matrices
                //  SpMat FES_Graph(CellNumbs.at(2) + 1,CellNumbs.at(1) + 1);
                // SpMat FES_Graph;
                // FES_Graph*FES_Graph.transpose()

                if (configuration.at(1)) { // Nodes types statistics, indices and configuration entropy
                    //#include "QPsLab.h"
                }
                if (configuration.at(2)) { /// Edges types statistics, indices and configuration entropy
                    EdgesStat(CellNumbs, numerator, SpecialCellMap, FES, odir, special_faces_fraction);
                }
                if (configuration.at(3)) { // Faces types statistics and structural indices
                    //#include "FacesLab.h"
                }
                if (configuration.at(4)) { // Grains types statistics, indices and configuration entropy
                    //#include "GrainsLab.h"
                }
                if (configuration.at(5)) { // Nodes Laplacian with its spectrum for the nodes network
                    //#include"NodeLaplacians.h"
                }
                if (configuration.at(6)) { // Edges Laplacian with its spectrum for the edges network
                    //#include"EdgeLaplacians.h"
                }
                if (configuration.at(7)) { // Face Laplacian with its spectrum for the special faces network
                }
                if (configuration.at(8)) { // Betti numbers (FB0, FB1 and FB2) of the special faces network
                    //#include"BettiNumbers.h"
                }
                if (configuration.at(9)) { // Grain Laplacian with its spectrum for the grain network
                    //#include"GrainLaplacians.h"
                }
                if (configuration.at(10)) { // Tutte polynomial for the special faces network
                    //#include"FaceTutte.h"
                }
                if (configuration.at(11)) { // Statistical physics module
                    //#include"DCCStatistical.h"
                }

                cout << "Fraction of special faces:" << 1.0 - ordinary_faces_fraction << ";\t Size of S-Face Adjacency matrix \t" << SAM_FacesGraph.nonZeros() << endl;
            } // End of analysis and output iterator ( IF: iterator % X == 0 )
        }while(ordinary_faces_fraction > (1.0 - max_sFaces_fraction)); /// End of the Random generation process

/// Closing and deleting
        OutFaceFile.close(); // 2-Cells related matrices
        OutTJsFile.close(); // 1-Cells statistics output
        OutSFile.close(); // 1-Cells Configurational Entropy output
        //Remove all elements anf free the memory from the probe SpecialCells and vector
        SpecialCells.clear();
        OrdinaryCells.clear();
        SpecialCells.shrink_to_fit();
        OrdinaryCells.shrink_to_fit();

    } ///End of 'S' type simulations

    else if (stype == 'I') {

    }

    else if (stype == 'E') {

    }
    else { cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl; return;}


    return;
} /// The end of HAGBsProbability3D()

/// 2D version (DCC with max 2-Cells without 3-Cells and higher dimensions) of the code
void HAGBsProbability2D(char* AN, char* AE, char* AF, char* MEN, char* MFE, char* CNumbs, char* config, char* odir, char stype) {
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// odir - source path
// stype - simulation type ('R', 'S', 'I',....)
} /// The end of HAGBsProbability2D()


/// ================================== Related functions ==================================///
std::vector<int> confCout(char* config, vector<int> const& configuration) {
    std::string line;
    std:vector<int> res;
    ifstream inConf(config);
    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n')) {
            if(line.at(0) == '!') { res.push_back(line.at(2)-'0');
            }// Number of types from the first line of the file
            if(line.at(0) == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
            if(line.at(0) == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
        }
    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "The number of special face types:\t" << configuration.at(0) << endl;
    cout << "Nodes types statistics, indices and configuration entropy:     "; if (configuration.at(1) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges types statistics, indices and configuration entropy:     "; if (configuration.at(2) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Faces types statistics and structural indices:                 "; if (configuration.at(3) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain types statistics, indices and configuration entropy:     "; if (configuration.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Nodes Laplacian with its spectrum for the nodes network:       "; if (configuration.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges Laplacian with its spectrum for the edges network:       "; if (configuration.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Face Laplacian with its spectrum for the special faces network:"; if (configuration.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Betti numbers (FB0, FB1 and FB2) of the special faces network: "; if (configuration.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain Laplacian with its spectrum for the grain network:       "; if (configuration.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Tutte polynomial for the special faces network:                "; if (configuration.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Statistical physics module:                                    "; if (configuration.at(11) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;

    return res;
}

Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols) {

    Eigen::SparseMatrix<double> res(Rows,Cols);
    typedef Triplet<double> Tr; // Eigen library class
    std::vector<Tr> tripletList; // Probe vector of triplets

    int i = 0, j = 0, value = 0, t_length = 0;
    ifstream inAN(SMpath);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << SMpath << " cannot be read" << endl; //If something goes wrong
//Sparse AN matrix
    res.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements anf free the memory from the probe vector
    tripletList.clear();
    tripletList.shrink_to_fit();

    return res;
}

std::vector<unsigned int> VectorReader(char* FilePath) {
    std::vector<unsigned int> res;
    unsigned int i=0;
    ifstream inCellNumbers(FilePath);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) res.push_back(i);
    } else std:cout << "The file " << FilePath << " cannot be read" << endl; //If something goes wrong

    return res;
}

/// Heap
// (1)    cout << std::numeric_limits<long unsigned int>::max() << endl;

// (2) ofstream NONStructGlCoordOut; NONStructGlCoordOut.open("C:VAnRes\\(MFE21)GNC_NONStructured.txt", ios::trunc);

//(3) ARCHIVE testing output  //   cout << "The number of columns in the created matrices:\t" << "\tANS --\t" << ANS.cols() << "\tAES --\t" << AES.cols()<< "\tAFS --\t" << AFS.cols() << "\tAGS --\t" << AGS.cols() << "\tENS --\t" << ENS.cols() << "\tFES --\t" << FES.cols() << "\tGFS --\t" << GFS.cols() << endl;
