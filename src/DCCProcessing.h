///================================ DCC Processing module =============================================================///
///=======================================================================================================================///
/** The function in this library generate different quasi-random process for generation of the special elements on the    ///
*   elements of the pre-constructed discrete sell complex (DCC) with the whole set of incidence and adjacency martices    **/
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

/// Declaration of functions, see the function bodies at the end of file.
void confCout(vector<bool> configuration, int number_of_types); // Configuration's of the problem output read from file "config"
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, Eigen::SparseMatrix<double> const& SM1);

////////////////////// Matrices initialisation and reading from file //////////////////////
//Triplet T=T(i,j,value), where i and j element's indices in the corresponding matrix
typedef Triplet<double> Tr; // Eigen library class
typedef SparseMatrix<double> SpMat; // Declares a column-major sparse matrix type of double
std::vector<Tr> tripletList, SFaces_Triplet_list, ISFaces_Triplet_list; // Probe vector of triplets

//Reading from file by the in_streams creation:
/// Reading vector from the file "number_of_cells" the numbers of cells od different types :: vector components: [0] - Nodes, [1] - Edges, [2] - Faces, [3] - Grains
    long unsigned int i = 0, j = 0, value = 0, t_length = 0;
    std::vector<unsigned int> CellNumbs; //number_of_cells.out
    ifstream inCellNumbers(CNumbs);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) CellNumbs.push_back(i);
    } else cout << "The file " << CNumbs << " cannot be read" << endl; //If something goes wrong
    //Screen output for the numbers of cells in the DCC
    cout << "The number of different k-cells in the DCC:" << endl;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;

////////=============================== Reading from files of all the adjacency matrices of the DCC ===========================////
/// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
    SpMat ANS(CellNumbs.at(0) + 1,CellNumbs.at(0) + 1), AES(CellNumbs.at(1) + 1,CellNumbs.at(1) + 1),
            AFS(CellNumbs.at(2) + 1,CellNumbs.at(2) + 1), AGS(CellNumbs.at(3) + 1,CellNumbs.at(3) + 1);
    SMatrixReader(AN, ANS);
    SMatrixReader(AE, AES);
    SMatrixReader(AF, AFS);
    SMatrixReader(AG, AGS);
/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
    SpMat ENS(CellNumbs.at(1) + 1,CellNumbs.at(0) + 1), FES(CellNumbs.at(2) + 1,CellNumbs.at(1) + 1), GFS(CellNumbs.at(3) + 1,CellNumbs.at(2) + 1);
    SMatrixReader(MEN, ENS);
    SMatrixReader(MFE, FES);
    SMatrixReader(MGF, GFS);

////////////////////////////////////// KINETIC PROCESS ///////////////////////////////////////////////
    ofstream OutFaceFile; // Output stream for variables  related with special 2-Cells in DCC
    OutFaceFile.open("/Users/user/Dropbox/OFFICE/CProjects/Voro3D/test/SpecialFaceAdjacencyM.txt", ios::trunc);
    ofstream OutTJsFile; // Output stream for variables  related with special 2-Cells in DCC
    OutTJsFile.open("/Users/user/Dropbox/OFFICE/CProjects/Voro3D/test/TJsTypes.txt", ios::trunc);
    OutTJsFile << "J0:\t" << "\t J1:\t" << "\t J2:\t" << "\t J3:\t" << endl;
    OutTJsFile.close();
    ofstream OutSFile;// Output stream for configurational entropy related with 1-Cells
    OutSFile.open("/Users/user/Dropbox/OFFICE/CProjects/Voro3D/test/ConTJsEntropy.txt", ios::trunc);
    OutSFile << "\t Configurational [1-Cells] Entropy \t" << endl;
    OutSFile.close();
/////////////////////////////       Random generation case     ///////////////////////////
    if (stype == 'R') {

/// Read simulation configuration from file :: the number of special face types and calculating parameters
        std::vector<bool> configuration;
        std::string line;
        int number_of_types = 1;
        ifstream inConf(config);
        if (inConf) { //If the file was successfully open, then
            while(getline(inConf, line, '\n')) {
                if(line.at(0) == '!') { number_of_types = line.at(2)-'0';
                }// Number of types from the first line of the file
                if(line.at(0) == '#') configuration.push_back(1); // 1 and # means accept - the parameter will be calculated
                if(line.at(0) == '%') configuration.push_back(0); // 0 and % means ignore - the parameter will not be calculated
            }
        } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

/// ============    Output of the current configuration to the screen   ======================///
        confCout(configuration,number_of_types);

///=========================================================================================////
/// ==========================>  Random generation process   <==============================////
        std::vector<bool> SpecialCells(CellNumbs.at(2), 0), OrdinaryCells(CellNumbs.at(2),1);
        std::vector<unsigned int> SpecialCellNumbs(CellNumbs.at(2), 0);
        map<unsigned int,unsigned int> SpecialCellMap;
        for (unsigned int lit = 0; lit < SpecialCellNumbs.size(); lit++) SpecialCellNumbs[lit] = lit;

        double ordinary_edges_fraction = 1.0, special_edges_fraction = 0.0;
        std::vector<vector<bool>> VerySpecialFaces(number_of_types,SpecialCells); //Very Special Faces Vector -- vector of all special faces listed according to their global numeration in DCC
        std::vector<vector<unsigned int>> NewFacesList;

        long unsigned int NewFaceType = 0;
        map<unsigned int, unsigned int>::iterator sit;
        unsigned int itSC=0;

        srand (time(NULL)); // The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
        long unsigned int numerator = 0;
        /// Maximal fraction for simulation loop max_sFaces_fraction = [0,1]
        double max_sFaces_fraction = 0.3;
       /// Step for output (proportional to the number of iterations made)
        int output_step = 70;

/// ================= Loop over all possible fractions of special cells [0,1] =======================>
        do {
            long unsigned int New2CellNumb = 0, NewFaceType = 1;
            // Do while the loop did not find and original (still not converted) face

            New2CellNumb = rand() % (SpecialCellNumbs.size() - 1); // Random generation of the boundary number
            if (number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types

            SpecialCellMap[SpecialCellNumbs.at(New2CellNumb)] = numerator;
            long unsigned int sit = SpecialCellMap[SpecialCellNumbs[New2CellNumb]];
            for(unsigned int j = 0; j < CellNumbs.at(2); j++) if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1));
            numerator++;

            OrdinaryCells.at(SpecialCellNumbs.at(New2CellNumb)) = 0;
    //        cout << SpecialCellNumbs.at(New2CellNumb) << "\t\t"<<SpecialCellNumbs.size()<<endl;

            // Special Faces fraction calculation considering 1 - "ordinary cells" fraction
            ordinary_edges_fraction = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1) / (double) CellNumbs.at(2);
            //      cout << ordinary_edges_fraction << endl;

            ///=================================== Characterisation module =============================///
            unsigned int OCellAmount = std::count (OrdinaryCells.begin(), OrdinaryCells.end(), 1);
            if( OCellAmount % output_step == 0 && ordinary_edges_fraction > 0.05) { // Periodicity of characterisation output
                //cout << endl << numerator <<"\t"<< CellNumbs.at(2) - OCellAmount << "\t" << CellNumbs.at(2)  << endl; //exit(0);

                // cout << SFaces_Triplet_list.size() << endl;
                /// Creation of the sparse adjacency (SAM_FacesGraph) matrix for special Faces
                /// Sparse adjacency matrix from the corresponding triplet list of special faces
                SpMat SAM_FacesGraph(numerator, numerator);
                SAM_FacesGraph.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());

                /// Degree vector and matrix
                vector<int> SFace_degrees;
                for (unsigned int i = 0; i < numerator; i++) {
                    int degree_Fcounter = 0;
                    for (unsigned int j = 0; j < numerator; j++) {
                        SAM_FacesGraph.coeff(j, i);
                        if (j == 1 && i != j) degree_Fcounter++;
                    }
                    SFace_degrees.push_back(degree_Fcounter);
                    //   if (degree_Fcounter > 1) cout << "DC\t" << degree_Fcounter << endl;
                }

                SpMat SFDegree(SFace_degrees.size(), SFace_degrees.size());
                SFDegree.setIdentity();
                /// Creation of the S-Face degree matrix
                unsigned int numerator_degree = 0;
                for (double num: SFace_degrees) {
                SFDegree.coeffRef(numerator_degree,numerator_degree) = num;
                numerator_degree++;
            }

        ///Special Faces Adjacency matrix to file output
                OutFaceFile << endl << "Special Faces fraction is equal to\t " << (1.0 - ordinary_edges_fraction)  << endl << endl;

                OutFaceFile << MatrixXd(SAM_FacesGraph) << endl;
                OutFaceFile << SAM_FacesGraph  << endl;
//                SpMat FES(CellNumbs.at(2) + 1,CellNumbs.at(1) + 1);

                /// Creation of the sparse incidence matrix (FES_Graph) for special Faces sparse FES matrices
               // SpMat FES_Graph;
               // FES_Graph*FES_Graph.transpose()

                if (configuration.at(0)) { // Nodes types statistics, indices and configuration entropy
                    //#include "QPsLab.h"
                }
                if (configuration.at(1)) { /// Edges types statistics, indices and configuration entropy
                    EdgesStat(CellNumbs, numerator, SAM_FacesGraph, FES, odir);
                }
                if (configuration.at(2)) { // Faces types statistics and structural indices
                    //#include "FacesLab.h"
                }
                if (configuration.at(3)) { // Grains types statistics, indices and configuration entropy
                    //#include "GrainsLab.h"
                }
                if (configuration.at(4)) { // Nodes Laplacian with its spectrum for the nodes network
                    //#include"NodeLaplacians.h"
                }
                if (configuration.at(5)) { // Edges Laplacian with its spectrum for the edges network
                    //#include"EdgeLaplacians.h"
                }
                if (configuration.at(6)) { // Face Laplacian with its spectrum for the special faces network
                    //#include"FaceLaplacians.h"
                    //temp
                    SpMat Laplacian_SFaces(numerator,numerator);
                    //temp
                    SFDegree.setIdentity();
                    Laplacian_SFaces = SAM_FacesGraph + 2.0*SFDegree;
                    /// Eigenvalues of the S-Faces laplacian
                    //AFS(CellNumbs.at(2) + 1,CellNumbs.at(2) + 1);
 //                   EigenVals(AFS, Laplacian_SFaces); /// Important solver method
                    Eigen::VectorXcd eigval;
                    Eigen::MatrixXcd eigvec;
                    //SpMat eigvec;
                    //vector<vector<double>>
                    Eigen::SparseMatrix<double> Asym = 0.5*(SAM_FacesGraph+Eigen::SparseMatrix<double>(SAM_FacesGraph.transpose()));

                    // Construct matrix operation object using the wrapper class SparseGenMatProd
                    SparseGenMatProd<double> op(SAM_FacesGraph);
                    // Construct eigen solver object, requesting the largest three eigenvalues
                    GenEigsSolver<SparseGenMatProd<double>> eigs(op, numerator-2, numerator);
                    // Initialize and compute
                    eigs.init();
                    int nconv = eigs.compute(SortRule::LargestMagn);
                    // Retrieve results
                    if(eigs.info() == CompInfo::Successful){
                        eigval = eigs.eigenvalues();
                        eigvec = eigs.eigenvectors();
                    }

                    //ordinary_edges_fraction = std::count(OrdinaryCells.begin(), OrdinaryCells.end(), 1)
                 //   std::cout << "Eigenvalues found:\n" << eigval << std::endl;
                    std::cout << "\n Eigenvectors nonZeros :\n" << eigvec.nonZeros() << "\n Eigenvectors size :\n" << eigvec.size() << std::endl;


//                        cout << Laplacian_specialFaces << endl;
                }
                if (configuration.at(7)) { // Betti numbers (FB0, FB1 and FB2) of the special faces network
                    //#include"BettiNumbers.h"
                }
                if (configuration.at(8)) { // Grain Laplacian with its spectrum for the grain network
                    //#include"GrainLaplacians.h"
                }
                if (configuration.at(9)) { // Tutte polynomial for the special faces network
                    //#include"FaceTutte.h"
                }
                if (configuration.at(10)) { // Statistical physics module
                    //#include"DCCStatistical.h"
                }
                /// Function which output the results in files
                /// DCCWriter();

                /// End of the Ftype loop

                cout << "Fraction of special faces:" << 1.0 - ordinary_edges_fraction << ";\t Size of S-Face Adjacency matrix \t" << SAM_FacesGraph.nonZeros() << endl;
            } // End of IF: iterator % X == 0
        }while(ordinary_edges_fraction > (1.0 - max_sFaces_fraction));

        OutFaceFile.close(); // 2-Cells related matrices
        OutTJsFile.close(); // 1-Cells statistics output
        OutSFile.close(); // 1-Cells Configurational Entropy output
        //Remove all elements anf free the memory from the probe SpecialCells vector
        SpecialCells.clear();
        SpecialCells.shrink_to_fit();

    } ///End of 'R' type simulations

    else if (stype == 'S') {

    }
    else if (stype == 'I') {

    }
    else if (stype == 'E') {

    }
    else { cout << "ERROR [HAGBsProbability3D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl; return;}


    return;
} /// The end of HAGBsProbability3D()

void HAGBsProbability2D(char* AN, char* AE, char* AF, char* MEN, char* MFE, char* CNumbs, char* config, char* odir, char stype) {
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// odir - source path
// stype - simulation type ('R', 'S', 'I',....)
////////////////////// Matrices initialisation and reading from file //////////////////////
//Triplet T=T(i,j,value), where i and j element's indices in the corresponding matrix
    typedef Triplet<double> Tr; // Eigen library class
    typedef SparseMatrix<int> SpMat; // Declares a column-major sparse matrix type of int
    std::vector<Tr> tripletList; // probe vector of triplets

//Reading from file by the in_streams creation:
/// Reading vector from the file "number_of_cells" the numbers of cells od different types:: vector components: [0] - Nodes, [1] - Edges, [2] - Faces, [3] - Grains
    long unsigned int i = 0, j = 0, value = 0, t_length = 0;
    std::vector<unsigned int> CellNumbs; //number_of_cells.out
    ifstream inCellNumbers(CNumbs);
    if (inCellNumbers.is_open()) { //If the file was successfully open, then
        while(inCellNumbers >> i) CellNumbs.push_back(i);
    } else cout << "The file " << CNumbs << " cannot be read" << endl; //If something goes wrong
    //Screen output for the numbers of cells in the DCC
    cout << "The number of different k-cells in the DCC:" << endl;
    for (int j : CellNumbs) cout << t_length++ << "-cells #\t" << j << endl;
/// Adjacency sparse matrix for nodes
    i = 0, j = 0, value = 0, t_length = 0;
    ifstream inAN(AN);
    if (inAN.is_open()) { //If the file was successfully open, then
        while(!inAN.eof()) {
            inAN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << AN << " cannot be read" << endl; //If something goes wrong
//Sparse AN matrix
    SpMat ANS(CellNumbs.at(0) + 1,CellNumbs.at(0) + 1);
    ANS.setFromTriplets(tripletList.begin(), tripletList.end());

/// Adjacency sparse matrix for edges
    i = 0, j = 0, value = 0; tripletList.clear();
    ifstream inAE(AE);
    if (inAE.is_open()) { //If the file was successfully open, then
        while(!inAE.eof()) {
            inAE >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << AE << " cannot be read" << endl; //If something goes wrong
//Sparse AE matrix
    SpMat AES(CellNumbs.at(1) + 1,CellNumbs.at(1) + 1);
    AES.setFromTriplets(tripletList.begin(), tripletList.end());

/// Adjacency sparse matrix for faces
    i = 0, j = 0, value = 0; tripletList.clear();
    ifstream inAF(AF);
    if (inAF.is_open()) { //If the file was successfully open, then
        while(!inAF.eof()) {
            inAF >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << AF << " cannot be read" << endl; //If something goes wrong
//Sparse AF matrix
    SpMat AFS(CellNumbs.at(2) + 1,CellNumbs.at(2) + 1);
    AFS.setFromTriplets(tripletList.begin(), tripletList.end());

/// Incidence sparse matrix for Edges and Nodes
    i = 0, j = 0, value = 0; tripletList.clear();
    ifstream inMEN(MEN);
    if (inMEN.is_open()) { //If the file was successfully open, then
        while(!inMEN.eof()) {
            inMEN >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << MEN << " cannot be read" << endl; //If something goes wrong
//Sparse MEN matrix
    SpMat ENS(CellNumbs.at(1) + 1,CellNumbs.at(0) + 1);
    ENS.setFromTriplets(tripletList.begin(), tripletList.end());

/// Incidence sparse matrix for Faces and Edges
    i = 0, j = 0, value = 0; tripletList.clear();
    ifstream inMFE(MFE);
    if (inMFE.is_open()) { //If the file was successfully open, then
        while(!inMFE.eof()) {
            inMFE >> i >> j >> value;
            tripletList.push_back(Tr(i, j, value));
        }
    } else cout << "The file " << MFE << " cannot be read" << endl; //If something goes wrong
//Sparse MFE matrix
    SpMat FES(CellNumbs.at(2) + 1,CellNumbs.at(1) + 1);
    FES.setFromTriplets(tripletList.begin(), tripletList.end());

//Remove all elements anf free the memory from the probe vector
    tripletList.clear();
    tripletList.shrink_to_fit();

//////////////////////////////// KINETIC PROCESS //////////////////////////////////////////
    using namespace std;
/////////////////////////////       Random generation       ///////////////////////////////

    if (stype == 'R') { // Pure random generation
        srand (time(NULL)); // initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
        for (int l = 0; l < 1000; l++) {
            long unsigned int New2CellNumb = rand() % CellNumbs.at(2); //Random generation of the boundary number
        }


    } /// End of the Random generation case

    else if (stype == 'S') { // Maximum entropy production

    }
    else if (stype == 'I') { // Ising-like model

    }
    else { cout << "ERROR [HAGBsProbability2D] : unknown simulation type - please replace with 'R', 'S' or 'I'..!" << endl; return;}

    return;
} /// The end of HAGBsProbability2D()

/// ================================== Related functions ==================================///
void confCout(vector<bool> configuration, int number_of_types) {
    cout << "The number of special face types:\t" << number_of_types << endl;
    cout << "Nodes types statistics, indices and configuration entropy:     "; if (configuration.at(0) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges types statistics, indices and configuration entropy:     "; if (configuration.at(1) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Faces types statistics and structural indices:                 "; if (configuration.at(2) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain types statistics, indices and configuration entropy:     "; if (configuration.at(3) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Nodes Laplacian with its spectrum for the nodes network:       "; if (configuration.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Edges Laplacian with its spectrum for the edges network:       "; if (configuration.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Face Laplacian with its spectrum for the special faces network:"; if (configuration.at(6) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Betti numbers (FB0, FB1 and FB2) of the special faces network: "; if (configuration.at(7) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain Laplacian with its spectrum for the grain network:       "; if (configuration.at(8) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Tutte polynomial for the special faces network:                "; if (configuration.at(9) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Statistical physics module:                                    "; if (configuration.at(10) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    return;
}
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, Eigen::SparseMatrix<double> const& SM1) {
    Eigen::SparseMatrix<double> res = SM1;
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

/// Heap
// (1)    cout << std::numeric_limits<long unsigned int>::max() << endl;

// (2) ofstream NONStructGlCoordOut; NONStructGlCoordOut.open("C:VAnRes\\(MFE21)GNC_NONStructured.txt", ios::trunc);

//(3) ARCHIVE testing output  //   cout << "The number of columns in the created matrices:\t" << "\tANS --\t" << ANS.cols() << "\tAES --\t" << AES.cols()<< "\tAFS --\t" << AFS.cols() << "\tAGS --\t" << AGS.cols() << "\tENS --\t" << ENS.cols() << "\tFES --\t" << FES.cols() << "\tGFS --\t" << GFS.cols() << endl;


