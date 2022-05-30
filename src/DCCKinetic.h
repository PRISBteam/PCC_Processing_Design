///================================ DCC Kinetic module =============================================================///
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
//#include "Spectra/GenEigsSolver.h"
//#include "Spectra/SymEigsSolver.h"
//#include "Spectra/MatOp/SparseGenMatProd.h"
///-------------------------------------

using namespace std;
using namespace Eigen;
//using namespace Spectra;

void HAGBsKinetic3D(char* AN, char* AE, char* AF, char* AG, char* MEN, char* MFE, char* MGF, char* CNumbs, char* config, char* odir, char stype) {
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
//    ofstream OutFaceFile; // Output stream for variables related with special 2-Cells in DCC
//    OutFaceFile.open(odir + "FaceAdjacency.txt"s, ios::trunc);
//    OutFaceFile<< "Adjacency Matrix of all Special Faces" << endl;
//    OutFaceFile.close();

////////////////////////////////////// KINETIC PROCESS ///////////////////////////////////////////////
    if (stype == 'F') { // Random generation case

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
        double max_sFaces_fraction = 0.5;
        /// Step for output (structural analysis and data output will be performed after each output_step calculation steps (= number of newly converted elements))
        int output_step = 300;
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

                if (configuration.at(1)) {  // External force effect
                }
                if (configuration.at(2)) {  // Internal gradient force effect
                }
                if (configuration.at(3)) { // Chi-factor - orientation effect
                }
                if (configuration.at(4)) { // GBs and TJs types effect
                }
                if (configuration.at(5)) { // Grain boundary dislocations effect
                }

//                cout << "Average Strain:" << MStrain << ";\t Size of S-Face Adjacency matrix \t" << SAM_FacesGraph.nonZeros() << endl;
            } // End of analysis and output iterator ( IF: iterator % X == 0 )
        }while(ordinary_faces_fraction > (1.0 - max_sFaces_fraction)); /// End of the Random generation process

/// Closing and deleting

    } ///End of 'F' type simulations

    else { cout << "ERROR [HAGBsKinetic3D] : unknown simulation type - please replace with 'F' " << endl; return;}


    return;
} /// The end of HAGBsKinetic3D()

/// ================================== Related functions ==================================///
std::vector<int> confCout(char* config) {
    std::string line;
    std:vector<int> res;
    ifstream inConf(config);
    if (inConf) { //If the file was successfully open, then
        while(getline(inConf, line, '\n')) {
            if(line.at(0) == '#') res.push_back(1); // 1 and # means accept - the parameter will be calculated
            if(line.at(0) == '%') res.push_back(0); // 0 and % means ignore - the parameter will not be calculated
        }
    } else cout << "The file " << config << " cannot be read" << endl; // If something goes wrong

    cout << "External force effect:                                         "; if (res.at(1) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Internal gradient force (between adjacent grains) effect:      "; if (res.at(2) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Chi-factor (Xi) effect:                                        "; if (res.at(3) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "GBs and TJs types effect:                                      "; if (res.at(4) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;
    cout << "Grain boundary dislocations effect:                            "; if (res.at(5) == 1) cout << "\t[" << "On" << "]\t" << endl; else cout << "\t[" << "Off" << "]\t" << endl;

    return res;
}
/**
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
**/
/// Heap