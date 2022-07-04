/// Attached user defined C++ libraries:
///-------------------------------------
///-------------------------------------
#include "TJsLab.h"
#include "FaceLaplacians.h"
///-------------------------------------

///Structure characterisation tool
//unsigned int DCC_StructureCharacterisation(std::vector<unsigned int> const& S_Vector, std::vector<unsigned int> const& s_faces_sequence, std::vector<double> const& configuration, std::vector<unsigned int> const& CellNumbs, std::vector<char*> const paths, char* odir) {
int DCC_StructureCharacterisation(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, std::vector<double> configuration, std::vector<unsigned int> &CellNumbs, std::vector<char*> const paths, char* odir) {

// cout << "START of DCC Structure Characterisation Module" << endl;
    SpMat SpAM_SpecFaces(s_faces_sequence.size(), s_faces_sequence.size()), SFace_Laplacian(s_faces_sequence.size(),
                                                                                            s_faces_sequence.size()),
            Sym_SFace_Laplacian(s_faces_sequence.size(), s_faces_sequence.size()), RW_SFace_Laplacian(
            s_faces_sequence.size(), s_faces_sequence.size()),
            Id(s_faces_sequence.size(), s_faces_sequence.size());
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values // odir - source path
    /// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
//    SpMat ANS(CellNumbs.at(0), CellNumbs.at(0)), AES(CellNumbs.at(1), CellNumbs.at(1)), AFS(CellNumbs.at(2), CellNumbs.at(2)), AGS(CellNumbs.at(3), CellNumbs.at(3));
    SpMat AFS(CellNumbs.at(2), CellNumbs.at(2));

//    ANS = SMatrixReader(paths.at(0), (CellNumbs.at(0)), (CellNumbs.at(0))); //all Nodes
//    ANS = 0.5 * (ANS + SparseMatrix<double>(ANS.transpose())); // Full matrix instead of triagonal
//    AES = SMatrixReader(paths.at(1), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
//    AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose())); // Full matrix instead of triagonal
    AFS = SMatrixReader(paths.at(2), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal
//    AGS = SMatrixReader(paths.at(3), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
//    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
    SpMat FES(CellNumbs.at(1), CellNumbs.at(2));
//    SpMat ENS(CellNumbs.at(0), CellNumbs.at(1)), FES(CellNumbs.at(1), CellNumbs.at(2)), GFS(CellNumbs.at(2), CellNumbs.at(3));
    //    ENS = SMatrixReader(paths.at(4), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    FES = SMatrixReader(paths.at(5), (CellNumbs.at(1)), (CellNumbs.at(2))); //all Edges-Faces
//    GFS = SMatrixReader(paths.at(6), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
/// Set simulation configuration from configuration file :: the number of special face types and calculating parameters.
    int number_of_types = (int) configuration.at(0);

/// ===== Creation of the Special Faces triplet list =====>
// Another key point: creation of a triplet list SFaces_Triplet_list for the sparse matrix of special Faces
    vector<Tr> SFaces_Triplet_list;
    // REPAIR cout << AFS.cols() <<"  "<< AFS.rows() << "  " << CellNumbs.at(2) << endl;
    // REPAIR for (auto sit : s_faces_sequence) cout << "  "<< sit; cout << endl;

    unsigned int num = 0;
    for (auto sit: s_faces_sequence) {
        for (unsigned int k = 0; k < CellNumbs.at(2); k++) //  Loop over all the Faces in the DCC
            if (k != sit && AFS.coeff(sit, k) == 1) { // if (1) not the same element
                auto itr = std::find(s_faces_sequence.begin(), s_faces_sequence.end(), k);
                // (2) we find an adjacent element (neighbour)
                if (itr != s_faces_sequence.end()) {
                    unsigned int neighbour = std::distance(s_faces_sequence.begin(), itr);
                    SFaces_Triplet_list.push_back(Tr(num, neighbour, 1));
                } // if (s_faces_sequence)
            } // if ASF.coeff
        ++num;
    } // for(s_faces_sequence) | then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}
    // Sparse Face adjacency matrix
    SpAM_SpecFaces.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());
    SpAM_SpecFaces = 0.5 * (SpAM_SpecFaces +
                            SparseMatrix<double>(SpAM_SpecFaces.transpose())); // Full matrix instead of triagonal
    SFaces_Triplet_list.clear();

/// Degree vector and matrix
    vector<unsigned int> SFace_degrees, columns(s_faces_sequence.size());
    for (unsigned int i = 0; i < s_faces_sequence.size(); i++) { //columns
        unsigned int degree_Fcounter = 0;
        for (unsigned int j = 0; j < s_faces_sequence.size(); j++) { // rows
            if (SpAM_SpecFaces.coeff(j, i) == 1 && i != j)
                degree_Fcounter++;
        }
        SFace_degrees.push_back(degree_Fcounter);
    }

//Creation of the S-Face degree matrix
    SpMat SFDegree(SFace_degrees.size(), SFace_degrees.size());
    SFDegree.setIdentity();
    unsigned int numerator_degree = 0;
    for (double num: SFace_degrees) {
        SFDegree.coeffRef(numerator_degree, numerator_degree) = SFace_degrees.at(numerator_degree);
        numerator_degree++;
    }

// Identity matrix
    Id.setIdentity();
/// Special Faces Laplacian matrix
    SFace_Laplacian = SFDegree - SpAM_SpecFaces; // L = D - A

// Symmetric S-Faces Laplacian
/// Left Random-Walk normalized Laplacian matrix
//RW_SFace_Laplacian = Id - SFDegree.cwiseInverse() * SpAM_SpecFaces; // Ls = D^-1/2 L

/// Generation of the output streams
    ofstream OutFaceFile; // Output stream for variables related with special 2-Cells in DCC
    ofstream OutFLfile; // Special Face Laplacian
    ofstream OutRWFLfile; // Random Walker Special Face Laplacian

/// MAX special_faces_fraction
    double MAX_sfaces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);

/// Output to file Special Faces Adjacency matrix
    OutFaceFile.open(odir + "FaceAdjacency.txt"s, ios::trunc);
    if (OutFaceFile) {
        OutFaceFile << "Adjacency Matrix of all Special Faces" << endl;
        OutFaceFile << endl << "max Special Faces fraction is equal to\t " << MAX_sfaces_fraction << endl
                    << endl;                 // Output of special_faces_fraction to files
        OutFaceFile << SpAM_SpecFaces << endl; //  Sparse matrix output
        //OutFaceFile << MatrixXd(SpAM_SpecFaces) << endl; //  Dense matrix output
        OutFaceFile.close();
    } else cout << "Error: No such a directory for\t" << odir + "FaceAdjacency.txt"s << endl;

/// Output to file Special Faces Laplacian matrix
    OutFLfile.open(odir + "FaceLaplacian.txt"s, ios::trunc);
    if (OutFLfile) {
        OutFLfile << "Laplacian Matrix of all Special Faces" << endl;
        OutFLfile << endl << "Special Faces fraction is equal to\t " << MAX_sfaces_fraction << endl << endl;
                OutFLfile << SFace_Laplacian << endl; //  Sparse matrix output
//        OutFLfile << MatrixXd(SFace_Laplacian) << endl; // Dense matrix output
        OutFLfile.close();
    } else cout << "Error: No such a directory for\t" << odir + "FaceLaplacian.txt"s << endl;

/// Output to file Special Faces order
    OutRWFLfile.open(odir + "SpecialGrainBoundaries.txt"s, ios::trunc);
    if (OutRWFLfile) {
    OutRWFLfile << "Global numbers (in DCC) of special grain boundaries with the fraction " << MAX_sfaces_fraction
                << endl;
// for (auto sit : SpecialCellMap) OutRWFLfile << sit.second << endl;
    for (auto vit: s_faces_sequence) OutRWFLfile << vit << endl;
// Dense matrix output:         OutRWFLfile << MatrixXd(RW_Face_Laplacian) << endl;
    OutRWFLfile.close();
    } else cout << "Error: No such a directory for\t" << odir + "SpecialGrainBoundaries.txt"s << endl;

/// =========== Analysis tools ==============>
if (configuration.at(1)) { // Nodes types statistics, indices and configuration entropy
//#include "QPsLab.h"
}
if (configuration.at(2)) { /// Edges types statistics, indices and configuration entropy
    //map<unsigned int, unsigned int> Edges_Types_Map; // Map: [Edge number] --> [Enge type]
    /// Statistics of Edges
    EdgesStat(S_Vector, s_faces_sequence, CellNumbs, FES, odir);
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
    /// Vectors deletion
    SFaces_Triplet_list.clear();
    SFaces_Triplet_list.shrink_to_fit();
    SFace_degrees.clear();

    SFace_degrees.shrink_to_fit();
    columns.clear();
    columns.shrink_to_fit();
    SFace_Laplacian.makeCompressed();
    Sym_SFace_Laplacian.makeCompressed();
    RW_SFace_Laplacian.makeCompressed();
    Id.makeCompressed();
    //SpAM_SpecFaces.makeCompressed();

    return SpAM_SpecFaces.nonZeros();
} /// End of the Structure_Characterisation function




/// Archive
/*
/// Declaration of FUNCTIONS, see the function bodies at the end of file.
Eigen::SparseMatrix<double> SMatrixReader(char* SMpath, unsigned int Rows, unsigned int Cols); // The function read any matrices from lists of triplets and create corresponding sparse matrices
std::vector<unsigned int> VectorReader(char* FilePath); // The function read any matrices from lists of triplets and create corresponding sparse matrices
//std::vector<double> confCount(char* config); // Configuration's of the problem output read from file "config"
//vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, map<unsigned int, unsigned int> Edges_TypesMap);
unsigned int NewCellNumb_R(unsigned int SCellsNumb); // Random generation of a 2-Cell number
//    unsigned int Structure_Characterisation(unsigned long const& numerator, vector<unsigned int> const& CellNumbs, vector<Eigen::Triplet<double>> const& SFaces_Triplet_list, map<unsigned int, unsigned int> const& SpecialCellMap, vector<unsigned int> const& special_faces_sequence, Eigen::SparseMatrix<double> const& FES, char* odir, double special_faces_fraction, vector<double> const& configuration); // The whole characterisation module with the program output

//////////////////////=============== Matrices initialisation part =============== //////////////////////
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
OutSFile << "Strain\t" <<"\t Configurational [1-Cells] Entropy \t" <<"\t Median Face Entropy \t" <<"\t Skrew Face Entropy \t" << endl;

/// Another key point: creation of a triplet list SFaces_Triplet_list for the sparse matrix of special faces
for(unsigned int j = 0; j < CellNumbs.at(2); j++) //  Loop over all the Faces in the DCC
if(j != sit && AFS.coeff(sit,j) == 1 && (SpecialCellMap.find(j) != SpecialCellMap.end())) // if (1) not the same element (2) we find an adjacent element (neighbour) (3)! this element already added in the special faces map (= was randomly chosen)
{ SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); } // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}
//         { SFaces_Triplet_list.push_back(Tr(sit,SpecialCellMap[j],1)); SFaces_Triplet_list.push_back(Tr(SpecialCellMap[j], sit,1));} // then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}


*/