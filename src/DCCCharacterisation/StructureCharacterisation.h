/// Attached user defined C++ libraries:
///-------------------------------------
///-------------------------------------
#include "TJsLab.h"
#include "LaplaciansLab.h"
///-------------------------------------

///Structure characterisation tool
//unsigned int DCC_StructureCharacterisation(std::vector<unsigned int> const& S_Vector, std::vector<unsigned int> const& s_faces_sequence, std::vector<double> const& configuration, std::vector<unsigned int> const& CellNumbs, std::vector<char*> const paths, char* odir) {
int DCC_StructureCharacterisation(std::vector<unsigned int> &S_Vector, std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> &c_faces_sequence, std::vector<double> configuration, std::vector<unsigned int> &CellNumbs, std::vector<char*> const paths, char* odir, std::ofstream& OutFLfile, std::ofstream& OutElCondfile, std::ofstream&  OutSFile)
{
    vector<double> jcF_types_fractions, jF_types_fractions;
    double SFace_Entropy_Median = 0, SFace_Entropy_Skrew = 0, cFace_Entropy_Median = 0, cFace_Entropy_Skrew = 0;
    double SFinformativeness = 0, CFinformativeness = 0;
    double Svn = 0.0, Cvn = 0.0; // Von-Neumann entropy for Special Faces (Svn) and Cracked Faces (Cvn)
// cout << "START of DCC Structure Characterisation Module" << endl;
    SpMat SpAM_SpecFaces(s_faces_sequence.size(), s_faces_sequence.size()), SFace_Laplacian(s_faces_sequence.size(), s_faces_sequence.size()),
          SpAM_CrackFaces(c_faces_sequence.size(), c_faces_sequence.size()), CFace_Laplacian(c_faces_sequence.size(),c_faces_sequence.size()),
          Sym_SFace_Laplacian(s_faces_sequence.size(), s_faces_sequence.size()), RW_SFace_Laplacian(
            s_faces_sequence.size(), s_faces_sequence.size()),
            Ids(s_faces_sequence.size(), s_faces_sequence.size());
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

    /// "cracked" Adjacency matrix of special faces
    unsigned int num = 0;
    SFaces_Triplet_list.clear();

    for (auto sit: s_faces_sequence) { // take each special Face one by one
        for (unsigned int k = 0; k < CellNumbs.at(2); ++k) { //  Loop over all the Faces in the DCC
            if (k != sit && AFS.coeff(sit, k) == 1) { // if (1) not the same element and (2) is a neighbour
                auto itr = std::find(s_faces_sequence.begin(), s_faces_sequence.end(), k);
                if (itr != s_faces_sequence.end()) { // if the neighbour is also a special Face
                    unsigned int neighbour = std::distance(s_faces_sequence.begin(), itr); // receive neighbour number
//cout << s_faces_sequence.at(std::distance(s_faces_sequence.begin(), itr)) << " ==== " << std::distance(s_faces_sequence.begin(), itr) << endl;

                    /// Fracture effect: if two special faces contain the third cracked neighbour then they have no connection between each other
                    bool accept = 1;
                    for (auto m: c_faces_sequence) {
                             if (m != sit && m != k && AFS.coeff(sit, m) == 1 && AFS.coeff(k, m) == 1) { // if the new third Face is the neighbour of two others in the junction
                                 accept = 0;  break;
                             }
                    } // end  for (auto m: c_faces_sequence)

                    /// NEW ELEMENT TO THE TRIPLET LIST
                    if (accept == 1) {
                                 SFaces_Triplet_list.push_back(Tr(num, neighbour, 1));
                             } // end if (accept)

                } // if (itr != s_faces_sequence.end())
            } // end if (k != sit && AFS.coeff(sit, k) == 1)
        } // end for (unsigned int k = 0; k < CellNumbs.at(2); ++k)
        ++num;
    } // end for ( s_faces_sequence) | then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

    cout << "SFaces_Triplet_list size =  \t" << SFaces_Triplet_list.size() << endl;
    cout << "special faces fraction = \t" << s_faces_sequence.size()/ (double) CellNumbs.at(2) << "\t crack faces fraction = \t" << c_faces_sequence.size()/ (double) CellNumbs.at(2) << endl;

    /// Sparse Face adjacency matrix (already with cracks)
    if (SFaces_Triplet_list.size() > 50) {
        SpAM_SpecFaces.setFromTriplets(SFaces_Triplet_list.begin(), SFaces_Triplet_list.end());
        SpAM_SpecFaces = 0.5 * (SpAM_SpecFaces +
                                SparseMatrix<double>(SpAM_SpecFaces.transpose())); // Full matrix instead of triagonal

//        cout << "SpAM_SpecFaces size =  \t" << SpAM_SpecFaces.size() << endl;

/// Degree vector and matrix
        vector<unsigned int> SFace_degrees;
        for (unsigned int i = 0; i < SpAM_SpecFaces.cols(); i++) { //columns
            unsigned int degree_Fcounter = 0;
            for (unsigned int j = 0; j < SpAM_SpecFaces.rows(); j++) { // rows
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
//        Ids.setIdentity();
/// Special Faces Laplacian matrix
        SFace_Laplacian = SFDegree - SpAM_SpecFaces; // L = D - A

/// === Density matrix of states =======>
/*
double a_sum = 0.0;
for (unsigned int k = 0; k < SpAM_SpecFaces.rows(); ++k)
    for (unsigned int l = 0; l < SpAM_SpecFaces.cols(); ++l)
        a_sum += SpAM_SpecFaces.coeff(k,l);

    SpMat Face_Rho = SFace_Laplacian/ a_sum;
    /// Von-Neumann entropy for Special Faces (cracked) Laplacian
    MatrixXd Face_Rho_Dense;
        Face_Rho_Dense = MatrixXd(Face_Rho);
        Face_Rho_Dense = Face_Rho_Dense*Face_Rho_Dense.array().log().matrix();
// REPAIR       cout << "Face_Rho_Dense after == " << Face_Rho_Dense.nonZeros() << endl;

        Svn = Face_Rho_Dense.diagonal().sum();
        Face_Rho_Dense.resize(0,0);
*/
Svn = 0; ////////////////-
// Symmetric S-Faces Laplacian
/// Left Random-Walk normalized Laplacian matrix
//RW_SFace_Laplacian = Id - SFDegree.cwiseInverse() * SpAM_SpecFaces; // Ls = D^-1/2 L
        SFace_degrees.clear();
        SFace_degrees.shrink_to_fit();

    } // if (SFaces_Triplet_list.size() > M)

/// ====== Adjacency matrix and Laplacian for cracked boundaries ONLY =================================>
///====================================================================
    vector<Tr> CFaces_Triplet_list;

    /// Adjacency matrix of cracked Faces
    unsigned int cnum = 0;
    CFaces_Triplet_list.clear();
    for (auto sit: c_faces_sequence) { // take each special Face one by one
        for (unsigned int k = 0; k < CellNumbs.at(2); ++k) { //  Loop over all the Faces in the DCC
            if (k != sit && AFS.coeff(sit, k) == 1) { // if (1) not the same element and (2) is a neighbour
                auto itr = std::find(c_faces_sequence.begin(), c_faces_sequence.end(), k);
                if (itr != c_faces_sequence.end()) { // if the neighbour is also a special Face
                    unsigned int neighbour = std::distance(c_faces_sequence.begin(), itr); // receive neighbour number

                    /// NEW ELEMENT TO THE TRIPLET LIST
                        CFaces_Triplet_list.push_back(Tr(cnum, neighbour, 1));
                } // if (itr != c_faces_sequence.end())
            } // end if (k != sit && AFS.coeff(sit, k) == 1)
        } // end for (unsigned int k = 0; k < CellNumbs.at(2); ++k)
        ++cnum;
    } // end for ( s_faces_sequence) | then we add this element to the new Special Faces Matrix with the number {numerator, new number of the neighbour}

    /// Sparse CRACKED Face adjacency matrix
    if (CFaces_Triplet_list.size() > 50) {
        SpAM_CrackFaces.setFromTriplets(CFaces_Triplet_list.begin(), CFaces_Triplet_list.end());
        SpAM_CrackFaces = 0.5 * (SpAM_CrackFaces +
                                SparseMatrix<double>(SpAM_CrackFaces.transpose())); // Full matrix instead of triagonal
/// Degree vector and matrix
        vector<unsigned int> CFace_degrees;
        for (unsigned int i = 0; i < SpAM_CrackFaces.cols(); i++) { //columns
            unsigned int degree_CFcounter = 0;
            for (unsigned int j = 0; j < SpAM_CrackFaces.rows(); j++) { // rows
                if (SpAM_CrackFaces.coeff(j, i) == 1 && i != j)
                    degree_CFcounter++;
            }
            CFace_degrees.push_back(degree_CFcounter);
        }

//Creation of the S-Face degree matrix
        SpMat SFCDegree(CFace_degrees.size(), CFace_degrees.size());
        SFCDegree.setIdentity();
        unsigned int cnumerator_degree = 0;
        for (double cnum: CFace_degrees) {
            SFCDegree.coeffRef(cnumerator_degree, cnumerator_degree) = CFace_degrees.at(cnumerator_degree);
            cnumerator_degree++;
        }

/// Special Faces Laplacian matrix
        CFace_Laplacian = SFCDegree - SpAM_CrackFaces; // L = D - A

        /// === Density matrix of Cracked Faces states =======>
/*
        double ac_sum = 0.0;
        for (unsigned int k = 0; k < SpAM_CrackFaces.rows(); ++k)
            for (unsigned int l = 0; l < SpAM_CrackFaces.cols(); ++l)
                ac_sum += SpAM_CrackFaces.coeff(k,l);

        SpMat Face_cRho = CFace_Laplacian/ ac_sum;
        /// Von-Neumann entropy for Special Faces (cracked) Laplacian
        MatrixXd Face_Rho_cDense;
        Face_Rho_cDense = MatrixXd(Face_cRho);
        Face_Rho_cDense = Face_Rho_cDense*Face_Rho_cDense.array().log().matrix();

        Cvn = Face_Rho_cDense.diagonal().sum();
        Face_Rho_cDense.resize(0,0);
        */
Cvn = 0; ////////////////-


        // Symmetric S-Faces Laplacian
/// Left Random-Walk normalized Laplacian matrix
// Identity matrix
//        Id.setIdentity();
//RW_SFace_Laplacian = Id - SFCDegree.cwiseInverse() * SpAM_CrackFaces; // Ls = D^-1/2 L
        CFace_degrees.clear();
        CFace_degrees.shrink_to_fit();

    } // if (SFaces_Triplet_list.size() > M)


/// ============== Generation of the output streams ==============================================>>
    ofstream OutFaceFile; // Output stream for variables related with special 2-Cells in DCC
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
if (configuration.at(4)) { // Nodes types statistics, indices and configuration entropy
//#include "QPsLab.h"
}
if (configuration.at(5)) { /// Edges types statistics, indices and configuration entropy
    //map<unsigned int, unsigned int> Edges_Types_Map; // Map: [Edge number] --> [Enge type]
    /// Statistics of Edges
    // Special
    EdgesStat( s_faces_sequence, CellNumbs, FES, odir, SFace_Entropy_Median, SFace_Entropy_Skrew, SFinformativeness, jF_types_fractions);
    // Crack
    EdgesStat( c_faces_sequence, CellNumbs, FES, odir, cFace_Entropy_Median, cFace_Entropy_Skrew, CFinformativeness,jcF_types_fractions);

    OutSFile << s_faces_sequence.size() / (double) CellNumbs.at(2) << " " << c_faces_sequence.size() / (double) CellNumbs.at(2) << " " << Svn << " " << jF_types_fractions.at(0) << " " << jF_types_fractions.at(1) << " " << jF_types_fractions.at(2) << " " << jF_types_fractions.at(3) << " " << ( jF_types_fractions.at(1) + 2.0*jF_types_fractions.at(2) + 3.0*jF_types_fractions.at(3))/ 3.0 << " " << SFace_Entropy_Median + SFace_Entropy_Skrew << " " << SFace_Entropy_Median << " " << SFace_Entropy_Skrew << " " << SFinformativeness << " "
                << " " << Cvn << " " << jcF_types_fractions.at(0)<< " " << jcF_types_fractions.at(1)<< " " << jcF_types_fractions.at(2)<< " " << jcF_types_fractions.at(3) << " " << ( jcF_types_fractions.at(1) + 2.0*jcF_types_fractions.at(2) + 3.0*jcF_types_fractions.at(3))/ 3.0 << " " << cFace_Entropy_Median + cFace_Entropy_Skrew << " " << cFace_Entropy_Median << " " << cFace_Entropy_Skrew << " " << CFinformativeness << endl;
}
if (configuration.at(6)) { // Faces types statistics and structural indices
//#include "FacesLab.h"
}
if (configuration.at(7)) { // Grains types statistics, indices and configuration entropy
//#include "GrainsLab.h"
}
if (configuration.at(8)) { // Nodes Laplacian with its spectrum for the nodes network
//#include"NodeLaplacians.h"
}
if (configuration.at(9)) { // Edges Laplacian with its spectrum for the edges network
//#include"EdgeLaplacians.h"
}
if (configuration.at(10)) { // Face Laplacian with its spectrum for the special faces network

        vector<double> Cface_spectrum, inverse_Cface_spectrum;
        // if (CFace_Laplacian.rows() > 50) {
        if (CFaces_Triplet_list.size() > 100) {
                Cface_spectrum = FaceLaplacian(CFace_Laplacian, CellNumbs);
            inverse_Cface_spectrum = Cface_spectrum;
            for (unsigned int i = 0; i < Cface_spectrum.size(); ++i) if ( Cface_spectrum.at(i) >= pow(10,-5)) inverse_Cface_spectrum.at(i) = 1.0 / Cface_spectrum.at(i);
        } // end if (CFace_Laplacian.rows() > 50)

    /// (cracked) Face electrical conductivity
    vector<double> Lface_spectrum, inverse_Lface_spectrum;
        // if (SFace_Laplacian.rows() > 50) {
        if (SFaces_Triplet_list.size() > 100) {
                Lface_spectrum = FaceLaplacian(SFace_Laplacian, CellNumbs);
                inverse_Lface_spectrum = Lface_spectrum;
                for (unsigned int i = 0; i < Lface_spectrum.size(); ++i) if ( Lface_spectrum.at(i) >= pow(10,-5)) inverse_Lface_spectrum.at(i) = 1.0 / Lface_spectrum.at(i);

                cout << "Laplacian size :: " << SFace_Laplacian.cols() << "\tNumber of Faces\t" << CellNumbs.at(2) << endl;
                //system("pause");

                /// Output to file Special Face Laplacian eigenvalues
                if (OutFLfile) {
                    OutFLfile << endl << "Special Faces fraction is equal to\t " << MAX_sfaces_fraction << endl << endl;
                    OutFLfile << SFace_Laplacian.cols() / CellNumbs.at(2) << " ";
                    for (unsigned int i = 0; i < Lface_spectrum.size(); ++i) OutFLfile << Lface_spectrum.at(i) << " "; OutFLfile << endl;

                } else cout << "Error: No such a directory for\t" << odir + "FaceLaplacian.txt"s << endl;
        } // if (SFace_Laplacian.rows()

    OutElCondfile << s_faces_sequence.size() / (double) CellNumbs.at(2) << " " << c_faces_sequence.size() / (double) CellNumbs.at(2) << " " << std::accumulate(inverse_Lface_spectrum.begin(), inverse_Lface_spectrum.end(), 0) << " " << 1.0 / ( Lface_spectrum.size() * std::accumulate(inverse_Lface_spectrum.begin(), inverse_Lface_spectrum.end(), 0) ) << " " << std::count(Lface_spectrum.begin(), Lface_spectrum.end(), 0.0) << " " << SFinformativeness << " "
                    << std::accumulate(inverse_Cface_spectrum.begin(), inverse_Cface_spectrum.end(), 0) << " " << 1.0 / ( Cface_spectrum.size() * std::accumulate(inverse_Cface_spectrum.begin(), inverse_Cface_spectrum.end(), 0) ) << " " << std::count(Cface_spectrum.begin(), Cface_spectrum.end(), 0.0) << " " << CFinformativeness << endl;
} /// if (configuration.at(10))
if (configuration.at(11)) { // Grain Laplacian with its spectrum for the grain network
//#include"GrainLaplacians.h"
}
if (configuration.at(12)) { // Statistical physics module
//#include"DCCStatistical.h"
}
    /// Vectors deletion
    SFaces_Triplet_list.clear();
    SFaces_Triplet_list.shrink_to_fit();
    SFace_Laplacian.makeCompressed();
    Sym_SFace_Laplacian.makeCompressed();
    RW_SFace_Laplacian.makeCompressed();
    Ids.makeCompressed();
    //SpAM_SpecFaces.makeCompressed();

    return SpAM_SpecFaces.nonZeros();
} /// End of the Structure_Characterisation function

/// Archive