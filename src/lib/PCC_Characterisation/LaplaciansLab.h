///================================ PCC Laplacian laboratory =============================================================///
///=======================================================================================================================///
/** This sub file contain functions allowing to calculate the combinatorial Laplacian of the special k-cells graph with its spectrum and eigenvectors **/
///=======================================================================================================================///

/// # 0 # reduced laplacians (vector of their sparse matrices)
// std::vector<Eigen::SparseMatrix<double>> ReducedFaceLaplacians(CellsDesign &cells_design, std::vector<unsigned int> const &CellNumbs) { // OSM - operator's sparse matrix
std::vector<Eigen::SparseMatrix<double>> ReducedFaceLaplacians(CellsDesign &cells_design) { // OSM - operator's sparse matrix

std::vector<Eigen::SparseMatrix<double>> laplacians_vector; // function output: { ENM_laplacian, FEM_laplacian, GFM_laplacian }
// AN - Nodes (Vertices or 0-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AE - Edges (Triple Junctions or 1-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AF - Faces (Grain Boundaries or 2-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// AG - Grains (Volumes or 3-Cells) sparse adjacency matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MEN - Edges - Nodes (1-Cells to 0-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MFE - Faces - Edges (2-Cells to 1-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values
// MGF - Grains - Faces (3-Cells to 2-Cells) sparse incidence matrix  in the form of three column {i, j, value}, where i and j are the indices of elements with non-zero values // odir - source path

/// Adjacency sparse matrix for nodes /// Adjacency sparse matrix for edges /// Adjacency sparse matrix for faces /// Adjacency sparse matrix for grains
//    SpMat ANS(CellNumbs.at(0), CellNumbs.at(0)), AES(CellNumbs.at(1 + (dim - 3)), CellNumbs.at(1 + (dim - 3))), AFS(CellNumbs.at(2 + (dim - 3)), CellNumbs.at(2 + (dim - 3))), AGS(CellNumbs.at(3 + (dim - 3)), CellNumbs.at(3 + (dim - 3)));
//    ANS = SMatrixReader(paths.at(0), (CellNumbs.at(0)), (CellNumbs.at(0))); //all Nodes
//    ANS = 0.5 * (ANS + SparseMatrix<double>(ANS.transpose())); // Full matrix instead of triagonal
//    AES = SMatrixReader(paths.at(1 + (dim - 3)), (CellNumbs.at(1 + (dim - 3))), (CellNumbs.at(1 + (dim - 3)))); //all Edges
//    AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose())); // Full matrix instead of triagonal
//    AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2 + (dim - 3))), (CellNumbs.at(2 + (dim - 3)))); //all Faces
//    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose())); // Full matrix instead of triagonal
//    AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3 + (dim - 3))), (CellNumbs.at(3 + (dim - 3)))); //all Volumes
//    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal

/// Incidence sparse matrix for Edges and Nodes /// Incidence sparse matrix for Faces and Edges /// Incidence sparse matrix for Grains and Faces
    SpMat ENS(CellNumbs.at(0), CellNumbs.at(1)), FES(CellNumbs.at(1 + (dim - 3)), CellNumbs.at(2 + (dim - 3))), GFS(CellNumbs.at(2 + (dim - 3)), CellNumbs.at(3 + (dim - 3)));

    ENS = SMatrixReader(paths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    FES = SMatrixReader(paths.at(5 + (dim - 3)), (CellNumbs.at(1 + (dim - 3))), (CellNumbs.at(2 + (dim - 3)))); //all Edges-Faces
    GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2 + (dim - 3))), (CellNumbs.at(3 + (dim - 3)))); //all Faces-Grains

/// Triplet lists
    vector<Tr> sEN_triplet_list, sFE_triplet_list, sGF_triplet_list;

/// Sparse matrices (Eigen/ Spectra)
    SpMat ENM_SpecialCells, FEM_SpecialCells, GFM_SpecialCells, L0_laplacian, L1_laplacian, L2_laplacian;

/// Special face-related PCC elements
    std::vector<unsigned int> node_fsequence_vector, edge_fsequence_vector, face_sequence_vector, grain_fsequence_vector; // number of different sequences based on the face design vector
    std::vector<int> node_fstate_vector(CellNumbs.at(0),0), edge_fstate_vector(CellNumbs.at(1 + (dim - 3)),0),
                            face_state_vector(CellNumbs.at(2 + (dim - 3)),0), grain_fstate_vector(CellNumbs.at(3 + (dim - 3)),0); // number of different sequences based on the face design vector

    unsigned int node_numerator = 0, edge_numerator = 0, face_numerator = 0, grain_numerator = 0;
    face_sequence_vector = cells_design.Get_f_sequence(); // face sequence from cells_design after Processing
    /// faces loop
    for (auto  itf = face_sequence_vector.begin(); itf != face_sequence_vector.end(); ++itf) {
        // grains
        for (unsigned int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g) //  Loop over all the Grains in the PCC
                if (GFS.coeff(*itf, g) != 0) {
                    grain_fsequence_vector.push_back(g);

                    if (GFS.coeff(*itf, g) == 1) {
                        grain_fstate_vector.at(g) = 1;
                        sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf),grain_numerator++,+1));
                    }
                    else if (GFS.coeff(*itf, g) == -1) {
                        grain_fstate_vector.at(g) = -1;
                        sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf), grain_numerator++, -1));
                    }
                } // end if (GFS.coeff(*itf, g) == 1)

        for (unsigned int e = 0; e < CellNumbs.at(1 + (dim - 3)); ++e) //  Loop over all the Edges in the PCC
            if (FES.coeff(e, *itf) != 0) {
                edge_fsequence_vector.push_back(e);

                if (FES.coeff(e, *itf) == 1) {
                    edge_fstate_vector.at(e) = 1;
                    sFE_triplet_list.push_back(Tr(edge_numerator++,distance(face_sequence_vector.begin(), itf),+1));
                }
                else if (FES.coeff(e, *itf) == -1) {
                    edge_fstate_vector.at(e) = -1;
                    sFE_triplet_list.push_back(Tr(edge_numerator++,distance(face_sequence_vector.begin(), itf),-1));
                }
            } // end of  if (FES.coeff(e, fseq) != 0)

    } // end for (auto  itf = face_sequence_vector.begin(); itf != face_sequence_vector.end(); ++itf)

    if (dim == 3) { // nodes - only for the 3D case
        for (auto  ite = edge_fsequence_vector.begin(); ite != edge_fsequence_vector.end(); ++ite) {

            for (unsigned int n = 0; n < CellNumbs.at(0); ++n) //  Loop over all the Nodes in the PCC
                if (ENS.coeff(n, *ite) != 0) {
                    node_fsequence_vector.push_back(n);

                    if (edge_fstate_vector.at(distance(edge_fsequence_vector.begin(), ite)) == 1) {
                        node_fstate_vector.at(n) = 1;
                        sEN_triplet_list.push_back(
                                Tr(node_numerator++, distance(edge_fsequence_vector.begin(), ite), +1));
                    } else if (edge_fstate_vector.at(distance(edge_fsequence_vector.begin(), ite)) == -1) {
                        node_fstate_vector.at(n) = -1;
                        sEN_triplet_list.push_back(
                                Tr(node_numerator++, distance(edge_fsequence_vector.begin(), ite), -1));
                    }
                } // end of if (ENS.coeff(n, *ite) != 0)
        } // end of  for (auto  ite = edge_fsequence_vector.begin(); ite != edge_fsequence_vector.end(); ++ite)
    } // end of if (dim == 3)

/// Forming sparse Incidence matrices: B0, B1 and B2
    ENM_SpecialCells.setFromTriplets(sEN_triplet_list.begin(), sEN_triplet_list.end());
    FEM_SpecialCells.setFromTriplets(sFE_triplet_list.begin(), sFE_triplet_list.end());
    GFM_SpecialCells.setFromTriplets(sGF_triplet_list.begin(), sGF_triplet_list.end());

    /// Identity matrices
    SpMat Id_NSize(ENM_SpecialCells.rows(), ENM_SpecialCells.rows()), Id_ESize(ENM_SpecialCells.cols(),ENM_SpecialCells.cols()), Id_FSize(FEM_SpecialCells.cols(), FEM_SpecialCells.cols());
        for (auto i = 0; i < ENM_SpecialCells.rows(); ++i)
            for (auto j = 0; j < ENM_SpecialCells.rows(); ++j)
                if (i == j) Id_NSize.coeffRef(i,j) = 1;
                    else Id_NSize.coeffRef(i,j) = 0;
        for (auto i = 0; i < ENM_SpecialCells.cols(); ++i)
            for (auto j = 0; j < ENM_SpecialCells.cols(); ++j)
                if (i == j) Id_ESize.coeffRef(i,j) = 1;
                    else Id_ESize.coeffRef(i,j) = 0;
        for (auto i = 0; i < FEM_SpecialCells.cols(); ++i)
            for (auto j = 0; j < FEM_SpecialCells.cols(); ++j)
                if (i == j) Id_FSize.coeffRef(i,j) = 1;
                    else Id_FSize.coeffRef(i,j) = 0;

    L0_laplacian = ENM_SpecialCells * ENM_SpecialCells.transpose();
    L1_laplacian = ENM_SpecialCells.transpose()*Id_NSize*ENM_SpecialCells*Id_ESize + Id_ESize*FEM_SpecialCells*Id_FSize*FEM_SpecialCells.transpose();
    L2_laplacian = FEM_SpecialCells.transpose()*Id_ESize*FEM_SpecialCells*Id_FSize;

    laplacians_vector = { L0_laplacian, L1_laplacian, L2_laplacian };

/// Adjacency matrices A0, A1, A2, A3
//    SpAM_SpecFaces = 0.5 * (SpAM_SpecFaces + SparseMatrix<double>(SpAM_SpecFaces.transpose())); // Full matrix instead of triagonal

    return laplacians_vector;
} // end of std::vector<Eigen::SparseMatrix<double>> ReducedLaplacians( )

/// # 1 # All eigenvalues (spectrum) of a Linear Operator OSM
std::vector<double> OperatorSpectrum(Eigen::SparseMatrix<double> const &OSM, std::vector<unsigned int> &CellNumbs) { // OSM - operator's sparse matrix
std::vector<double> operator_spectrum; // function's output

Eigen::VectorXcd eigval;
//Eigen::MatrixXcd eigvec;

/// Eigenvalues of the S-Faces laplacian
// Construct matrix operation object using the wrapper class SparseGenMatProd
SparseGenMatProd<double> op(OSM);
// Construct eigen solver object, requesting all the eigenvalues
GenEigsSolver<SparseGenMatProd<double>> eigs(op, OSM.rows() - 2, OSM.rows());
//nev	Number of eigenvalues requested. This should satisfy [1 ≤ 𝑛𝑒𝑣 ≤ 𝑛−2], where 𝑛 is the size of matrix.
//ncv	Parameter that controls the convergence speed of the algorithm. Typically a larger ncv means faster convergence, but it may also result in greater memory use and more matrix operations in each iteration. This parameter must satisfy (𝑛𝑒𝑣 + 2 ≤ 𝑛𝑐𝑣 ≤ 𝑛), and is advised to take 𝑛𝑐𝑣 ≥2⋅(𝑛𝑒𝑣 + 1).
// Initialize and compute
eigs.init();
eigs.compute(SortRule::LargestMagn);

// Retrieve results
if (eigs.info() == CompInfo::Successful) {
eigval = eigs.eigenvalues();
// eigvec = eigs.eigenvectors();
}

// new eigs vector
for (unsigned int i = 0; i < OSM.rows() - 2; ++i)
    if(real(eigval(i,0)) >= pow(10,-5))
        operator_spectrum.push_back(real(eigval(i,0)));

for (unsigned int j = 0; j < OSM.rows() - 2; ++j)
    if(real(eigval(j,0)) < pow(10,-5))
        operator_spectrum.push_back(0.0);

return operator_spectrum;
} //END of OperatorSpectrum()

/// # 2 # All eigenvectors of a Linear Operator OSM
Eigen::MatrixXcd OperatorEigvectors(Eigen::SparseMatrix<double> const &OSM, std::vector<unsigned int> &CellNumbs) { // OSM - operator's sparse matrix
//std::vector<double> operator_spectrum; // function's output

Eigen::MatrixXcd eigvec;

/// Eigenvalues of the S-Faces laplacian
// Construct matrix operation object using the wrapper class SparseGenMatProd
SparseGenMatProd<double> op(OSM);
// Construct eigen solver object, requesting all the eigenvalues
GenEigsSolver<SparseGenMatProd<double>> eigs(op, OSM.rows() - 2, OSM.rows());
//nev	Number of eigenvalues requested. This should satisfy [1 ≤ 𝑛𝑒𝑣 ≤ 𝑛−2], where 𝑛 is the size of matrix.
//ncv	Parameter that controls the convergence speed of the algorithm. Typically a larger ncv means faster convergence, but it may also result in greater memory use and more matrix operations in each iteration. This parameter must satisfy (𝑛𝑒𝑣 + 2 ≤ 𝑛𝑐𝑣 ≤ 𝑛), and is advised to take 𝑛𝑐𝑣 ≥2⋅(𝑛𝑒𝑣 + 1).
// Initialize and compute
eigs.init();
eigs.compute(SortRule::LargestMagn);

// Retrieve results
if (eigs.info() == CompInfo::Successful) {
eigvec = eigs.eigenvectors();
}

return eigvec;
} //END of OperatorSpectrum()