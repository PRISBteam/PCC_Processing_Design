///================================ PCC Laplacian laboratory =============================================================///
///=======================================================================================================================///
/** This sub file contain functions allowing to calculate the combinatorial Laplacian of the special k-cells graph with its spectrum and eigenvectors **/
///=======================================================================================================================///

/// # 0 # reduced laplacians (vector of their sparse matrices)
std::vector<Eigen::SparseMatrix<double>> ReducedFaceLaplacians(std::vector<unsigned int> &face_sequence_vector, SpMat &ENS, SpMat &FES, SpMat &GFS) {
std::vector<Eigen::SparseMatrix<double>> laplacians_vector; // function output: { ENM_laplacian, FEM_laplacian, GFM_laplacian }

/// Triplet lists
    vector<Tr> sEN_triplet_list, sFE_triplet_list, sGF_triplet_list;

/// Special face-related PCC elements
    std::vector<unsigned int> node_fsequence_vector, edge_fsequence_vector, grain_fsequence_vector; // number of different sequences based on the face design vector
    std::vector<int> node_ori_vector(CellNumbs.at(0),0), edge_ori_vector(CellNumbs.at(1 + (dim - 3)),0),
                            face_ori_vector(CellNumbs.at(2 + (dim - 3)),0), grain_ori_vector(CellNumbs.at(3 + (dim - 3)),0); // number of different sequences based on the face design vector

//    unsigned int node_numerator = 0, edge_numerator = 0, grain_numerator = 0;
    /// faces loop
    for (auto  itf = face_sequence_vector.begin(); itf != face_sequence_vector.end(); ++itf) {
        // grains
        for (unsigned int g = 0; g < CellNumbs.at(3 + (dim - 3)); ++g) //  Loop over all the Grains in the PCC
                if (GFS.coeff(*itf, g) != 0) {
                    std::vector<unsigned int>::iterator pointer_to_grain = std::find(grain_fsequence_vector.begin(), grain_fsequence_vector.end(),g);

//                    if(grain_fsequence_vector.size() == 0)
//                        grain_fsequence_vector.push_back(g);

                        if(pointer_to_grain == grain_fsequence_vector.end()) {

                        if (GFS.coeff(*itf, g) == 1) {
                            grain_ori_vector.at(g) = 1;
                            sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf), distance(grain_fsequence_vector.begin(), grain_fsequence_vector.end()), +1));
                        } else if (GFS.coeff(*itf, g) == -1) {
                            grain_ori_vector.at(g) = -1;
                            sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf), distance(grain_fsequence_vector.begin(), grain_fsequence_vector.end()), -1));
                        }
                        /// push back
                            grain_fsequence_vector.push_back(g);
                        } // end if(pointer_to_grain == grain_fsequence_vector.end())
                    else if(pointer_to_grain != grain_fsequence_vector.end()) {
                        if (GFS.coeff(*itf, g) == 1) {
                            grain_ori_vector.at(g) = 1;
                            sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf), distance(grain_fsequence_vector.begin(), pointer_to_grain), +1));
                        } else if (GFS.coeff(*itf, g) == -1) {
                            grain_ori_vector.at(g) = -1;
                            sGF_triplet_list.push_back(Tr(distance(face_sequence_vector.begin(), itf), distance(grain_fsequence_vector.begin(), pointer_to_grain), -1));
                        }
                    } // else if(pointer_to_grain != grain_fsequence_vector.end()) {
                } // end if (GFS.coeff(*itf, g) == 1)

        // edges
        for (unsigned int e = 0; e < CellNumbs.at(1 + (dim - 3)); ++e) //  Loop over all the Edges in the PCC
                if (FES.coeff(e, *itf) != 0) {
                    std::vector<unsigned int>::iterator pointer_to_edge = std::find(edge_fsequence_vector.begin(), edge_fsequence_vector.end(),e);

//                    if(edge_fsequence_vector.size() == 0)
//                        edge_fsequence_vector.push_back(e);

                    if(pointer_to_edge == edge_fsequence_vector.end()) {

                    if (FES.coeff(e, *itf) == 1) {
                        edge_ori_vector.at(e) = 1;
                        sFE_triplet_list.push_back(Tr(distance(edge_fsequence_vector.begin(), edge_fsequence_vector.end()), distance(face_sequence_vector.begin(), itf),1));
                    } else if (FES.coeff(e, *itf) == -1) {
                        edge_ori_vector.at(e) = -1;
                        sFE_triplet_list.push_back(Tr(distance(edge_fsequence_vector.begin(), edge_fsequence_vector.end()), distance(face_sequence_vector.begin(), itf),-1));
                    }
                    /// push back
                        edge_fsequence_vector.push_back(e);
                    } // end of if(pointer_to_edge == edge_fsequence_vector.end())
                else if(pointer_to_edge != edge_fsequence_vector.end()) {
                    if (FES.coeff(e, *itf) == 1) {
                        edge_ori_vector.at(e) = 1;
                        sFE_triplet_list.push_back(Tr(distance(edge_fsequence_vector.begin(), pointer_to_edge), distance(face_sequence_vector.begin(), itf),1));
                    } else if (FES.coeff(e, *itf) == -1) {
                        edge_ori_vector.at(e) = -1;
                        sFE_triplet_list.push_back(Tr(distance(edge_fsequence_vector.begin(), pointer_to_edge), distance(face_sequence_vector.begin(), itf),-1));
                    }
                } // end of if(pointer_to_edge == edge_fsequence_vector.end())

                } // end of  if (FES.coeff(e, fseq) != 0)

    } // end for (auto  itf = face_sequence_vector.begin(); itf != face_sequence_vector.end(); ++itf)

// REPAIR    for (auto sge : sGF_triplet_list)
//        cout << "sfe: " << sge.row() << "  " << sge.col() << "  "  << sge.value() << endl;
// cout << endl;
// REPAIR    for (auto sfe : sFE_triplet_list)
//        cout << "sfe: " << sfe.row() << "  " << sfe.col() << "  "  << sfe.value() << endl;

    /// Sparse matrices (Eigen/ Spectra)
    SpMat FEM_SpecialCells(edge_fsequence_vector.size(), face_sequence_vector.size()),
            GFM_SpecialCells(face_sequence_vector.size(),grain_fsequence_vector.size());

//    cout << "F_G" << endl;
//    for (auto gtl : sGF_triplet_list)
//        cout << gtl.col() << " " << gtl.row() << " " << gtl.value() << endl;
//cout << endl;
//    cout << "E_F" << endl;
//    for (auto ftl : sFE_triplet_list)
//        cout << ftl.col() << " " << ftl.row() << " " << ftl.value() << endl;
//    cout << endl;
//    cout << edge_fsequence_vector.size() << " " << face_sequence_vector.size() << " " << grain_fsequence_vector.size() << endl;

//    for (auto fss : face_sequence_vector)
//        cout << fss << " ";
//        cout << endl;

    /// Forming sparse Incidence matrices: B1 and B2
    GFM_SpecialCells.setFromTriplets(sGF_triplet_list.begin(), sGF_triplet_list.end());
    FEM_SpecialCells.setFromTriplets(sFE_triplet_list.begin(), sFE_triplet_list.end());

    if (dim == 3) { // nodes - only for the 3D case
    /// faces loop
        for (auto  ite = edge_fsequence_vector.begin(); ite != edge_fsequence_vector.end(); ++ite) {

            for (unsigned int n = 0; n < CellNumbs.at(0); ++n) {//  Loop over all the Nodes in the PCC
//                cout << "HERE3:  " << *ite << "  " << n << endl;

                if (ENS.coeff(n, *ite) != 0) {
                    std::vector<unsigned int>::iterator pointer_to_node = std::find(node_fsequence_vector.begin(), node_fsequence_vector.end(),n);

//                    if(node_fsequence_vector.size() == 0)
//                        node_fsequence_vector.push_back(n);

                    if(pointer_to_node == node_fsequence_vector.end()) {
                        if (ENS.coeff(n, *ite) == 1) {
//                        cout << "HERE33:  " << *ite << "  " << n << endl;
                            node_ori_vector.at(n) = 1;
                            //                       cout << "HERE34:  " << *ite << "  " << n << endl;
                            sEN_triplet_list.push_back(Tr(distance(node_fsequence_vector.begin(), node_fsequence_vector.end()), distance(edge_fsequence_vector.begin(), ite),1));

                        } else if (ENS.coeff(n, *ite) == -1) {
                            node_ori_vector.at(n) = -1;
                            sEN_triplet_list.push_back(Tr(distance(node_fsequence_vector.begin(), node_fsequence_vector.end()), distance(edge_fsequence_vector.begin(), ite),-1));
                            //                        cout << "HERE5:  " << *ite << "  " << n << endl;
                        }
                    /// push back
                        node_fsequence_vector.push_back(n);
                    } // end of if(pointer_to_node == edge_fsequence_vector.end())
                    else if(pointer_to_node != node_fsequence_vector.end()) {
                        if (ENS.coeff(n, *ite) == 1) {
                            node_ori_vector.at(n) = 1;
//                       cout << "HERE34:  " << *ite << "  " << n << endl;
                            sEN_triplet_list.push_back(Tr(distance(node_fsequence_vector.begin(), pointer_to_node), distance(edge_fsequence_vector.begin(), ite),1));
                        } else if (ENS.coeff(n, *ite) == -1) {
                            node_ori_vector.at(n) = -1;
                            sEN_triplet_list.push_back(Tr(distance(node_fsequence_vector.begin(), pointer_to_node), distance(edge_fsequence_vector.begin(), ite),-1));
//                        cout << "HERE5:  " << *ite << "  " << n << endl;
                        }
                    } // end of else if(pointer_to_node != node_fsequence_vector.end())

                } // end of if (ENS.coeff(n, *ite) != 0)
            } // end of for (unsigned int n = 0; n < CellNumbs.at(0); ++n)

        } // end of  for (auto  ite = edge_fsequence_vector.begin(); ite != edge_fsequence_vector.end(); ++ite)
    } // end of if (dim == 3)
    SpMat ENM_SpecialCells(node_fsequence_vector.size(),edge_fsequence_vector.size());

    /// Forming sparse Incidence matrix: B0
    ENM_SpecialCells.setFromTriplets(sEN_triplet_list.begin(), sEN_triplet_list.end());

    /// Identity matrices
    std::vector<Tr> NSize_tr_list, ESize_tr_list, FSize_tr_list;
    SpMat Id_NSize(ENM_SpecialCells.rows(), ENM_SpecialCells.rows()), Id_ESize(ENM_SpecialCells.cols(),ENM_SpecialCells.cols()), Id_FSize(FEM_SpecialCells.cols(), FEM_SpecialCells.cols());
//        cout << "HERE1" << endl;

        for (auto i = 0; i < ENM_SpecialCells.rows(); ++i)
            for (auto j = 0; j < ENM_SpecialCells.rows(); ++j)
                if (i == j) NSize_tr_list.push_back(Tr(i, j, 1));
                else NSize_tr_list.push_back(Tr(i, j, 0));

        Id_NSize.setFromTriplets(NSize_tr_list.begin(), NSize_tr_list.end());
//    cout << "HERE2" << endl;

        for (auto i = 0; i < ENM_SpecialCells.cols(); ++i)
            for (auto j = 0; j < ENM_SpecialCells.cols(); ++j)
                if (i == j) ESize_tr_list.push_back(Tr(i, j, 1));
                    else ESize_tr_list.push_back(Tr(i, j, 0));

        Id_ESize.setFromTriplets(ESize_tr_list.begin(), ESize_tr_list.end());
//    cout << "HERE3" << endl;

        for (auto i = 0; i < FEM_SpecialCells.cols(); ++i)
            for (auto j = 0; j < FEM_SpecialCells.cols(); ++j)
                if (i == j) FSize_tr_list.push_back(Tr(i, j, 1));
                    else FSize_tr_list.push_back(Tr(i, j, 0));

        Id_FSize.setFromTriplets(FSize_tr_list.begin(), FSize_tr_list.end());
//    cout << "HERE4" << endl;

//    cout << "HERE and Now 9 !  " << FEM_SpecialCells.cols() << "  " << ENM_SpecialCells.cols() << "  " << ENM_SpecialCells.rows() << endl;
//    cout << "HERE and Now 9 !  " << Id_FSize.rows() << "  " << Id_ESize.cols() << "  " << Id_NSize.rows() << endl;
//    cout << "HERE and Now 9 !  " << Id_FSize.nonZeros() << "  " << Id_ESize.nonZeros() << "  " << Id_NSize.nonZeros() << endl;

    SpMat L0_laplacian, L1_laplacian, L2_laplacian;
//    cout << "HERE5" << endl;

    L0_laplacian = ENM_SpecialCells * ENM_SpecialCells.transpose();
//    cout << "HERE51" << endl;
    SpMat L11_laplacian = ENM_SpecialCells.transpose()*Id_NSize*ENM_SpecialCells*Id_ESize;
//    cout << "HERE52" << endl;
    SpMat L12_laplacian = Id_ESize*FEM_SpecialCells*Id_FSize*FEM_SpecialCells.transpose();
//    cout << "HERE53" << endl;
    L1_laplacian = L11_laplacian + L12_laplacian;
//    cout << "HERE4" << endl;
    L2_laplacian = FEM_SpecialCells.transpose()*Id_ESize*FEM_SpecialCells*Id_FSize;
//    cout << "HERE6" << endl;

    laplacians_vector = { L0_laplacian, L1_laplacian, L2_laplacian };

    return laplacians_vector;
} // end of std::vector<Eigen::SparseMatrix<double>> ReducedLaplacians( )

/// # 1 # All eigenvalues (spectrum) of a Linear Operator OSM
std::vector<double> OperatorSpectrum(Eigen::SparseMatrix<double> const &OSM) { // OSM - operator's sparse matrix
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
Eigen::MatrixXcd OperatorEigvectors(Eigen::SparseMatrix<double> const &OSM) { // OSM - operator's sparse matrix
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

/// # 3 # Betti numbers

std::vector<double> OperatorsBetti(std::vector<unsigned int> &face_sequence_vector, SpMat &ENS, SpMat &FES, SpMat &GFS) {
    std::vector<double> Betti_vector;
    std::vector<vector<double>> Laplace_spectra;

    std::vector<SpMat> Laplacians_vector; // [0] -> L0, [1] -> L1, [2] -> L2
    cout << "[OperatorsBetti] Start of the creation Reduced Face Laplacians: " << endl;
    Laplacians_vector = ReducedFaceLaplacians(face_sequence_vector, ENS, FES, GFS);

    cout << "[OperatorsBetti] face sequence vector size: " << face_sequence_vector.size() << endl;
    cout << "[OperatorsBetti] Reduced Laplacians row sizes: " << Laplacians_vector.at(0).rows() << "  " << Laplacians_vector.at(0).cols() << "  " << Laplacians_vector.at(1).rows() << "  " << Laplacians_vector.at(1).cols() << "  " << Laplacians_vector.at(2).rows() << "  " << Laplacians_vector.at(2).cols() << endl;
    Out_logfile_stream << "[OperatorsBetti] face sequence vector size: " << face_sequence_vector.size() << endl;
    Out_logfile_stream << "[OperatorsBetti] Reduced Laplacians row sizes: " << Laplacians_vector.at(0).rows() << "  " << Laplacians_vector.at(0).cols() << "  " << Laplacians_vector.at(1).rows() << "  " << Laplacians_vector.at(1).cols() << "  " << Laplacians_vector.at(2).rows() << "  " << Laplacians_vector.at(2).cols() << endl;

    cout << "[OperatorsBetti] Start of the Laplacians' spectra calculation: " << endl;
    Out_logfile_stream << "[OperatorsBetti] Start of the Laplacians' spectra calculation: " << endl;

#pragma omp parallel for // parallel execution by OpenMP
    for (auto LVS : Laplacians_vector) {
        Laplace_spectra.push_back(OperatorSpectrum(LVS));
        cout << "[OperatorsBetti] Reduced Laplacians' spectra size: " << Laplace_spectra.back().size() << "  ";
        Out_logfile_stream << "[OperatorsBetti] Reduced Laplacians' spectra size: " << Laplace_spectra.back().size() << "  ";
// REPAIR        for(auto sp : Laplace_spectra.back()) //        cout << sp << " ";  cout << endl;
    }

    cout << endl << "[OperatorsBetti] Start of the Betti numbers calculation: " << endl;
    Out_logfile_stream << endl << "[OperatorsBetti] Start of the Betti numbers calculation: " << endl;
    int Betticounter = 0;
    for (auto ls : Laplace_spectra) { // normally only 3 vectors
        Betti_vector.push_back(std::count(ls.begin(), ls.end(), 0));
        cout << "[OperatorsBetti] Betti " << Betticounter++ << " is equal to " << Betti_vector.back() << endl;
        Out_logfile_stream << "[OperatorsBetti] Betti " << Betticounter << " is equal to " << Betti_vector.back() << endl;
    }
    return Betti_vector;
} // END of OperatorsBetti()

double gb_Resistivity(std::vector<int> &face_states_vector, Eigen::SparseMatrix<double> const &WLMatrix){
    double resistivity = 0;
    std::vector<double> operator_spectrum = OperatorSpectrum(WLMatrix);

//    for (double os : operator_spectrum)
        for (unsigned int  os = 0; os < operator_spectrum.size() - 2; ++os) {
//            cout << operator_spectrum.at(os) << endl;
            if (operator_spectrum.at(os) != 0.0) resistivity += 1.0 / operator_spectrum.at(os);
//            else resistivity = 1000000.0; // = infinity
        }
    return resistivity;
} // end of function gb_Resistivity()