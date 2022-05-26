///================================ DCC Face Laplacian module =============================================================///
///=======================================================================================================================///
/** This subfile calculates the combinatorial Laplacian of the special Faces graph with its spectrum **/
///=======================================================================================================================///

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
