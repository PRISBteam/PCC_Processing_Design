///================================ DCC Face Laplacian module =============================================================///
///=======================================================================================================================///
/** This subfile calculates the combinatorial Laplacian of the special Faces graph with its spectrum **/
///=======================================================================================================================///

Laplacian_specialFaces = FES * FES.transpose();
EigenSolver<SpMat> Laplacian_Eigenvalues(Laplacian_specialFaces);
VectorXcd Face_eigenValues = Laplacian_Eigenvalues.eigenvalues();
MatrixXcd Face_eigenVectors = Laplacian_Eigenvalues.eigenvectors();
for (int k : Face_eigenValues) cout << k << endl;