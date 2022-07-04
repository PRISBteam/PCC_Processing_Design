///================================ DCC Face Laplacian module =============================================================///
///=======================================================================================================================///
/** This subfile calculates the combinatorial Laplacian of the special Faces graph with its spectrum **/
///=======================================================================================================================///

/*
 * foo() {
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
* return ;
 */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * 	//Eigenvalues for DENSE matrices
	MatrixXd LapFEig(Facenumb, Facenumb);
//	MatrixXd LapFEig(LimFaceNumb, LimFaceNumb);

	for (int i = 0; i < Facenumb; i++) for (int j = 0; j < Facenumb; j++) LapFEig(i,j)=LapF[i][j];
	for (int i = 0; i < LimFaceNumb; i++) for (int j = 0; j < LimFaceNumb; j++) LapFEig(i,j)=LapF[i][j];
	EigenSolver<MatrixXd> es(LapFEig);
	cout<<"es(LapFEig) just calculated"<<"    size = "<< es.eigenvalues().size() <<endl;
	MatrixXcd EigF = es.eigenvalues();

	//REsistivity and Conductivity
	RG=0.0;
//	for (int ig = 0; ig < (es.eigenvalues().size()-1); ig++) cout<< real(EigF(ig,0))<<endl;
	for (int ig = 0; ig < (es.eigenvalues().size()-1); ig++) { RG += pow(real(EigF(ig,0)),-1);}
		RG = RG*es.eigenvalues().size();
		SG = pow(RG,-1);

//	ofstream AMFStream;		//	AMFStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)LapF_HAGBs.txt", ios::trunc);
			AMFStream.open("C:\\Users\\v94623eb\\Dropbox\\Projects\\Current simmulation\\Voronois\\Voronois\\VAnRes\\(HAGBs)Conduct.txt", ios::app);
	//		for (int i = 0; i < Facenumb; i++) { for (int j = 0; j < Facenumb; j++) {
	//				if(AMF[i][j]>=0) AMFStream << AMF[i][j] << "\t"; else j++; }			AMFStream << "\n";		}
//			for (int i = 0; i < Facenumb; i++) {for (int j = 0; j < Facenumb; j++) AMFStream << LapF[i][j] << "\t"; AMFStream << "\n"; }
//			for (int i = 0; i < Facenumb; i++) {for (int j = 0; j < Facenumb; j++) AMFStream << LapFEig(i,j) << "\t"; AMFStream << "\n"; }
//			AMFStream << endl << EigF << endl << endl;
			AMFStream << HAGBsFunc[i][1] << "\t" << RG<< "\t" << SG << endl;
			AMFStream.close();	//for(int kii = 0; kii < Facenumb; kii++)	{for(int l = 0; l < 15; l++) if (NewFaceNeighbours[kii][l]>=0) cout<<"   "<<NewFaceNeighbours[kii][l];} cout<<endl;
//	cout<<"EigF file is written"<<endl;  //	system("pause");
 */