
void DCC_Plane_cut (std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> const   &c_faces_sequence) {

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file paths.at(5 + (dim - 3))
    SpMat ENS = SMatrixReader(paths.at(4 + (dim - 3)), (CellNumbs.at(0)), (CellNumbs.at(1))); //all Nodes-Edges
    //SpMat AES = SMatrixReader(paths.at(1 + (dim - 3)), (CellNumbs.at(1)), (CellNumbs.at(1))); //all Edges
    ///  Full symmetric AES matrix instead of triagonal
    //AES = 0.5 * (AES + SparseMatrix<double>(AES.transpose()));

    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces
    //SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(2))); //all Faces
    ///  Full symmetric AFS matrix instead of triagonal
    //AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose()));

    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose()));

    /// A new vector<vector<int>> star3vector: grain number (g#) -> grain's nodes numbers (gn#) ::
    vector<vector<int>> star3list; vector<int> nodes_list;
    star3list.clear();
    /// GFS -> FES -> ENS
     for(unsigned int m = 0; m < CellNumbs.at(3); m++) {// each Grain
         nodes_list.clear();
         for(unsigned int l = 0; l < CellNumbs.at(2); l++) {// over all Faces (l)
              if (GFS.coeff(m, l) == 1) {
                      for(unsigned int j = 0; j < CellNumbs.at(1); j++) // over all Edges (j)
                          if (FES.coeff(l, j) == 1) { // at the chosen Face with ID = 'l'
                              for(unsigned int i = 0; i < CellNumbs.at(0); i++) // over all Nodes
                                  if (ENS.coeff(j, i) == 1) nodes_list.push_back(i); // at the chosen Face with ID = 'l'
                          } // end of if (FES.coeff(l, j) == 1)
              } // end of (GFS.coeff(m, l) == 1)

         } // end of for(unsigned int l = 0; l < CellNumbs.at(2); l++) - Faces
    star3list.push_back(nodes_list);
     } // end of for(unsigned int m = 0; m < CellNumbs.at(3); m++) - Grains

    /// Read node coordinates from file

    /// ? new class Grains

    return;
}