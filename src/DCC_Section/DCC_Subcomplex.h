/// Attached user defined C++ libraries:
///-------------------------------------
#include "Planecut_Functions.h"
///-------------------------------------

//std::vector <unsigned int> DCC_Subcomplex(subcomplex new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> const &c_faces_sequence) {
subcomplex DCC_Subcomplex(subcomplex &new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> const &c_faces_sequence) {
// sub_grains_sequence - all grains in the subcomplex, sub_faces_sequence - all faces in the subcomplex, common_faces_sequence - all faces common for two grains in the subcomplex, s_sub_faces_sequence - special faces, c_sub_faces_sequence - induced (fractured, for instance) faces
std::vector <unsigned int>  sub_grains_sequence, sub_faces_sequence, common_faces_sequence, s_sub_faces_sequence, c_sub_faces_sequence;

    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose()));

    /// Vertex coordinates reader into triplet double vector

    vector<tuple<double, double, double>> vertex_coordinates = vertex_coordinates_vector;
    vector<tuple<double, double, double>> common_faces_coordinates;
    vector<tuple<double, double, double>> subcomplex_face_coordinates;
    vector<tuple<double, double, double>> subcomplex_grain_coordinates;
    vector<tuple<double, double, double>> grain_coordinates = grain_coordinates_vector;

/// All subcomplex grains (subcomplex_grain_sequence)
double a_coeff = 0.0, b_coeff = 0.0, c_coeff = 1.0, D_coeff = 0.5;
/// Grains in a plane
sub_grains_sequence = DCC_Plane_cut(a_coeff, b_coeff, c_coeff, D_coeff);

/// Common grain coordinates
    for(auto subgc : sub_grains_sequence)
        subcomplex_grain_coordinates.push_back(grain_coordinates.at(subgc));

/// All subcomplex faces (sub_faces_sequence)
if (sub_grains_sequence.size() > 0) {
    for (auto grain_id : sub_grains_sequence) {

        for (unsigned int l = 0; l < CellNumbs.at(2); l++)
            if (GFS.coeff(l, grain_id) == 1)
                sub_faces_sequence.push_back(l);
    } // end for (auto grain_id : sub_grains_sequence)
}// end of if(sub_grains_sequence.size() > 0)

/// Common subcomplex faces (common_faces_sequence)
if (sub_faces_sequence.size() > 0)
    for (auto face_id : sub_faces_sequence)
        if(count(sub_faces_sequence.begin(), sub_faces_sequence.end(),face_id) == 2)
            common_faces_sequence.push_back(face_id);

    /// Common face coordinates
    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber) {
        if(std::find(common_faces_sequence.begin(), common_faces_sequence.end(), fnumber) != common_faces_sequence.end())
            common_faces_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates)); // tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);
        //NewSeedsStream << fnumber << "\t" << get<0>(NewSeed_coordinates) << "\t" << get<1>(NewSeed_coordinates) << "\t" << get<2>(NewSeed_coordinates) << endl;
    }
    /// Special faces
    if (s_faces_sequence.size() > 0)
        for (unsigned int face_number : s_faces_sequence)
            if (std::find(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_number) != sub_faces_sequence.end())
                s_sub_faces_sequence.push_back(face_number);

    /// Induced faces
    if (c_faces_sequence.size() > 0)
        for (unsigned int face_number : c_faces_sequence)
            if (std::find(sub_faces_sequence.begin(), sub_faces_sequence.end(), face_number) != sub_faces_sequence.end())
                c_sub_faces_sequence.push_back(face_number);

    /// Setting all quantities to the subcomplex new_subPCC with id = 0
    new_cut.Set_grains_sequence(sub_grains_sequence);
    new_cut.Set_faces_sequence(sub_faces_sequence);
    new_cut.Set_common_faces_coordinates(common_faces_coordinates);
    new_cut.Set_sub_grain_coordinates(subcomplex_grain_coordinates);
    new_cut.Set_sfaces_sequence(s_sub_faces_sequence); //special faces
    new_cut.Set_cfaces_sequence(c_sub_faces_sequence); //cracked (induced) faces

return new_cut;
}

/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence