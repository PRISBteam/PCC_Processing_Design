/// Attached user defined C++ libraries:
///-------------------------------------
#include "Planecut_Functions.h"
///-------------------------------------

std::vector <unsigned int> DCC_Subcomplex(subcomplex new_cut, std::vector<unsigned int> const &s_faces_sequence, std::vector<unsigned int> const &c_faces_sequence) {
// sub_grains_sequence - all grains in the subcomplex, sub_faces_sequence - all faces in the subcomplex, common_faces_sequence - all faces common for two grains in the subcomplex, s_sub_faces_sequence - special faces, c_sub_faces_sequence - induced (fractured, for instance) faces
std::vector <unsigned int>  sub_grains_sequence, sub_faces_sequence, common_faces_sequence, s_sub_faces_sequence, c_sub_faces_sequence;

    SpMat GFS = SMatrixReader(paths.at(6 + (dim - 3)), (CellNumbs.at(2)), (CellNumbs.at(3))); //all Faces-Grains
    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Volumes
    ///  Full symmetric AGS matrix instead of triagonal
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose()));

    /// Vertex coordinates reader into triplet double vector
    string VCpath_string = input_folder + "voro_seeds.txt"s;
    string GCpath_string = input_folder + "grain_seeds.txt"s;
    char* VCpath = const_cast<char*>(VCpath_string.c_str());
    char* GCpath = const_cast<char*>(GCpath_string.c_str());
    vector<tuple<double, double, double>> vertex_coordinates = TuplesReader(VCpath);
    vector<tuple<double, double, double>> face_coordinates;
    vector<tuple<double, double, double>> grain_coordinates = TuplesReader(VCpath);

/// All subcomplex grains (sub_grains_sequence)
double a_coeff = 0.0, b_coeff = 0.0, c_coeff = 1.0, D_coeff = 0.5;
sub_grains_sequence = DCC_Plane_cut(a_coeff, b_coeff, c_coeff, D_coeff, vertex_coordinates);
int subcomplex_id = 0;
double crack_length_ratio = 0.1;
new_cut.Get_half_plane(subcomplex_id, sub_grains_sequence, crack_length_ratio, 1);


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
            face_coordinates.push_back(find_aGBseed(fnumber, paths, CellNumbs, grain_coordinates)); // tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);
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


return s_sub_faces_sequence;
}

/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence