/// Attached user defined C++ libraries:
///-------------------------------------
#include "Planecut_Functions.h"
///-------------------------------------

std::vector <unsigned int> DCC_Subcomplex(double a_coeff, double b_coeff, double c_coeff, double D_coeff, std::vector<unsigned int> const &s_faces_sequence) {
std::vector <unsigned int>  sub_grains_sequence, sub_faces_sequence, common_faces_sequence, s_sub_faces_sequence;

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
sub_grains_sequence = DCC_Plane_cut (a_coeff, b_coeff, c_coeff, D_coeff, vertex_coordinates);

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

//s_faces_sequence

return s_sub_faces_sequence;
}

/// Heap ///
//std::vector<unsigned int> const   &c_faces_sequence