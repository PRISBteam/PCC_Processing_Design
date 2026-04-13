/// FACE barycentre coordinates output to file 'new_face_barycentre_coordinates.txt'
/// Dr Elijah Borodin, 2024

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <set>
#include <vector>

using namespace std;
extern string output_dir;

#include "../../src/lib/PCC_Objects.h"

#include "calculation_face_edge_barycentre_coordinates.h"

void calculate_face_barycentre_coordinates(void) {
    ofstream Out_face_barycentres_file; // Agglomeration sequences output
    string fbc_odir = output_dir + "new_face_barycentre_coordinates.txt"s; // output directory
    Out_face_barycentres_file.open(fbc_odir, ios::trunc);

    PCC current_PCC;
    current_PCC.Set_face_barycentre_coordinates();

    for (auto fbc: current_PCC.Get_face_barycentre_coordinates()) {
// REPAIR cout << get<0>(fbc) << "\t" << get<1>(fbc) << "\t" << get<2>(fbc) << endl;
        Out_face_barycentres_file << get<0>(fbc) << "\t" << get<1>(fbc) << "\t" << get<2>(fbc) << endl;
    }
    Out_face_barycentres_file.close();
}
//exit(0);
void calculate_edge_barycentre_coordinates(void) {
    ofstream Out_edge_barycentres_file; // Agglomeration sequences output
    string ebc_odir = output_dir + "new_edge_barycentre_coordinates.txt"s; // output directory
    Out_edge_barycentres_file.open(ebc_odir, ios::trunc);

    PCC current_PCC;
    current_PCC.Set_edge_barycentre_coordinates();

    for (auto ebc: current_PCC.Get_edge_barycentre_coordinates()) {
// REPAIR cout << get<0>(fbc) << "\t" << get<1>(fbc) << "\t" << get<2>(fbc) << endl;
        Out_edge_barycentres_file << get<0>(ebc) << "\t" << get<1>(ebc) << "\t" << get<2>(ebc) << endl;
    }
    Out_edge_barycentres_file.close();
}

/*
     #include "tasks/calculation_face_edge_barycentre_coordinates.h"
        calculate_edge_barycentre_coordinates();
        calculate_face_barycentre_coordinates();

 */