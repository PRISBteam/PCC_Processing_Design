/// DCC_Processing
string P_type_here = "R"s;
special_face_design = DCC_Processing(special_faces_sequence, State_sVector, P_type_here);

// special_faces_sequence cut
double special_faces_fraction = 0.5;
vector<unsigned int> sfaces_cut_sequence;

#pragma omp parallel for // parallel execution by OpenMP
for (auto sface_counter = 0; sface_counter < special_faces_fraction*(double) CellNumbs.at(2); ++sface_counter)
sfaces_cut_sequence.push_back(special_faces_sequence.at(sface_counter));

double JEntropy = 0;
//    cout << "R type is here now" << endl;

JEntropy = Get_TJsEntropy(sfaces_cut_sequence);
cout << "TJs configuration entropy: " << JEntropy << endl;

/// DCC_Writer
vector<unsigned int> kinetic_faces_sequence;
int w_counter = 0;
DCC_sequences_Writer(special_faces_sequence, kinetic_faces_sequence, w_counter);
