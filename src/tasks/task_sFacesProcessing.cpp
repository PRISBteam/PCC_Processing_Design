/// DCC_Processing
string P_type_here = "R"s;
cout << "------" << endl; cout << "START of the DCC Processing module" << endl; Out_logfile_stream << "------" << endl; Out_logfile_stream << "START of the DCC Processing module" << endl;
//special_cells_design = DCC_Processing(special_faces_sequence, State_sVector, P_type_here);

vector<unsigned int> kinetic_faces_sequence;
int w_counter = 0;

for(int j=0; j < 1; ++j){
for( int i = 0; i < 1000000; ++i){
unsigned int New2CellNumb = 0;
New2CellNumb = NewCellNumb_R(1000000000); // advanced random generator of  pecial faces
special_faces_sequence.push_back(New2CellNumb);
}

ofstream OutSfaces_file; // Special and Ordinary faces sequence output
OutSfaces_file.open(output_dir + "RandomSpecialFaces.txt"s, ios::trunc);
if (OutSfaces_file) { //        OutSGBfile << "Global numbers (in DCC) of special grain boundaries with the fraction " << special_faces_sequence.size()/ CellNumbs.at(2) << endl;
for (auto rvit: special_faces_sequence) OutSfaces_file << rvit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
OutSfaces_file.close();
} else cout << "Error: No such a directory for\t" << output_dir << "RandomSpecialFaces.txt"s << endl;


/// ===== Elapsing time Processing ================
unsigned int Processing_time = clock();
P_time = (double) Processing_time;

cout << "Processing time is equal to  " << P_time / ((j+1)*1000000.0*pow(10.0, 6.0)) << "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl;
cout << "Processing time is equal to  " << P_time / ((j+1)*pow(10.0, 6.0)) << "  seconds" << endl; cout << "-------------------------------------------------------------------------" << endl;
P_time = 0;
}

exit(0);

Out_logfile_stream << "Processing time is equal to  " << P_time / pow(10.0, 6.0) << "  seconds" << endl; Out_logfile_stream << "-------------------------------------------------------------------------" << endl;

//exit(0);

// special_faces_sequence cut
double special_faces_fraction = 0.99;
vector<unsigned int> sfaces_cut_sequence;

//#pragma omp parallel for // parallel execution by OpenMP
///for (auto sface_counter = 0; sface_counter < special_faces_fraction*(double) CellNumbs.at(2); ++sface_counter)
for (auto sface_counter = 0; sface_counter < special_faces_sequence.size(); ++sface_counter)
    sfaces_cut_sequence.push_back(special_faces_sequence.at(sface_counter));

double JEntropy = 0;
//    cout << "R type is here now" << endl;

JEntropy = Get_TJsEntropy(sfaces_cut_sequence);
cout << "TJs configuration entropy: " << JEntropy << endl;

/// DCC_Writer
cout << "START of the DCC Writer module" << endl; Out_logfile_stream << "START of the DCC Writer module" << endl;
DCC_sequences_Writer(special_faces_sequence, kinetic_faces_sequence, w_counter);

/// ===== Elapsing time Writer ================
unsigned int Writer_time = clock();
W_time = (double) Writer_time - P_time;
cout << "------" << endl; cout << "Writer time is equal to  " << W_time / pow(10.0, 6.0) << "  seconds" << endl; cout << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
Out_logfile_stream << "------" << endl; Out_logfile_stream << "Writer time is equal to  " << W_time / pow(10.0, 6.0) << "  seconds" << endl; Out_logfile_stream << "--------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
