///================================ A part of the DCC Writer module ====================================================================///
///=====================================================================================================================================///
/** The library contains functions for the formatted output of various kind of data generated by other modules                         **/
///===================================================================================================================================///
//Output directory
string odir = output_folder;
//char* odir = const_cast<char*>(output_folder.c_str());

int DCC_sequences_Writer(std::vector<unsigned int> const  &special_faces_sequence, std::vector<unsigned int> const &ordinary_faces_sequence) {
/// output special_Face_sequence to file
    ofstream OutSfaces_file, OutOfaces_file; // Special and Ordinary faces sequence output
/// Output to file Special and Ordinary faces order :: tess - means "numeration of Faces start with 1 instead of 0 like in the NEPER output"
    // Special
    OutSfaces_file.open(output_folder + "SpecialFaces.txt"s, ios::trunc); OutOfaces_file.open(output_folder + "OrdinaryFaces.txt"s, ios::trunc);
    if (OutSfaces_file) { //        OutSGBfile << "Global numbers (in DCC) of special grain boundaries with the fraction " << special_faces_sequence.size()/ CellNumbs.at(2) << endl;
        for (auto vit: special_faces_sequence) OutSfaces_file << vit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
        cout << "(1) Special faces_sequence has been successfully written in " << output_folder + "SpecialFaces.txt"s << endl;
        OutSfaces_file.close();
    } else cout << "Error: No such a directory for\t" << output_folder << "SpecialFaces.txt"s << endl;
    // Ordinary
    if (OutOfaces_file) {
    for (auto vit: ordinary_faces_sequence) OutOfaces_file << vit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
        cout << "(2) Special faces_sequence has been successfully written in " << output_folder << "OrdinaryFaces.txt"s << endl;
    OutOfaces_file.close();
    } else cout << "Error: No such a directory for\t" << output_folder << "OrdinaryFaces.txt"s << endl;

    return 0;
}