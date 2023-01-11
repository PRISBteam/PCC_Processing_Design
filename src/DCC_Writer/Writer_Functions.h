///================================ A part of the DCC Writer module ====================================================================///
///=====================================================================================================================================///
/** The library contains functions for the formatted output of various kind of data generated by other modules                         **/
///===================================================================================================================================///
//Output directory
string odir = output_folder;
//char* odir = const_cast<char*>(output_folder.c_str());

int DCC_sequences_Writer(std::vector<unsigned int> const  &special_faces_sequence, std::vector<unsigned int> const &kinetic_faces_sequence, int &output_counter) {
/// output special_Face_sequence to file
    ofstream OutSfaces_file, OutOfaces_file; // Special and Ordinary faces sequence output
/// Output to file Special and Ordinary faces order :: tess - means "numeration of Faces start with 1 instead of 0 like in the NEPER output"
    // Random
    OutSfaces_file.open(output_folder + "RandomSpecialFaces.txt"s, ios::trunc);
    if (OutSfaces_file) { //        OutSGBfile << "Global numbers (in DCC) of special grain boundaries with the fraction " << special_faces_sequence.size()/ CellNumbs.at(2) << endl;
        for (auto rvit: special_face_design.at(0)) OutSfaces_file << rvit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
        cout <<"("<<output_counter++<<")  "<< "Random special faces_sequence has been successfully written in " << output_folder + "RandomSpecialFaces.txt"s << endl;
        OutSfaces_file.close();
    } else cout << "Error: No such a directory for\t" << output_folder << "RandomSpecialFaces.txt"s << endl;

// Smax
    if (design_number > 0) {
        OutOfaces_file.open(output_folder + "LengthySpecialFaces.txt"s, ios::trunc);
        if (OutOfaces_file) {
            for (auto smvit: special_face_design.at(1))
                OutOfaces_file << smvit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
            cout << "(" << output_counter++ << ")  " << "Lengthy random faces_sequence has been successfully written in "
                 << output_folder
                 << "LengthySpecialFaces.txt"s << endl;
            OutOfaces_file.close();
        } else cout << "Error: No such a directory for\t" << output_folder << "SmaxSpecialFaces.txt"s << endl;
    }

    if(kinetic_faces_sequence.size() > 0) {
        /// output cracked_Face_sequence to file
        ofstream OutCFSfile; // Special Faces sequence output
        /// Output to file Cracked Faces order :: tess - means "numeration of Faces start with 1 instead of 0 like in the NEPER output"
        OutCFSfile.open(output_folder + "CrackedGrainBoundaries.txt"s, ios::trunc);
        if (OutCFSfile) {
//        OutSGBfile << "Global numbers (in DCC) of cracked grain boundaries with the fraction " << special_faces_sequence.size()/ CellNumbs.at(2) << endl;
            for (auto vit: kinetic_faces_sequence)
                OutCFSfile << vit + 1 << endl; /// vit + 1 !!! for compatibility with the Neper output
//        { if (unsigned int numerator = 0; numerator < max_cFaces_fraction*CellNumbs.at(2)) OutCFSfile << vit + 1 << endl; ++numerator; }
            cout << "("<<output_counter++<<")  " << "Kinetic generated special faces sequence has been successfully written in " << output_folder
                 << "CrackedGrainBoundaries.txt"s << endl;
            OutCFSfile.close();
        } else cout << "Error: No such a directory for\t" << output_folder << "CrackedGrainBoundaries.txt"s << endl;
    }
    return 0;
}

/// Fro 'L' processing type
//--------------------------------------------------------------
void Strips_distribution_Writer(int &output_counter) { // strip distribution output
    ofstream OutStripsDistribution; // Strip distributions output
    /// Output to file Cracked Faces order :: tess - means "numeration of Faces start with 1 instead of 0 like in the NEPER output"
    OutStripsDistribution.open(output_folder + "P_StripsDistribution.txt"s, ios::trunc);
    if (OutStripsDistribution) {
        for (auto sdis : face_strip_distribution)
            OutStripsDistribution << sdis << endl;

        cout <<"("<<output_counter++<<")  " << " Special face strips distribution has been successfully written in " << output_folder
             << "P_StripsDistribution.txt"s << endl;
        cout << "The total sum of all the values in the distribution is equal to:  " << std::accumulate(face_strip_distribution.begin(), face_strip_distribution.end(),
                                                                                                        decltype(face_strip_distribution)::value_type(0)) << endl;

        OutStripsDistribution.close();
    } else cout << "Error: No such a directory for\t" << output_folder << "P_StripsDistribution.txt"s << endl;

    return;
}

void RWseries_Writer(vector<vector<int>> const &RW_series_vector, int &output_counter) {
    ofstream RWseriesOutput; // Strip distributions output
    /// Output to file Cracked Faces order :: tess - means "numeration of Faces start with 1 instead of 0 like in the NEPER output"
    RWseriesOutput.open(output_folder + "P_RWseries.txt"s, ios::trunc);
    if (RWseriesOutput) {
        for (auto RWsfv : RW_series_vector) {
            for (auto RWsf: RWsfv) {
                RWseriesOutput << RWsf << " ";
            }
            RWseriesOutput << endl;
        }
        cout <<"("<<output_counter++<<")  " << " Strips of special faces has been successfully written in " << output_folder <<
             "P_RWseries.txt"s << endl;

        RWseriesOutput.close();
    } else cout << "Error: No such a directory for\t" << output_folder << "P_RWseries.txt"s << endl;

    return;
}

void Agglomerations_Writer(int &output_counter, vector<vector<int>> const &RW_series_vector, double const &mu_f, double const &sigm_f) {
    /// vector containing agglomerations (DCC_Objects_class objects) as its elements
    vector<agglomeration> agglomerations_vector;
    agglomerations_vector.clear();

    ///Agglomerations: fraction and average value
    vector<unsigned int> a_State_Vector(CellNumbs.at(2),0);

    /// form preliminary State Vector of agglomerations
    for (auto RWsfv : RW_series_vector)
        for (auto RWsf: RWsfv)
            a_State_Vector.at(RWsf) += 1; // change element of the State Vector

    /// form a new Vector of agglomerations
    for (auto  a_vector_itr = a_State_Vector.begin(); a_vector_itr != a_State_Vector.end(); ++a_vector_itr)
    {
        if (*a_vector_itr > 1)
            // agglomeration(unsigned int AFace, int AglPower)
            agglomerations_vector.push_back( agglomeration(distance(a_State_Vector.begin(), a_vector_itr), *a_vector_itr) );

    } // end of for (auto  a_vector_itr = a_State_Vector.begin(); a_vector_itr != a_State_Vector.end(); ++a_vector_itr)
    a_State_Vector.clear();

    /// Agglomerations  output:
    OutAgglStatistics.open(output_folder + "Stat_Agglomerations.txt"s, ios::app);
    OutAvLengthsADistributions.open(output_folder + "AvLengths_Distr_Agglomerations.txt"s, ios::trunc);
    OutPowersADistributions.open(output_folder + "Powers_Distr_Agglomerations.txt"s, ios::trunc);

    if (OutAgglStatistics && OutAvLengthsADistributions && OutPowersADistributions) {
        /// Calculation of an average length of strips related to the agglomeration
        int AgglAvLength = 0, AgglAvPower = 0;

        /// 1. Agglomeration distributions
        for (auto aggl: agglomerations_vector) {
            AgglAvLength += aggl.GetAvLength(RW_series_vector); //int GetAvLength(vector<vector<int>> const &RW_series_vector)
            OutAvLengthsADistributions << aggl.GetAvLength() << endl; //file output

            AgglAvPower += aggl.GetAPower(RW_series_vector); //    int GetAPower(vector<vector<int>> const &RW_series_vector)
            OutPowersADistributions << aggl.GetAPower() << endl; //file output
        } // end of for (auto aggl: agglomerations_vector)

        /// 2. Agglomeration statistics
        AgglAvLength /= agglomerations_vector.size();
        AgglAvPower /= agglomerations_vector.size();

        /// output to the Stat_Agglomerations.txt file
        // (1) agglomerations fraction (2) mu_f (3) sigm_f (4) agglomeration average power (5) average length of strips related to the agglomeration
        OutAgglStatistics << (double) agglomerations_vector.size() / CellNumbs.at(2) << "   " << mu_f << "   " << sigm_f << "   " << AgglAvLength << "   " << AgglAvPower << endl;
        //                OutAgglomerations << (double) special_faces_sequence.size()/ CellNumbs.at(2) << "   " << mu_f << "   " << sigm_f << "   " << (double) 100.0* agglomerations_Vector.size()/ CellNumbs.at(2) << "   " << std::accumulate(agglomerations_Vector.begin(), agglomerations_Vector.end(), decltype(agglomerations_Vector)::value_type(0))/ agglomerations_Vector.size() << endl; //        cout << "(" << output_counter++ << ")  " << " Agglomerations has been successfully written in " << output_folder << "P_Agglomerations.txt"s << endl;

    } else cout << "Error: No such a directory for\t" << output_folder << "Stat_Agglomerations.txt or Powers_Distr_Agglomerations.txt or AvLengths_Distr_Agglomerations"s << endl;

    OutAgglStatistics.close();
    OutAvLengthsADistributions.close();
    OutPowersADistributions.close();

    return;
} // end of void Agglomerations_Writer(...);

//void OfStreams_trancator() {return;}

//EntIndices(sfaces_sequence, CellNumbs, FES, Face_Entropy_Median, Face_Entropy_Skrew)
/*vector<int> EntIndices( std::vector<unsigned int> &s_faces_sequence, std::vector<unsigned int> const& CellNumbs, Eigen::SparseMatrix<double> const& FES, double &Face_Entropy_Median, double &Face_Entropy_Skrew)
{
    vector<int> TJsTypes(CellNumbs.at(1),0);
    map<unsigned int, unsigned int> res; // Here 100 is an arbitrary number of Edge types
    map<unsigned int, unsigned int>::iterator sit; // Special iterator for this map
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configurational_Face_Entropy = 0;
    Face_Entropy_Median = 0; Face_Entropy_Skrew = 0;

    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
//    for (auto sfe : TJsTypes) cout << sfe << "\t" ; cout << endl;

    J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
    J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
    J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
    J0 = CellNumbs.at(1) - J1 - J2 - J3;
    Jall = CellNumbs.at(1);
// Conversion from numbers to fractions
// (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
    j0 = J0/Jall; j1 = J1/Jall; j2 = J2/Jall; j3 = J3/Jall;
    double j0s = j0, j1s = j1, j2s = j2, j3s = j3;
    /// using values with pow(10,-10) instead of 0s!
//    if (j0s == 0) j0s = pow(10,-30); if (j1s == 0) j1s = pow(10,-30); if (j2s == 0) j2s = pow(10,-30); if (j3s == 0) j3s = pow(10,-30); //Gives 0 in entropy!
    if (j0s != 0) j0s = j0* log2(j0); if (j1s != 0) j1s = j1* log2(j1); if (j2s != 0) j2s = j2* log2(j2); if (j3s != 0) j3s = j3* log2(j3); //Gives 0 in entropy!
    /// Configuration Entropy related with Faces
//    Configurational_Face_Entropy = - (j0s* log2(j0s) + j1s* log2(j1s) + j2s* log2(j2s) + j3s* log2(j3s));
    Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

    /// Median part in the entropy decomposition
    j_types_fractions = {j0, j1, j2, j3}; /// using values with pow(10,-10) instead of 0s!
    if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
        Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(j0 * j1 * j2 * j3);
    } else Face_Entropy_Median = 0.0;

    /// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
    for (int j = 0; j < j_types_fractions.size(); j++)
        for (int i = 0; i < j; i++)
            if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
                Face_Entropy_Skrew +=
                        -(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
                        log2(j_types_fractions[i] / j_types_fractions[j]);
            } else Face_Entropy_Skrew += 0.0;

    /// Complexity/ Informativeness
    // Srand and Smax reading from files
    string input_filename_SstrRand = "Random_Entropy_100.txt"s, input_filename_SstrMAX = "Maximum_Entropy_100.txt"s;
    string input_RandomEntropy_dir = output_dir + input_filename_SstrRand, input_MAXEntropy_dir = output_dir + input_filename_SstrMAX;
    char* RandomEntropy_dir = const_cast<char*>(input_RandomEntropy_dir.c_str()); char* MAXEntropy_dir = const_cast<char*>(input_MAXEntropy_dir.c_str()); // From string to char for the passing folder path to a function
    // Format of the elements:: special Face fraction _ Conf entropy value _ Mean entropy value (log(p1*p2*..*pn))
    vector<vector<double>>  RandomEntropy = VectorVectors4Reader(RandomEntropy_dir);
    vector<vector<double>>  MAXEntropy = VectorVectors4Reader(MAXEntropy_dir);
    //for ( auto tpl : MAXEntropy) cout << tpl[0] << " " << tpl[1] << " " << tpl[2] << endl;

    /// Index of complexity:
    double RandomEntropy_p = 0, MAXEntropy_p = 0;
    vector<double> Delta_MAX;
    double sff =(double) s_faces_sequence.size()/ CellNumbs.at(2); // special face fraction
    for ( auto tpl : MAXEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Smax = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    for ( auto tpl : RandomEntropy)  Delta_MAX.push_back(abs(sff - tpl.at(0)));
    auto numb_Srand = std::min_element(std::begin(Delta_MAX), std::end(Delta_MAX)) - std::begin(Delta_MAX); // gives index of the max element
    Delta_MAX.clear();
    /// Informativeness parameter
//REPAIR    cout << numb_Smax << " " << numb_Srand << endl;
//REPAIR    cout << Configurational_Face_Entropy << " " << RandomEntropy[numb_Srand][1] << " " << MAXEntropy[numb_Smax][1] << endl;
    informativeness =  ( Configurational_Face_Entropy - RandomEntropy[numb_Srand][1])/ ( MAXEntropy[numb_Smax][1] - RandomEntropy[numb_Srand][1]);
    if (informativeness > 1) informativeness = 1;
    /// ====== Data output =====================>
    /// Opening of the output streams
    string TJs_output_filename = "TJsLab_TJsTypes.txt"s, Entropy_output_filename = "TJsLab_ConTJsEntropy.txt"s,
            output_TJs_dir = output_dir + TJs_output_filename, output_Entropy_dir = output_dir + Entropy_output_filename;
    char* cTJs_dir = const_cast<char*>(output_TJs_dir.c_str()); char* cEntropy_dir = const_cast<char*>(output_Entropy_dir.c_str()); // From string to char for the passing folder path to a function

    ofstream OutTJsFile; OutTJsFile.open(cTJs_dir, ios::app);
    double special_faces_fraction =(double) s_faces_sequence.size()/ CellNumbs.at(2);
//    cout << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t"<< Face_Entropy_Median << endl;

    OutTJsFile << special_faces_fraction << "\t" << j0 << "\t" << j1 << "\t" << j2 << "\t" << j3 << Configurational_Face_Entropy << "\t" << Face_Entropy_Median << "\t\t" << Face_Entropy_Skrew << endl;
    OutTJsFile.close();

//    int b = 0;
//    for (unsigned int it : TJsTypes) res[b++] = it;

    return TJsTypes;
}

 */