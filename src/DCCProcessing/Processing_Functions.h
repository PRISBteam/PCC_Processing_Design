///================================ A part of the DCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifications for k-Cells               *
 *  in the DCC Processing module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                            **/
///================================================================================================================================///

/// (1) Totally random choice of element with # New2CellNumb from the list of numbers {0 to SCellsNumb}
unsigned int NewCellNumb_R(unsigned int OCellsNumb){ // Random generation of a 2-Cell number
    unsigned int New2CellNumb;
    New2CellNumb = rand() % (OCellsNumb-1); // Random generation of the boundary number in the range from 0 to OrdinaryCellNumbs.size()-1
    return New2CellNumb;
}

/// (2) The Random generation process function
int Processing_Random( std::vector<unsigned int> &S_Vector,  std::vector<unsigned int> &s_faces_sequence, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs) {
///================================================================= 'R' =======================================================================////
/// ====================================================>  Random generation process  <========================================================////
///===========================================================================================================================================////

/// Numerates newly created Faces during the random generation process
    long unsigned int numerator = 0;
///=============================================================================================================================================////
/// =====> Initial initialisation with the previous calculation step (if any) based on the _special_faces_sequence_
    // OrdinaryCellNumbs is just a tricky way to acceleration the random process
    std::vector<unsigned int> OrdinaryCellNumbs(CellNumbs.at(2), 1); // Vector of the size equal to the total number of faces in DCC initialised with '0's
    for( unsigned int lit = 0; lit < OrdinaryCellNumbs.size(); lit++) OrdinaryCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces
    for( unsigned int itr : S_Vector) if(itr != 0) OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + itr); // !!! Delete its element from the vector decreasing its size BUT

    // unsigned int OCellAmount = std::count(SpecialCellNumbs->begin(), SpecialCellNumbs->end(), 0);
    double ordinary_faces_fraction = OrdinaryCellNumbs.size()/ (double) CellNumbs.at(2);
    double special_faces_fraction = 1.0 - ordinary_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) return 0;

/// ================= Loop before max_sFaces_fraction fractions of special cells  =======================>
    /// The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
    srand((unsigned) time(NULL));

    do { // do{ ... }while(output_step) loop starting point
        int New2CellNumb = 0;
        New2CellNumb = NewCellNumb_R(OrdinaryCellNumbs.size());
        /// Random generation of types !!! with IDs < number_of_types
        int NewFaceType = 1;
        if(number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types

        /// Changes in vectors from Main function
        S_Vector.at(OrdinaryCellNumbs.at(New2CellNumb)) = NewFaceType;
        s_faces_sequence.push_back(OrdinaryCellNumbs.at(New2CellNumb));

        // It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
        OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size BUT

        // Special and Ordinary Faces fraction calculation
        ordinary_faces_fraction = OrdinaryCellNumbs.size() / (double) CellNumbs.at(2);
        special_faces_fraction = 1.0 - ordinary_faces_fraction;
// REPAIR        for (auto sfe : s_faces_sequence) cout << sfe << "\t"; cout << endl;
//        cout << special_faces_fraction << "\t" << max_sFaces_fraction << endl;
    }while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process

/// Closing and deleting
    //Remove all elements anf free the memory from the probe OrdinaryCellNumbs vector
    OrdinaryCellNumbs.clear();
    OrdinaryCellNumbs.shrink_to_fit();

    return 0;
} // end  of Random

/// (3) Maximum entropy generation process
int Processing_maxEntropy(std::vector<unsigned int>  &S_Vector, std::vector<unsigned int>  &s_faces_sequence, double max_sFaces_fraction, int number_of_types, std::vector<unsigned int> &CellNumbs, std::vector<char*> const paths) {
///=============================================================================================================================================////
///========================================================================= 'S' =================================================================////
/// ====================================================>  Maximum entropy production process   <========================================================////
///=============================================================================================================================================////
    typedef Eigen::SparseMatrix<double> SpMat;

    /// Local functions definition
    //vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, map<unsigned int, unsigned int> Edges_TypesMap);
    vector<double> GBIndex(unsigned int face_number, Eigen::SparseMatrix<double> const& FES, std::vector<int> Edges_TypesMap);
    vector<int> EdgesTypesCalc(std::vector<unsigned int> const& CellNumbs, vector<unsigned int> &SpecialCellMap, Eigen::SparseMatrix<double> const& FES);

    // Map and its iterator for Faces
//    map<unsigned int,unsigned int> SpecialCellMap; // Mapping [k]->[l] from the set of all 2-Cells in the initial DCC to the set of newly generated special cells
//    map<unsigned int, unsigned int>::iterator sit; // Special iterator for such map
    // Vectors for Edges and Edges-related configurational entropy
    vector<int> TJsTypes(CellNumbs.at(1) + 1,0);
    vector<double> EntropyIncreaseList;
    double J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
    double Configuration_Face_Entropy = 0, Face_Entropy_Median = 0, Face_Entropy_Skrew = 0;

    // Obtaining Faces (coloumns) - Edges (rows) incidence matrix from file (paths.at(5))
    SpMat FES(CellNumbs.at(1),CellNumbs.at(2));
    FES = SMatrixReader(paths.at(5), CellNumbs.at(1), CellNumbs.at(2)); // Edges-Faces

    /// Numerates newly created Faces during the generation process
    long unsigned int numerator = 0;

    double special_faces_fraction = s_faces_sequence.size()/ (double) CellNumbs.at(2);
    double ordinary_faces_fraction = 1.0 - special_faces_fraction;
    if (special_faces_fraction >= max_sFaces_fraction) { cout << "(!) Early exit: Fraction of special faces is already >= MAX values from config.txt" << endl; return 0; }
    else if (special_faces_fraction < 0.05) {
        /// ================> Initial Face seeds (initial state for the MAX entropy production algorithm)
        Processing_Random(S_Vector, s_faces_sequence, 0.05, number_of_types, CellNumbs);
    }

///=============================================================================================================================================////
/// ================= S_MAX loop over the fractions of special cells [0.05,1] =======================>
    do { // do{ ... }while(output_step) loop starting point
//            unsigned int OCellAmount = ordinary_faces_fraction*(double) CellNumbs.at(2);
//            special_faces_fraction = 1.0 - ordinary_faces_fraction;
        /// ========== Calculation configuration entropy at each calculation step ===============>
        // Zeroing coefficients
        std::fill(TJsTypes.begin(), TJsTypes.end(), 0);
        /// TJs type calculations
        TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
//REPAIR               for(auto kl : TJsTypes) cout << " " <<kl ; cout << endl; exit(10);
        J0 = 0, J1 = 0, J2 = 0, J3 = 0, Jall = 0, j0 = 0, j1 = 0, j2 = 0, j3 = 0;
        Configuration_Face_Entropy = 0;

        J1 = std::count(TJsTypes.begin(), TJsTypes.end(), 1);
        J2 = std::count(TJsTypes.begin(), TJsTypes.end(), 2);
        J3 = std::count(TJsTypes.begin(), TJsTypes.end(), 3);
        J0 = CellNumbs.at(1) - J1 - J2 - J3;
        Jall = CellNumbs.at(1);
// Conversion from numbers to fractions | (!) log2 means binary (or base-2) logarithm and we use "-" for fractions to make the value positive
        j0 = J0 / Jall;
        j1 = J1 / Jall;
        j2 = J2 / Jall;
        j3 = J3 / Jall;
        /// Configuration Entropy related with Faces
        Configuration_Face_Entropy = -(j0 * log2(j0) + j1 * log2(j1) + j2 * log2(j2) + j3 * log2(j3));

        /// Calculation of all Edges types :: Edges_TypesMap: [Edge number] --> [Edge type]
//        std::map<unsigned int, unsigned int> Edges_TypesMap; // Map: [Edge number] --> [Edge type]
//        Edges_TypesMap = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);

        /// ======== Loop over all Faces ===========>
        map<double,unsigned int> EList_to_FaceIDs; // Map: [Entropy increase] -> [Face number]
        EList_to_FaceIDs.clear(); EntropyIncreaseList.clear(); // Clearence of the Map and Vector of all configuration entropy increases
        for (unsigned int k = 0; k < CellNumbs.at(2); k++) { //loop over all Faces in DCC
            double J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;

            if (S_Vector.at(k) == 0) { // Loop over each still ORDINARY element neighbours
                J00 = 0, J0N = 0, J10 = 0, J1N = 0, J20 = 0, J2N = 0, J30 = 0, J3N = 0, CFace_EntropyIncrease = 0;
                /// GBIndex calculation :: gives a vector with the number of Edges of each type
                vector<double> j_types_neigh_fractions = GBIndex(k, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
                //REPAIR    for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;

                // Values before conversion
                J00 = j_types_neigh_fractions.at(0);
                J10 = j_types_neigh_fractions.at(1);
                J20 = j_types_neigh_fractions.at(2); // cout << " J00= " << J00<< " J10= " << J10 << " J20= " << J20 << endl;

                // Values after conversion
//Oldie//                J0N = J10;
                J1N = J00;
                J2N = J10;
                J3N = J20;
                // The entropy increase calculation for a given Face
//                CFace_EntropyIncrease =
                        abs(J0 - J00) * log2(abs(J0 - J00)+pow(10,-10)) + abs(J1 - J10 + J1N) * log2(abs(J1 - J10 + J1N)+pow(10,-10)) +
                        abs(J2 - J20 + J2N) * log2(abs(J2 - J20 + J2N)+pow(10,-10)) + abs(J3 + J3N) * log2(abs(J3 + J3N)+pow(10,-10)); //  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;
                CFace_EntropyIncrease =
                        (J0 - J00) * log2((J0 - J00)+pow(10,-10)) + (J1 - J10 + J1N) * log2((J1 - J10 + J1N)+pow(10,-10)) +
                        (J2 - J20 + J2N) * log2((J2 - J20 + J2N)+pow(10,-10)) + (J3 + J3N) * log2((J3 + J3N)+pow(10,-10)); //  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;

//Oldie//                CFace_EntropyIncrease = (J0 - J00 + J0N) * log2(J0 - J00 + J0N) + (J1 - J10 + J1N) * log2(J1 - J10 + J1N) + (J2 - J20 + J2N) * log2(J2 - J20 + J2N); //  cout  << "\t\t" <<  CFace_EntropyIncrease << "\t\t" << endl;

                // The result of the one iteration
                EntropyIncreaseList.push_back(CFace_EntropyIncrease); //        for(auto kl :EntropyIncreaseList) cout << kl << endl;
                EList_to_FaceIDs[CFace_EntropyIncrease] = k;
            } // if OrdinaryCells (S_Vector.at(Face) == 0)
        } // for (..k < CellNumbs.at(2)..)

        double New2CellNumb = 0, NewFaceType = 1; // Only one possible Face type (binary model)
        /// Number of element giving the maximum increase in configuration entropy at its conversion
        unsigned int dist = std::max_element(std::begin(EntropyIncreaseList), std::end(EntropyIncreaseList)) - std::begin(EntropyIncreaseList); // gives index of the max element
        New2CellNumb = EList_to_FaceIDs[EntropyIncreaseList.at(dist)];
//REPAIR        cout << s_faces_sequence.size() << "   " << New2CellNumb << endl;

        // Then all the corresponding maps chain
        S_Vector.at((unsigned int) New2CellNumb) = 1; // Replace the chosen element with 0 instead of 1 in the State Faces vector
//////////////////        OrdinaryCells.at(New2CellNumb) = 0; // Replace the chosen element with 0 instead of 1 in the Original Faces vector
        //SpecialCellMap[numerator++] = New2CellNumb ; // Assign new elements to the map and increase numerator
//       SpecialCellMap.insert( std::pair<unsigned int, unsigned int>(,New2CellNumb));
        if(find(s_faces_sequence.begin(),s_faces_sequence.end(),New2CellNumb) == s_faces_sequence.end())
            s_faces_sequence.push_back(New2CellNumb);
        //REPAIRS/////////
        //       for(auto kl :s_faces_sequence) cout << " " <<kl ; cout << endl; //exit(10);


        // Special and Ordinary Faces fraction calculation
        special_faces_fraction = s_faces_sequence.size() / (double) CellNumbs.at(2);
        ordinary_faces_fraction = 1.0 - special_faces_fraction;

        if ((int) special_faces_fraction*CellNumbs.at(2) % (int) 0.01*max_sFaces_fraction*CellNumbs.at(2) == 1) cout << special_faces_fraction << endl;
    }while(special_faces_fraction < max_sFaces_fraction); /// End of the Random generation process

//REPAIR    cout << "in_new:" <<endl; for (auto itd : s_faces_sequence) cout << itd << endl;
/// Closing and deleting
    //Remove all elements anf free the memory from the probe SpecialCells and vector
//    SpecialCellMap.clear();
    //   SpecialCellMap.shrink_to_fit();

    return 0;
} // end of S_max
//REPAIRS/////////        for(auto kl :S_Vector) cout << " " <<kl ; cout << endl; //exit(10);

/// (4) DDRX process
int Processing_DDRX(std::vector<unsigned int>  &State_Vector, std::vector<unsigned int>  &special_faces_sequence, double max_sFaces_fraction, std::vector<char*> const paths, int number_of_types, std::vector<unsigned int> &CellNumbs) {
///=============================================================================================================================================////
///========================================================================= 'D' =================================================================////
/// ==================================================================>  DDRX process   <========================================================////
///=============================================================================================================================================////
    /// Function
    tuple<double, double, double> find_aGBseed(unsigned int Facenumb, std::vector<char*> const paths, std::vector<unsigned int> & CellNumbs, vector<tuple<double, double, double>> & AllSeeds_coordinates);
    /// Streams
    ofstream NewSeedsStream;
    NewSeedsStream.open(paths.at(8), ios::trunc);
    NewSeedsStream << "New generated seeds\t" << endl;
    NewSeedsStream.close();

    /// Vector
//    vector<Triplet<double>> AllSeeds_coordinates, NewSeeds_coordinates;
    vector<tuple<double, double, double>> AllSeeds_coordinates;
    tuple<double, double, double> NewSeed_coordinates = make_tuple(0, 0, 0);

    /// Face Seeds reading from file Seeds.txt
    AllSeeds_coordinates = TuplesReader(paths.at(7));

/// ============== Constants and model parameters ===================
    double kb = 1.23* pow(10,-23);
    double GSav = 200.0*pow(10,-9), Temp = 300, Stress = 4.0*pow(10.0,6), Rcr = 0.5*GSav/10.0,
            Q_energy = 213.0*pow(10,1), R_const = 8.31;
/// Conversion probability
    double prob_seed = exp(-Q_energy/(R_const*Temp));
    //    double prob_seed = exp(- Stress*(13.0*pow(Rcr,3)/3.0)/(kb*Temp));
    cout << "DDRX module:: prob_seed =\t" << prob_seed << endl;

    srand((unsigned) time(NULL)); // seed for random
    // Output stream opening
    NewSeedsStream.open(paths.at(8), ios::app);

    for(unsigned int fnumber = 0; fnumber < CellNumbs.at(2); ++fnumber) {
        double rv = (rand() / (RAND_MAX + 1.0));
        //   cout << rv << "\t" << prob_seed << endl;
        if (rv <= prob_seed) {
            NewSeed_coordinates = find_aGBseed(fnumber, paths, CellNumbs, AllSeeds_coordinates);
            NewSeedsStream << fnumber << "\t" << get<0>(NewSeed_coordinates) << "\t" << get<1>(NewSeed_coordinates) << "\t" << get<2>(NewSeed_coordinates) << endl;
        }
    }
    NewSeedsStream.close();

    return 0;
}

/*
*     /// Maximal fraction for simulation loop max_sFaces_fraction = [0,1]
    max_sFaces_fraction = 0.3;
    /// Step for output (structural analysis and data output will be performed after each output_step calculation steps (= number of newly converted elements))
    int output_step = 300;

*/

/*
 * /// At first, random generation of small amount about 0.05 of special seeds
    do { // do{ ... }while(0.05) loop for seed points
        int New2CellNumb = 0;
        New2CellNumb = NewCellNumb_R(OrdinaryCellNumbs.size());
        /// Random generation of types !!! with IDs < number_of_types
        int NewFaceType = 1;
        if(number_of_types > 1) NewFaceType = rand() % number_of_types; // Random chose for the chosen boundary to be assigned over all special types

        /// Changes in vectors from Main function
        State_Vector->at(OrdinaryCellNumbs.at(New2CellNumb)) = NewFaceType;
        special_faces_sequence->push_back(OrdinaryCellNumbs.at(New2CellNumb));
        /// Special Faces map
        SpecialCellMap[OrdinaryCellNumbs.at(New2CellNumb)] = numerator++; // Assign new elements to the map and increase numerator
        long unsigned int sit = SpecialCellMap[OrdinaryCellNumbs.at(New2CellNumb)]; // Just a useful variable for the number of newly converted Face (= numerator--)

        // It is a key tricky point for the fast Random generation process: vector decreases after each turn from an ordinary to special face BUT we chose all the boundary one by one because their initial numeration stored in the SpecialCellNumbs[] vector !
        OrdinaryCellNumbs.erase(OrdinaryCellNumbs.begin() + New2CellNumb); // !!! Delete its element from the vector decreasing its size BUT

        // Special and Ordinary Faces fraction calculation
        ordinary_faces_fraction = OrdinaryCellNumbs.size() / (double) CellNumbs->at(2);
        special_faces_fraction = 1.0 - ordinary_faces_fraction;

    }while(ordinary_faces_fraction > 0.95); /// End of the Random generation process

 */