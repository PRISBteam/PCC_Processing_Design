///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new "secondary" identifications of k-Cells   **/
/**  in the PCC Processing module. It makes them "fractured" and can go alongside their special or ordinary primal types.         **/
///==============================================================================================================================///


/// #1# Kinetic function for multiple cracking
std::vector <unsigned int> PCC_Kinetic_cracking(std::vector <double> &face_elastic_energies, std::vector<unsigned int> &frac_sfaces_sequence, macrocrack &large_crack, Eigen::SparseMatrix<double> const& AFS, Eigen::SparseMatrix<double> const& FES) {

    std::vector <unsigned int> crack_faces_sequence, S_crackVector(CellNumbs.at(2), 0); // sequence of the cracked Faces and State vector for cracks
    vector<double> TJsTypes(CellNumbs.at(1), 0), TJsCrackTypes(CellNumbs.at(1), 0); // calculations of TJs related to inclusions and cracks
    vector<double> newCrack_neigh_TJs, newCrack_neigh_Faces; // only for neighbouring faces
    std::vector <double> Face_weight(CellNumbs.at(2),0), Face_rGO_index(CellNumbs.at(2),0), Face_crack_index(CellNumbs.at(2),0); // face weight (number of neighbours) and indexes for all GBs
    std::vector <double> lsc_rGO_energy(CellNumbs.at(2),0), lsc_crack_energy(CellNumbs.at(2),0); // energies related with the local stress concentrators Bl and Cl
    std::vector <double> Face_energy(CellNumbs.at(2),0), Face_current_energy(CellNumbs.at(2),0); // energies for all GBs
    double crack_fraction = 0.0;
    double total_crack_energies = 0.0;

    double sface_energy_matrix = 2.0, sface_energy_rGO = 1.0, sface_energy_aggl = 0.4;
    /// still arbitrarily chosen!
    double  lsc_rGO = 4.0*sface_energy_rGO, lsc_crack = 2.0*sface_energy_rGO; // energy values related with the stress concentrators //    double kB = 1.3807*pow(10,-23); // Boltzmann constant

    /// Initial energies of special faces
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) = sface_energy_matrix * face_areas_vector.at(i) * sample_dimensions.at(0) * sample_dimensions.at(1); // assignments of the matrix surface energies
    for (auto ks : frac_sfaces_sequence) Face_energy.at(ks) = sface_energy_rGO * face_areas_vector.at(ks) * sample_dimensions.at(0) * sample_dimensions.at(1); // assignments of the rGO defect surface energies

    ///Calculation of the TJs types
    TJsTypes = EdgesTypesCalc(CellNumbs, frac_sfaces_sequence, FES);

    /// GB_indices calculation
    for (unsigned int fn = 0; fn < CellNumbs.at(2); ++fn) {
        vector<double> j_types_neigh_fractions = GBIndex(fn, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
//REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;

        // rGO (Bl) index calculation
        Face_rGO_index.at(fn) = (j_types_neigh_fractions.at(0) + 2.0 * j_types_neigh_fractions.at(1) + 3.0 * j_types_neigh_fractions.at(2)) / 3.0;

    } // end of for ( fn < CellNumbs.at(2)) GB_indices calculation

    /// Local stress concentrators and related energies
    // Face weights calculation
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i)
        for (int l = 0; l < AFS.rows(); l++) // Loop over all Faces
            if (AFS.coeff( l, i) == 1) Face_weight.at(i)++;

    // Energies
    for (unsigned int j = 0; j < CellNumbs.at(2); ++j) lsc_rGO_energy.at(j) = lsc_rGO * face_areas_vector.at(j) * sample_dimensions.at(0) * sample_dimensions.at(1) * (Face_rGO_index.at(j) / Face_weight.at(j));
    for (unsigned int k = 0; k < CellNumbs.at(2); ++k) lsc_crack_energy.at(k) = lsc_crack * face_areas_vector.at(k) * sample_dimensions.at(0) * sample_dimensions.at(1) * (Face_crack_index.at(k) / Face_weight.at(k));

    /// Final surface energy - a single calculation of Bl concentrators effect on face energies
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) -= lsc_rGO_energy.at(i);

    /// Crack face sequence creation
    // (1) FAST case : cracks do NOT affect the process of the further fracture`
    //     for (unsigned int i = 0; i < CellNumbs.at(2); ++i){ ;}

    // (2) SLOW case : cracks affect the process of the further fracture`
    /// "Infinite" loop before there are new fractured faces
    bool is_fractured = 1; // id signifying that some GBs was fractured - the infinite loop stops when no one no crack is appeared during an iteration
    do {
        is_fractured = 0; // initial assumption

        double NewCrackNumb = 0;
        /// Loop over all faces in PCC
        for (unsigned int i = 0; i < CellNumbs.at(2); ++i) { // crack loop over all the GBs
            // Update for local stress concentrators (Cl) and related energies
            for (unsigned int k = 0; k < CellNumbs.at(2); ++k)
                lsc_crack_energy.at(k) = lsc_crack * face_areas_vector.at(k) * sample_dimensions.at(0) * sample_dimensions.at(1) * (Face_crack_index.at(k) / Face_weight.at(k));

            /// Current GB energies
            /// New surface energy - calculation of Cl concentrators effect on face energies
            for (unsigned int l = 0; l < CellNumbs.at(2); ++l)
                Face_current_energy.at(l) = Face_energy.at(l) - lsc_crack_energy.at(l);

            /// New surface energy - the effect of external elastic stresses (!)
            for (unsigned int m = 0; m < CellNumbs.at(2); ++m)
                Face_current_energy.at(m) -= face_elastic_energies.at(m);

            /// New fractured face of elements
            NewCrackNumb = i;
            /// if the face is not already fractured and its full energy less than 0 (including the "-" full external elastic energy)
            if (S_crackVector.at(NewCrackNumb) != 1 && Face_current_energy.at(NewCrackNumb) < 0) {
                crack_faces_sequence.push_back(NewCrackNumb); // if the element is not in the sequence - add
                is_fractured = 1; // there is at least one fractured GB
                S_crackVector.at(NewCrackNumb) = 1;
                if (find(frac_sfaces_sequence.begin(), frac_sfaces_sequence.end(), NewCrackNumb) == frac_sfaces_sequence.end()) // is matrix boundary without inclusion
                    Face_energy.at(NewCrackNumb) = 2.0 * sface_energy_matrix * face_areas_vector.at(NewCrackNumb) * sample_dimensions.at(0) * sample_dimensions.at(1);
                else Face_energy.at(NewCrackNumb) = 2.0 * sface_energy_rGO * face_areas_vector.at(NewCrackNumb) * sample_dimensions.at(0) * sample_dimensions.at(1);

                /// Update for total surface energy of all the fractured GBs
                total_crack_energies += Face_energy.at(NewCrackNumb);

                /// New cracked neighbours of the converted NewCrackNumb Face
                newCrack_neigh_TJs.clear();
                newCrack_neigh_Faces.clear();
                for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the Edges
                    if (FES.coeff(k, NewCrackNumb) == 1) newCrack_neigh_TJs.push_back(k);
                // update for cracked TJs types
                if (newCrack_neigh_TJs.size() > 0) for (auto itr: newCrack_neigh_TJs) TJsCrackTypes.at(itr)++;

                for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
                    // update for crack neighbour faces
                    if (AFS.coeff(NewCrackNumb, m) == 1 && S_crackVector.at(m) == 0) newCrack_neigh_Faces.push_back(m);

                // new GB index calculation
                if (newCrack_neigh_Faces.size() > 0) {
                    for (auto fn: newCrack_neigh_Faces) {
                        vector<double> j_types_crack_n_fractions = GBIndex(fn, FES,TJsCrackTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
                        /// New crack index
                        Face_crack_index.at(fn) = (j_types_crack_n_fractions.at(0) + 2.0 * j_types_crack_n_fractions.at(1) + 3.0 * j_types_crack_n_fractions.at(2)) / 3.0;
                    } // end for (auto fn: new_crack_neigh_Faces)
                } // end if (new_crack_neigh_Faces.size() > 0)

            } /// end if (S_crackVector.at(NewCrackNumb) != 1 && Face_current_energy.at(NewCrackNumb) < 0)

            crack_fraction = crack_faces_sequence.size() / (double) CellNumbs.at(2);
//REPAIR    cout << "crack_fraction = " << crack_fraction <<"  " << "Total crack surfaces energy = " << total_crack_energies << endl;
        } // end for ( i < CellNumbs.at(2) )

    } while (is_fractured);

    large_crack.Set_multiple_cracking_energy(total_crack_energies); //Set_multiple_cracking_energy(double total_energy)

    return crack_faces_sequence;
} /// end of Kinetic_cracking
/// -----------------------------------------------------------------------------------------------------------------///

std::vector <unsigned int> DCC_Kinematic_cracking(std::vector<unsigned int> &s_faces_sequence, Eigen::SparseMatrix<double> const& AFS, Eigen::SparseMatrix<double> const& FES) {

    std::vector <unsigned int> crack_faces_sequence, S_crackVector(CellNumbs.at(2), 0); // sequence of the cracked Faces
    vector<double> TJsTypes(CellNumbs.at(1), 0), TJsCrackTypes(CellNumbs.at(1), 0);
    vector<double> newCrack_neigh_TJs, newCrack_neigh_Faces; // only for neighbouring faces
    std::vector <double> Face_weight(CellNumbs.at(2),0), Face_rGO_index(CellNumbs.at(2),0), Face_crack_index(CellNumbs.at(2),0); // indexes for all GBs
    std::vector <double> lsc_rGO_energy(CellNumbs.at(2),0), lsc_crack_energy(CellNumbs.at(2),0); // energies related with the local stress concentrators
    std::vector <double> Face_energy(CellNumbs.at(2),0), Face_current_energy(CellNumbs.at(2),0); // energies for all GBs
    double crack_fraction = 0.0;
//old    double sface_energy_rGO = 0.5*pow(10,-3), sface_energy_matrix = 1.0*pow(10,-3);
    double sface_energy_matrix = 2.0, sface_energy_rGO = 1.0, sface_energy_aggl = 0.4;

    double  lsc_rGO = 4.0*sface_energy_rGO, lsc_crack = 2.0*sface_energy_rGO; // energy values related with the stress concentrators
    double kB = 1.3807*pow(10,-23); // Boltzmann constant

    /// Initial energies
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) = sface_energy_matrix; // assignments of the matrix surface energies
    for (auto ks : s_faces_sequence) Face_energy.at(ks) = sface_energy_rGO; // assignments of the rGO defect surface energies

    //Calculation of the TJs types
    TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
    /// GB_indices calculation
    for (unsigned int fn = 0; fn < CellNumbs.at(2); ++fn) {
        vector<double> j_types_neigh_fractions = GBIndex(fn, FES, TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
        //REPAIR   for(auto kl : j_types_neigh_fractions) cout << " " <<kl ; cout << endl;

        // rGO index calculation
        Face_rGO_index.at(fn) = (j_types_neigh_fractions.at(0) + 2.0 * j_types_neigh_fractions.at(1) + 3.0 * j_types_neigh_fractions.at(2)) / 3.0;

    } // end of for ( fn < CellNumbs.at(2))

    /// Local stress concentrators and related energies
    // Face weights calculation
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i)
        for (int l = 0; l < AFS.rows(); l++) // Loop over all Faces
            if (AFS.coeff( l, i) == 1) Face_weight.at(i)++;

    // Energies
    for (unsigned int j = 0; j < CellNumbs.at(2); ++j) lsc_rGO_energy.at(j) = lsc_rGO*(Face_rGO_index.at(j)/Face_weight.at(j));
    for (unsigned int k = 0; k < CellNumbs.at(2); ++k) lsc_crack_energy.at(k) = lsc_crack*(Face_crack_index.at(k)/Face_weight.at(k));

    /// Final surface energy
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) Face_energy.at(i) = Face_energy.at(i) -= lsc_rGO_energy.at(i);

    /// Crack face sequence creation
    /// (1) FAST case : cracks do NOT affect the process of the further fracture`
    //     for (unsigned int i = 0; i < CellNumbs.at(2); ++i){ ;}

    /// (2) SLOW case : cracks affect the process of the further fracture`
    double NewCrackNumb = 0;
    for (unsigned int i = 0; i < CellNumbs.at(2); ++i) { // crack loop over all the GBs
        /// New surface energy
        // Local stress concentrators and related energies
        for (unsigned int k = 0; k < CellNumbs.at(2); ++k) lsc_crack_energy.at(k) = lsc_crack*(Face_crack_index.at(k)/Face_weight.at(k));
        // Current surface energy
        for (unsigned int l = 0; l < CellNumbs.at(2); ++l) Face_current_energy.at(l) = Face_energy.at(l) - lsc_crack_energy.at(l);

        /// Number of element
        NewCrackNumb = std::min_element(std::begin(Face_energy), std::end(Face_energy)) - std::begin(Face_energy); // gives index of the max element
        if (find(crack_faces_sequence.begin(), crack_faces_sequence.end(), NewCrackNumb) == crack_faces_sequence.end()) {
            crack_faces_sequence.push_back(NewCrackNumb); // if the element is not in the sequence - add
            S_crackVector.at(NewCrackNumb) = 1;
            Face_energy.at(NewCrackNumb) = 100000.0; /// arbitrarily large value !!!
        } // end if

        /// Neighbours of the converted Face
        newCrack_neigh_TJs.clear();
        newCrack_neigh_Faces.clear();
        for (int k = 0; k < CellNumbs.at(1); ++k) // Loop over all the Edges
            if (FES.coeff(k, NewCrackNumb) == 1) newCrack_neigh_TJs.push_back(k);
        if (newCrack_neigh_TJs.size() > 0) for (auto itr: newCrack_neigh_TJs) TJsCrackTypes.at(itr)++;

        for (int m = 0; m < CellNumbs.at(2); ++m) // Loop over all the Faces
            if (AFS.coeff(NewCrackNumb, m) == 1 && S_crackVector.at(m) == 0) newCrack_neigh_Faces.push_back(m);

        if (newCrack_neigh_Faces.size() > 0) {
            for (auto fn: newCrack_neigh_Faces) {
                vector<double> j_types_crack_n_fractions = GBIndex(fn, FES,TJsCrackTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                /// New crack index
                Face_crack_index.at(fn) = (j_types_crack_n_fractions.at(0) + 2.0 * j_types_crack_n_fractions.at(1) +
                                           3.0 * j_types_crack_n_fractions.at(2)) / 3.0;
            } // end for (auto fn: new_crack_neigh_Faces)
        } // end if (new_crack_neigh_Faces.size() > 0)

        crack_fraction = crack_faces_sequence.size() / (double) CellNumbs.at(2);
//REPAIR cout << crack_fraction << endl;
    } // end for ( i < CellNumbs.at(2) )

    return crack_faces_sequence;
} /// end of Kinetic_cracking
/// -----------------------------------------------------------------------------------------------------------------///
