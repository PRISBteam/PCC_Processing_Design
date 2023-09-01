///================================ A part of the PCC Processing module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new "secondary" identifications of k-Cells   **/
/**  in the PCC Processing module. It makes them "fractured" and can go alongside their special or ordinary primal types.         **/
///==============================================================================================================================///
// #include "Processing_Induced_functions.h"

std::vector <unsigned int> PCC_Kinematic_cracking(int cell_type, std::vector<unsigned int> &s_faces_sequence, std::vector<vector<int>> &Configuration_cState, std::vector<vector<double>> const &max_cfractions_vectors) {
    std::vector<unsigned int> crack_faces_sequence;// output of the function: sequence of the cracked Faces
    std::vector<int> S_cVector(CellNumbs.at(cell_type), 0); // State Vector for fractured cells

    vector<double> TJsTypes(CellNumbs.at(cell_type - 1), 0), TJsCrackTypes(CellNumbs.at(cell_type - 1), 0);
    vector<double> newCrack_neigh_TJs, newCrack_neigh_Faces; // only for neighbouring faces
    std::vector<double> Face_weight(CellNumbs.at(cell_type), 0), Face_inclusion_index(CellNumbs.at(cell_type),0), Face_crack_index(CellNumbs.at(cell_type), 0); // indexes for all GBs
    std::vector<double> lsc_inclusion_energy(CellNumbs.at(cell_type), 0), lsc_crack_energy(CellNumbs.at(cell_type),
                                                                                           0); // energies related with the local stress concentrators (lsc)
    std::vector<double> Cell_energy(CellNumbs.at(cell_type), 0), Cell_current_energy(CellNumbs.at(cell_type),
                                                                                     0); // energies for all GBs
    double agglomeration_fraction; // a-fraction in the PCC
    std::vector<int> agglomeration_SVector(CellNumbs.at(cell_type),
                                           0); // numbers in this State vector are Agglomeration Powers
    std::vector<unsigned int> agglomerations_set; // (pre-initially empty) set of agglomeration numbers
    std::vector<unsigned int> o_faces_set(CellNumbs.at(cell_type), 0), NotFracturedCellNumbs(CellNumbs.at(cell_type),
                                                                                             0); // set of ordinary and not cells faces in a PCC [technical parameter]

    for (unsigned int lit = 0; lit < CellNumbs.at(cell_type); ++lit)
        o_faces_set[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    // deleting of all faces contain inclusions //OLD and WRONG    for (auto istr = S_cVector.begin(); istr != S_cVector.end(); ++istr) if (*istr != 0) //        o_faces_set.erase(o_faces_set.begin() + ks);
    for (auto ks: s_faces_sequence) /// !!! Delete its element from the vector decreasing its size
        o_faces_set.erase(std::find(o_faces_set.begin(), o_faces_set.end(), ks));

    double crack_fraction = 0.0;
/// PARTICULAR VALUES FOR FACE ENERGIES (SHOULD GO AFTER TO THE "kinematic_fracture_data" FILE)
// all surface energies in [J/m^2] (!)
    double sface_energy_matrix = 2.0, sface_energy_inclusion = 1.0, sface_energy_aggl = 0.4;
    double nrm = 1.0 * sface_energy_matrix, nrr = 0, ncm = 2.0 * sface_energy_matrix, ncr = 2.0 * sface_energy_inclusion; // nrm and nrr - coefficients for energies caused by INCLUSIONS in a matrix (m) or inclusion (r) boundary; ncm and ncr - coefficients for energies caused by CRACKS in a matrix (m) or inclusion (r) boundary
    // double  lsc_inclusion = 4.0 * sface_energy_inclusion, lsc_crack = 2.0 * sface_energy_inclusion; // energy values related with the stress concentrators
    double kB = 1.3807 * pow(10, -23); // Boltzmann constant

    /// Initial energies
    for (unsigned int i = 0; i < CellNumbs.at(cell_type); ++i)
        Cell_energy.at(i) = sface_energy_matrix; // assignments of the matrix surface energies

    for (auto ks: s_faces_sequence)
        Cell_energy.at(ks) = sface_energy_inclusion; // assignments of the rGO defect surface energies

    /// Sparse Face-Edge Incidence matrix - reading from the file of the considered PCC
    SpMat FES = SMatrixReader(paths.at(5 + (dim - 3)), CellNumbs.at(cell_type - 1),
                              CellNumbs.at(cell_type)); // Edges-Faces sparse incidence matrix
    /// Sparse Face Adjacency matrix - reading from the file of the considered PCC
    SpMat AFS = SMatrixReader(paths.at(2 + (dim - 3)), (CellNumbs.at(cell_type)),
                              (CellNumbs.at(cell_type))); //all Faces
    AFS = 0.5 * (AFS + SparseMatrix<double>(AFS.transpose()));     //  Full symmetric AFS matrix instead of triagonal

    /// Calculation of the TJs types
    if (dim == 3) TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
    else if (dim == 2) TJsTypes = NodesTypesCalc(CellNumbs, s_faces_sequence, FES); // kind of a State Vector

    /// GB_indices calculation
    for (unsigned int fn = 0; fn < CellNumbs.at(cell_type); ++fn) {

        vector<double> j_types_neigh_fractions = GBIndex(fn, FES,
                                                         TJsTypes); // Types (up to 100 kinds) of the edges incident to the considered Face
        /// inclusion index calculation
        Face_inclusion_index.at(fn) = (j_types_neigh_fractions.at(1) + 2.0 * j_types_neigh_fractions.at(2) +
                                       3.0 * j_types_neigh_fractions.at(3));
    } // end of for ( fn < CellNumbs.at(2))

    /// Local stress concentrators and related energies ///
    /// Face weights calculation
    for (unsigned int i = 0; i < CellNumbs.at(cell_type); ++i)
        for (int l = 0; l < AFS.rows(); l++) // Loop over all Faces
            if (AFS.coeff(l, i) == 1)
                Face_weight.at(i)++;

/// Inclusions' local elastic stress concentrations Energies  (constant during the fracture process)
    for (auto ks: s_faces_sequence) { // with inclusion (special faces)
        lsc_inclusion_energy.at(ks) = nrr * (Face_inclusion_index.at(ks) / (0.75 * Face_weight.at(ks))); // = 0 (!)
//        lsc_crack_energy.at(ks) = ncr * (Face_crack_index.at(ks)/ (0.75*Face_weight.at(ks)));
    }
    for (auto os: o_faces_set) { // without inclusions (ordinary faces)
        lsc_inclusion_energy.at(os) = nrm * (Face_inclusion_index.at(os) / Face_weight.at(os));
//        lsc_crack_energy.at(os) = ncm * (Face_crack_index.at(os)/ Face_weight.at(os));
    }

    /// ++ AGGLOMERATIONS OF INCLUSIONS HERE
    for (unsigned int f = 0; f < CellNumbs.at(cell_type); ++f)
        if (Face_inclusion_index.at(f) > 9.0) agglomeration_SVector.at(f) = 1;

    agglomeration_fraction = std::count(agglomeration_SVector.begin(), agglomeration_SVector.end(), 1) /
                             CellNumbs.at(cell_type); // fraction of all grain boundaries containing agglomerations

    for (unsigned int f = 0; f < CellNumbs.at(cell_type); ++f) // agglomerations_set
        if (agglomeration_SVector.at(f) == 1) agglomerations_set.push_back(f);

    for (unsigned int af: agglomerations_set) /// lower energies of the faces containing agglomerations of inclusions
        Cell_energy.at(af) = sface_energy_aggl;

    /// Adhesion energy including stress concentrators from the side of inclusions
    for (unsigned int f = 0; f < CellNumbs.at(cell_type); ++f)
        Cell_energy.at(f) -= lsc_inclusion_energy.at(f);

    /// Crack face sequence creation
    /// (1) FAST case : cracks do NOT affect the process of the further fracture`
    //     for (unsigned int i = 0; i < CellNumbs.at(2); ++i){ ;}

    /// (2) SLOW case : cracks affect the process of the further fracture`
    double NewCrackNumb = 0;
//    double fractured_cells_fraction = 0; // pre-initial condition - no cracks
/// Possible initial pre-fractured state --> S_cVector
    if (Configuration_cState.at(cell_type).size() > 0)
        S_cVector = Configuration_cState.at(cell_type); // initial predefined system, if exists
// pre-initial 0-fractured state
    for (unsigned int lit = 0; lit < CellNumbs.at(cell_type); lit++)
        NotFracturedCellNumbs[lit] = lit; // Then the vector with the sequence of integers 1,2,3,... #Faces

    for (auto istr = S_cVector.begin(); istr != S_cVector.end(); ++istr)
        if (*istr != 0) NotFracturedCellNumbs.erase(
                    std::find(NotFracturedCellNumbs.begin(), NotFracturedCellNumbs.end(),
                              distance(S_cVector.begin(), istr)));
///        for (auto ictr = S_cVector.begin(); ictr != S_cVector.end(); ++ictr)
///         if(*ictr != 0) NotFracturedCellNumbs.erase(NotFracturedCellNumbs.begin() + distance(S_cVector.begin(),ictr)); // !!! Delete its element from the vector decreasing its size BUT

// calculation of the total max cracked cell fraction
    double total_max_cCells_fraction = 0.0;
    for (int j = 0; j < max_cfractions_vectors[cell_type].size(); ++j)
        if (max_cfractions_vectors[cell_type][j] > 0)
            total_max_cCells_fraction += max_cfractions_vectors[cell_type][j];

    if (total_max_cCells_fraction > 1.0) {
        cout << "WARNING! [Processing_Random()]: "s << cell_type << " total_max_sCell_fraction of " << cell_type
             << "-cells in the processing.ini file = " << total_max_cCells_fraction
             << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
        Out_logfile_stream << "WARNING! [Processing_Random()]: "s << cell_type << " total_max_sCell_fraction of "
                           << cell_type << "-cells in the processing.ini file = " << total_max_cCells_fraction
                           << " that is GREATER than 1 (!) Please decrease the fractions accordingly." << endl;
    } else if (total_max_cCells_fraction == 0.0) return crack_faces_sequence;

// initial calculated fractions of CRACKED faces
    double undamaged_cells_fraction = NotFracturedCellNumbs.size() / (double) CellNumbs.at(cell_type);
    double crack_cells_fraction = 1.0 - undamaged_cells_fraction; // special face vecror definition based on the ordinary face vector

/////// START OF THE DO CYCLE OF KINEMATIC FRACTURE ////////
    do { /// fracture loop over all fractured faces before f < c_max value taken from processing.ini file

/// Cracks' local elastic stress concentrations Energies (updeting at each simulation step)
        for (auto ks: s_faces_sequence) { // with inclusion (special faces)
//            lsc_inclusion_energy.at(ks) = nrr * (Face_inclusion_index.at(ks) / (0.75*Face_weight.at(ks)));
            lsc_crack_energy.at(ks) = ncr * (Face_crack_index.at(ks) / (0.75 * Face_weight.at(ks)));
        }
        for (auto os: o_faces_set) { // without inclusions (ordinary faces)
//            lsc_inclusion_energy.at(os) = nrm * (Face_inclusion_index.at(os) / Face_weight.at(os));
            lsc_crack_energy.at(os) = ncm * (Face_crack_index.at(os) / Face_weight.at(os));
        }

        /// Adhesion energy including stress concentrators from the side of inclusions
        // Face_energy.at(f) - constant during the frature process energy of a cell + inclusions effect
        for (unsigned int f = 0; f < CellNumbs.at(cell_type); ++f) // updating current cell's local elastic energy
            Cell_current_energy.at(f) = Cell_energy.at(f) - lsc_crack_energy.at(f); // [J/m^2]

        /// Number of the next fractured cell (!)
        NewCrackNumb = std::min_element(std::begin(Cell_current_energy), std::end(Cell_current_energy)) -
                       std::begin(Cell_current_energy); // gives index of the min element
// REPAIR        cout << "NewCrackNumb:  " << NewCrackNumb << "   NotFracturedCellNumbs.size():  " << NotFracturedCellNumbs.size() << endl;

        if (std::find(crack_faces_sequence.begin(), crack_faces_sequence.end(), NewCrackNumb) == crack_faces_sequence.end()) {
            crack_faces_sequence.push_back(NewCrackNumb); // if the element is NOT in the crack faces sequence (already fractured) - add it
            S_cVector.at(NewCrackNumb) = 1;
            NotFracturedCellNumbs.erase(std::find(NotFracturedCellNumbs.begin(), NotFracturedCellNumbs.end(),
                                                  NewCrackNumb)); // !!! Delete its element from the vector decreasing its size BUT
            Cell_energy.at(NewCrackNumb) = 100000.0; /// arbitrarily large value >> 1 !!!

        } // end if (std::find(...))

        /// Recalculation of the new crack index (CL) ///
/*
        //Recalculation of the NEW TJs types
        if (dim == 3) TJsTypes = EdgesTypesCalc(CellNumbs, s_faces_sequence, FES);
        else if (dim == 2) TJsTypes = NodesTypesCalc(CellNumbs, s_faces_sequence, FES);
        // GB_indices calculation
        for (unsigned int fn = 0; fn < CellNumbs.at(cell_type); ++fn) {
            vector<double> j_types_neigh_fractions = GBIndex(fn, FES,
                                                             TJsTypes); //Types (up to 100 kinds) of the edges incident to the considered Face
            // inclusion index calculation
            Face_inclusion_index.at(fn) = (j_types_neigh_fractions.at(1) + 2.0 * j_types_neigh_fractions.at(2) +
                                           3.0 * j_types_neigh_fractions.at(3));
        } // end of for ( fn < CellNumbs.at(cell_type))
//////// ????
*/
        newCrack_neigh_TJs.clear();
        if (dim == 3) TJsCrackTypes = EdgesTypesCalc(CellNumbs,crack_faces_sequence, FES);
        else if (dim == 2) TJsCrackTypes = NodesTypesCalc(CellNumbs, crack_faces_sequence, FES);

        newCrack_neigh_Faces.clear();
        for (int m = 0; m < CellNumbs.at(cell_type); ++m) // Loop over all the Faces
            if (AFS.coeff(NewCrackNumb, m) == 1 && S_cVector.at(m) == 0)
                newCrack_neigh_Faces.push_back(m);

        if (newCrack_neigh_Faces.size() > 0) {
            for (auto fn: newCrack_neigh_Faces) {
                vector<double> j_types_crack_n_fractions = GBIndex(fn, FES,TJsCrackTypes); //Types (up to 100 kinds) of the edges incident to the considered Face

                /// New crack index
                Face_crack_index.at(fn) = (j_types_crack_n_fractions.at(1) + 2.0 * j_types_crack_n_fractions.at(2) +
                                           3.0 * j_types_crack_n_fractions.at(3));
            } // end for (auto fn: new_crack_neigh_Faces)
        } // end if (new_crack_neigh_Faces.size() > 0)

        // current fractions of cracked cells
             undamaged_cells_fraction = NotFracturedCellNumbs.size() / (double) CellNumbs.at(cell_type);
             crack_cells_fraction = 1.0 - undamaged_cells_fraction; // special face vecror definition based on the ordinary face vector

        // if(10000*std::ceil(crack_cells_fraction) % 2 == 1)
        cout << "fractured cells fraction:   " << crack_cells_fraction << " < total_max_cCells_fraction "
             << total_max_cCells_fraction << endl;

//    } while (crack_cells_fraction < total_max_cCells_fraction); /// End of the PCC_Kinematic_cracking
    } while(crack_cells_fraction < total_max_cCells_fraction);

/*

    for (unsigned int i = 0; i < CellNumbs.at(cell_type); ++i) { // crack loop over all the GBs
        /// New surface energy
        // Local stress concentrators and related energies
        for (unsigned int k = 0; k < CellNumbs.at(2); ++k)
            lsc_crack_energy.at(k) = lsc_crack*(Face_crack_index.at(k)/Face_weight.at(k));
        // Current surface energy
        for (unsigned int l = 0; l < CellNumbs.at(2); ++l)
            Cell_current_energy.at(l) = Cell_energy.at(l) - lsc_crack_energy.at(l);

/// Number of the next fractured cell
        NewCrackNumb = std::min_element(std::begin(Cell_energy), std::end(Cell_energy)) - std::begin(Cell_energy); // gives index of the max element
        if (find(crack_faces_sequence.begin(), crack_faces_sequence.end(), NewCrackNumb) == crack_faces_sequence.end()) {
            crack_faces_sequence.push_back(NewCrackNumb); // if the element is not in the sequence - add
            S_crackVector.at(NewCrackNumb) = 1;
            Cell_energy.at(NewCrackNumb) = 100000.0; /// arbitrarily large value !!!
        } // end if

/// Recalculation of the new crack index (CL)
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

        //REPAIR cout << crack_fraction << endl;
    } // end for ( i < CellNumbs.at(2) )
*/
    crack_fraction = crack_faces_sequence.size() / (double) CellNumbs.at(cell_type);

    return crack_faces_sequence;
} /// end of Kinematic_cracking

/// -----------------------------------------------------------------------------------------------------------------///
/*
/// #2# Kinetic function for multiple cracking
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
*/