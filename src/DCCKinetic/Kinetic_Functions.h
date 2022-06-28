///================================ A part of the DCC Kinetic module =============================================================///
///=================================================================================================================================///
/** The library contains random and non-random functions for choice and generation of new identifications for k-Cells               *
 *  in the DCC Kinetic module. It makes them "special" and takes out of the set of "ordinary" k-Cells.                            **/
///================================================================================================================================///

int DCC_Kinetic_Wear(double ShearStress, vector<Tup> &Grain_Orientations, Eigen::SparseMatrix<double> const& FES, std::vector<unsigned int> &CellNumbs, char* input_folder, char* output_dir) {
/// Vectors
    /// Normals - reading from file
    string Norm_path = input_folder + "Normal.txt"s;
    char* normals = const_cast<char*>(Norm_path.c_str());
    /// Calculation of the NormalForce for each Face
    std::vector<Tup> NormalForce;
    std::vector<Tup> Norms_tuple = TuplesReader(normals); // Reading from file Normal vectors of all Faces
    for (unsigned long ik = 0; ik < CellNumbs.at(2); ++ik) NormalForce.push_back( { ShearStress*get<0>(Norms_tuple[ik]), ShearStress*get<1>(Norms_tuple[ik]),ShearStress*get<2>(Norms_tuple[ik]) } );
    /// Set grain orientations (Random + From file)
    // I. Random case
    // The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)
    srand((unsigned) time(NULL));

    cout << "Hello there!" << endl;
    return 1;
}

///*============================================================================================*///
///*============================== DCC_Kinetic_Plasticity function =============================*///

vector<Tup> DCC_Kinetic_Plasticity(Eigen::SparseMatrix<double> const& FES, std::vector<unsigned int> &CellNumbs, char* input_folder, char* output_dir)  {
    //resultant tuple
    vector<Tup> fraction_stress_temperature;
/// Functions
    vector<unsigned int> Metropolis(double &ShearStress, double &Temperature, std::vector<unsigned int> &CellNumbs, long iteration_number, double s, double E0, double alpha);
/// Variables
    vector<unsigned int> State_Vector(CellNumbs.at(2),0); // State vector filling with zeros

    /// Iteration over ShearStresses :: Input parameters
    double ShearStress = 0.0;
    double ShearStress_min = 1.0*pow(10,6), ShearStress_max = 1.0*pow(10,9);
    int simulation_steps = 100000, dSS = 0;
    /// Iteration over Temperature :: Input parameters
    double Temperature = 300.0;
    /// Other parameters
    double s = 0.1; // nano-slip value
    double lambda = 0.0*46.0*pow(10,9), //Shear modulus
             alpha = 7.0*pow(10,7); // model coefficient alpha*s^2

    ShearStress = ShearStress_min;
    dSS = (ShearStress_max - ShearStress_min)/ (double) simulation_steps;
    for (long i = 0; i < simulation_steps; ++i) {
        ShearStress += dSS;
        /// METROPOLIS algorithm
        State_Vector.clear();
        /// Iteration number for Metropolis algorithm
        long iteration_number = 3.0 * CellNumbs.at(2);

        /// =========== Metropolis algorithm ============>>>
        State_Vector = Metropolis(ShearStress, Temperature, CellNumbs, iteration_number, s, alpha, lambda); // Metropolis() is defined below

    /// Analysis
    double slip_fraction = std::count(State_Vector.begin(), State_Vector.end(), 1)/ (double) CellNumbs.at(2);
    fraction_stress_temperature.push_back(make_tuple(slip_fraction, ShearStress, Temperature));
    } // end of for (< calculation_steps)

    return fraction_stress_temperature;
}

vector<unsigned int> Metropolis(double &ShearStress, double &Temperature, std::vector<unsigned int> &CellNumbs, long iteration_number, double s, double alpha, double lambda){
    /// I. Constant initial state initialisation
    //zero vectors required for Processing_Random() input
    std::vector <unsigned int> SlipState_Vector(CellNumbs.at(2),0), s_faces_sequence(CellNumbs.at(2),0);
    double Rc = 8.31; //gas constant

    // ================> Initial p = 0.5 Face seeds
    Processing_Random(SlipState_Vector, s_faces_sequence, 0.5, 1, CellNumbs);

    // 1. Loop over all the Faces
    srand((unsigned) time(NULL)); // The function initialize random seed from the computer time (MUST BE BEFORE THE FOR LOOP!)

    for (long i = 0; i < iteration_number; ++i) {
    // 2. Random choice of a new Face
    long NewSlipNumber = rand() % (CellNumbs.at(2) - 1) ; // Random generation of the boundary number in the range from 0 to CellNumbs.at(2)
    // 3. If dH < 0 energetically favourable -> Accept trial and change the type
        if (SlipState_Vector.at(NewSlipNumber) == 1) {
            SlipState_Vector.at(NewSlipNumber) = 0;
        } //if
     // 4. Else if dH > 0 consider the acceptance probability
        else if (SlipState_Vector.at(NewSlipNumber) == 0) {
            //double slips_number = std::count(SlipState_Vector.begin(), SlipState_Vector.end(), 1);
            double P_ac = exp(- (alpha * pow(s,2.0) - ShearStress * s - lambda*s)/(Rc * Temperature));
            if (P_ac > 1) P_ac = 1.0;
                //cout << "P_ac =\t" << alpha*pow(s,2.0) - ShearStress*s << "\trandom number\t" << (alpha*pow(s,2.0) - ShearStress*s - lambda*s)/(Rc * Temperature) << endl;

            double rv = (rand() / (RAND_MAX + 1.0)); // Generate random value in the range [0,1]
                if (rv <= P_ac) {
                //     if (P_ac > 0) cout << "P_ac =\t" << (alpha*pow(s,2.0) - ShearStress*s - lambda*s)/(Rc * Temperature) << "\trandom number\t" << rv << endl;
                    SlipState_Vector.at(NewSlipNumber) = 1;
                } else continue;
        }

    } // for loop (i < iteration_number)

    return SlipState_Vector;
} // end of Metropolis