#ifndef PCC_PROCESSING_DESIGN_ANALYTICAL_SOLUTIONS_H
#define PCC_PROCESSING_DESIGN_ANALYTICAL_SOLUTIONS_H

/// #include "../../PCC_Processing/functions/Processing_Assignment_functions.h"


std::vector<vector<double>> TJsAnalytics_random(unsigned int HGBnumber, std::vector<tuple<double, double>> &ConfEntropies) {
    std::vector<vector<double>> j_random_vector; // main function output {p, jj0, jj1, jj2, jj3}
    //    vector<tuple<double, double>> ConfEntropies: a vector of tuples containing [0] -> Face_Entropy_Median, [1] -> Face_Entropy_Skrew

    double P00=0, P11=0, P22=0, P10=0, P20=0, P21=0, p=0, q=0, J0=0, J1=0, J2=0, J3=0, Jall=0;
    double dp = 1.0/ (double) HGBnumber;

    for (unsigned int k = 0; k < HGBnumber; ++k) {
        P00 = p; P11 = p; P22 = p; P10 = p; P20 = p; P21 = p; P22 = p;

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;

        // Conversion from numbers to fractions
        double  jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;
        if (j0s != 0) j0s = jj0* log2(jj0); if (j1s != 0) j1s = jj1* log2(jj1); if (j2s != 0) j2s = jj2* log2(jj2); if (j3s != 0) j3s = jj3* log2(jj3); //Gives 0 in entropy!
        /// Configuration Entropy related with Faces
        // Median part in the entropy decomposition
        double Configurational_Face_Entropy = 0.0, Face_Entropy_Median = 0.0, Face_Entropy_Skrew = 0.0;

//        Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

        std::vector<double> j_types_fractions = {jj0, jj1, jj2, jj3};
        /// add to the output vector
        j_random_vector.push_back({p, jj0, jj1, jj2, jj3});

        if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
            Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(jj0 * jj1 * jj2 * jj3);
        } else Face_Entropy_Median = 0.0;

        /// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
        for (int j = 0; j < j_types_fractions.size(); j++)
            for (int i = 0; i < j; i++)
                if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
                    Face_Entropy_Skrew +=
                            -(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
                            log2(j_types_fractions[i] / j_types_fractions[j]);
                } else Face_Entropy_Skrew += 0.0;

        /// add to the output
        ConfEntropies.push_back(make_tuple(Face_Entropy_Median, Face_Entropy_Skrew));

        /// Increment of p fraction
        p += dp;

    } // end of for ( k < HGBnumb )

    return j_random_vector;
} // END of TJsAnalytics_random()

std::vector<vector<double>> TJsAnalytics_crystallography(unsigned int HGBnumber, std::vector<tuple<double, double>> &ConfEntropies) {
    std::vector<vector<double>> j_crystallography_vector; // main function output {p, jj0, jj1, jj2, jj3}
    //    vector<tuple<double, double>> ConfEntropies: a vector of tuples containing [0] -> Face_Entropy_Median, [1] -> Face_Entropy_Skrew

    double P00=0, P11=0, P22=0, P10=0, P20=0, P21=0, p=0, q=0, J0=0, J1=0, J2=0, J3=0, Jall=0;
    double dp = 1.0/ (double) HGBnumber;

    for (unsigned int k = 0; k < HGBnumber; ++k) {
        q = 1.0 - p;

        P00 = p;
        if (p <= 0.75) {
            P10 = (1.0 - 6.0 * pow(q, (1.0 / 2.0)) + 15.0 * q - 10.0 * pow(q, (3.0 / 2.0))) / (3.0 * q);
            P20 = (2.0 - 12.0 * pow(q, (1.0 / 2.0)) + 24.0 * q - 14.0 * pow(q, (3.0 / 2.0))) /
                  (-1.0 + 6.0 * pow(q, (1.0 / 2.0)) - 12.0 * q + 10.0 * pow(q, (3.0 / 2.0)));
            P11 = (2.0 + 8.0 * pow(q, (1.0 / 2.0)) - 10.0 * q) / (3.0 + 3.0 * pow(q, (1.0 / 2.0)));
            P21 = (1.0 - 5.0 * pow(q, (1.0 / 2.0)) + 4.0 * q) / (-1.0 + 5.0 * pow(q, (1.0 / 2.0)) - 10.0 * q);
            P22 = (3.0 + 6.0 * pow(q, (1.0 / 2.0))) / (2.0 + 10.0 * pow(q, (1.0 / 2.0)));
        }
        if (p > 0.75) {
            P10 = (3.0 - 2.0 * pow(q, (1.0 / 2.0))) / 3.0;
            P20 = 1.0;
            P11 = (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0))) / (3.0 - 3.0 * q);
            P21 = (3.0 - 4.0 * pow(q, (1.0 / 2.0))) / (3.0 - 2.0 * pow(q, (1.0 / 2.0)));
            P22 = (3.0 - 9.0 * q + 6.0 * pow(q, (3.0 / 2.0))) / (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0)));
        }

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;

        // Conversion from numbers to fractions

        /// Median part in the entropy decomposition
        double Configurational_Face_Entropy = 0.0, Face_Entropy_Median = 0.0, Face_Entropy_Skrew = 0.0;

// Conversion from numbers to fractions
        double jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;
        Configurational_Face_Entropy = 0.0;
        if (j0s != 0) j0s = jj0* log2(jj0); if (j1s != 0) j1s = jj1* log2(jj1); if (j2s != 0) j2s = jj2* log2(jj2); if (j3s != 0) j3s = jj3* log2(jj3); //Gives 0 in entropy!
        /// Configuration Entropy related with Faces
        Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);

        /// Median part in the entropy decomposition
        Face_Entropy_Median = 0.0, Face_Entropy_Skrew = 0.0;

        std::vector<double> j_types_fractions = {jj0, jj1, jj2, jj3}; /// using values with pow(10,-10) instead of 0s!
/// add to the output vector
        j_crystallography_vector.push_back({p, jj0, jj1, jj2, jj3});

        if (j0s!=0 && j1s!=0 && j2s!=0 && j3s!=0) {
            Face_Entropy_Median = -(1.0 / j_types_fractions.size()) * log2(jj0 * jj1 * jj2 * jj3);
        } else Face_Entropy_Median = 0.0;

        /// Screw part (divergence from the uniform distribution -> S_max) in the entropy decomposition
        for (int j = 0; j < j_types_fractions.size(); j++)
            for (int i = 0; i < j; i++)
                if (j_types_fractions[i]!=0 && j_types_fractions[j]!=0) {
                    Face_Entropy_Skrew +=
                            -(1.0 / j_types_fractions.size()) * (j_types_fractions[i] - j_types_fractions[j]) *
                            log2(j_types_fractions[i] / j_types_fractions[j]);
                } else Face_Entropy_Skrew += 0.0;

        /// add to the output
        ConfEntropies.push_back(make_tuple(Face_Entropy_Median, Face_Entropy_Skrew));

        /// Increment of p fraction
        p += dp;

    } // end of for ( k < HGBnumb )

    return j_crystallography_vector;
}

std::vector<vector<double>> TJDsAnalytics_random(unsigned int HGBnumber) {
    std::vector<vector<double>> d_random_vector; // main function output {p, dd1, dd2, dd3}

    double P00=0, P11=0, P22=0, P10=0, P20=0, P21=0, p=0, q=0, J0=0, J1=0, J2=0, J3=0, Jall=0;
    double dp = 1.0/ (double) HGBnumber;

    for (unsigned int k = 0; k < HGBnumber; ++k) {
        P00 = p; P11 = p; P22 = p; P10 = p; P20 = p; P21 = p; P22 = p;

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;

        // Conversion from numbers to fractions
        double  jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;

        std::vector<double> d_types_fractions;
        if (jj0 != 1.0) d_types_fractions = {jj0/(1.0 - jj0), jj1/(1.0 - jj0), jj2/(1.0 - jj0), jj3/(1.0 - jj0)};
            else d_types_fractions = {0, 0, 0, 0};

        /// add to the output vector
        d_random_vector.push_back({p, d_types_fractions.at(0), d_types_fractions.at(1),d_types_fractions.at(2)});

        /// Increment of p fraction
        p += dp;

    } // end of for ( k < HGBnumb )

    return d_random_vector;
}

std::vector<vector<double>> TJDsAnalytics_crystallography(unsigned int HGBnumber) {
    std::vector<vector<double>> d_crystallography_vector; // main function output

    double P00=0, P11=0, P22=0, P10=0, P20=0, P21=0, p=0, q=0, J0=0, J1=0, J2=0, J3=0, Jall=0;
    double dp = 1.0/ (double) HGBnumber;

    for (unsigned int k = 0; k < HGBnumber; ++k) {
        q = 1.0 - p;

        P00 = p;
        if (p <= 0.75) {
            P10 = (1.0 - 6.0 * pow(q, (1.0 / 2.0)) + 15.0 * q - 10.0 * pow(q, (3.0 / 2.0))) / (3.0 * q);
            P20 = (2.0 - 12.0 * pow(q, (1.0 / 2.0)) + 24.0 * q - 14.0 * pow(q, (3.0 / 2.0))) /
                  (-1.0 + 6.0 * pow(q, (1.0 / 2.0)) - 12.0 * q + 10.0 * pow(q, (3.0 / 2.0)));
            P11 = (2.0 + 8.0 * pow(q, (1.0 / 2.0)) - 10.0 * q) / (3.0 + 3.0 * pow(q, (1.0 / 2.0)));
            P21 = (1.0 - 5.0 * pow(q, (1.0 / 2.0)) + 4.0 * q) / (-1.0 + 5.0 * pow(q, (1.0 / 2.0)) - 10.0 * q);
            P22 = (3.0 + 6.0 * pow(q, (1.0 / 2.0))) / (2.0 + 10.0 * pow(q, (1.0 / 2.0)));
        }
        if (p > 0.75) {
            P10 = (3.0 - 2.0 * pow(q, (1.0 / 2.0))) / 3.0;
            P20 = 1.0;
            P11 = (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0))) / (3.0 - 3.0 * q);
            P21 = (3.0 - 4.0 * pow(q, (1.0 / 2.0))) / (3.0 - 2.0 * pow(q, (1.0 / 2.0)));
            P22 = (3.0 - 9.0 * q + 6.0 * pow(q, (3.0 / 2.0))) / (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0)));
        }

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;
// Conversion from numbers to fractions
        // Conversion from numbers to fractions
        double  jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;
        std::vector<double> d_types_fractions;
        if (jj0 != 1.0) d_types_fractions = {jj0/(1.0 - jj0), jj1/(1.0 - jj0), jj2/(1.0 - jj0), jj3/(1.0 - jj0)};
            else d_types_fractions = {0, 0, 0, 0};

    /// add to the output vector
        d_crystallography_vector.push_back({p, d_types_fractions.at(0), d_types_fractions.at(1),d_types_fractions.at(2)});

    /// Increment of p fraction
        p += dp;

    } // end of for ( k < HGBnumb )

    return d_crystallography_vector;
}

/// # ? # plot Mackenzie distribution for grain boundary disorientations
std::vector<double> Mackenzie_distribution(void) {
    std::vector<double> disorientations;

    SpMat AGS = SMatrixReader(paths.at(3 + (dim - 3)), (CellNumbs.at(3)), (CellNumbs.at(3))); //all Faces
    AGS = 0.5 * (AGS + SparseMatrix<double>(AGS.transpose())); // Full matrix instead of triagonal
/// ==============================================

    std::vector<vector<double>> grain_quaternions(CellNumbs.at(3 + (dim - 3)), std::vector<double>(4)),
            grain_triangle_quaternions(CellNumbs.at(3 + (dim - 3)), std::vector<double>(3)); // grain orientations 2-quaternions and 3-quternions
    double HAGBs_threshold = 15.0; // treshold for the definition of the "high-angle" disorientations

/// Initial TRIANGLE LATTICE
// q1 \in [0, 1], q2 \in [0, 1], q3 \in [0, 1]
// n_grains - number of points, dq - "lattice parameter"
    int n_lattice_points = 100; // is an arbitrary user-defined parameter here (!)
    double  dq = 1.0 / (double) n_lattice_points;

    std::vector<vector<double>> q_coord_vector; // grain rotation lattice
    for(int qi = 0; qi <= n_lattice_points; ++qi) // 1
        for(int qk = 0; qk <= (n_lattice_points - qi); ++qk) // 2
            q_coord_vector.push_back({qi*dq, qk*dq, 1.0 - qi*dq - qk*dq});
// REPAIR    for(auto gtq : q_coord_vector)
//        cout << "q_coord_vector " << gtq.at(0) << " " << gtq.at(1) << " " << gtq.at(2) << " "<< gtq.at(0) + gtq.at(1) + gtq.at(2) << endl;

/// Assigning of the INITIAL grain orientations for each grain in the PCC
    double q0_coord = 0;
    unsigned int new_grain_triangle_coords = 0;
    grain_triangle_quaternions.clear();
    grain_quaternions.clear();

    for (int i = 0; i < CellNumbs.at(3 + (dim - 3)); ++i) {

// random choice of a point in a triangle coordinate space
        new_grain_triangle_coords = NewCellNumb_R(q_coord_vector.size());

        std::vector<double> grain_q_quaternion(3), grain_full_quaternion(4);
/// triangle Q3-quaternions
        grain_q_quaternion.at(0) = q_coord_vector[new_grain_triangle_coords][0];
        grain_q_quaternion.at(1) = q_coord_vector[new_grain_triangle_coords][1];
        grain_q_quaternion.at(2) = q_coord_vector[new_grain_triangle_coords][2];

        grain_triangle_quaternions.push_back({grain_q_quaternion.at(0), grain_q_quaternion.at(1), grain_q_quaternion.at(2)});
// REPAIR        cout << "grain_triangle_quaternions " << grain_triangle_quaternions.back().at(0) << " " << grain_triangle_quaternions.back().at(1) << " " << grain_triangle_quaternions.back().at(2) << " "<< grain_triangle_quaternions.back().at(0) + grain_triangle_quaternions.back().at(1) + grain_triangle_quaternions.back().at(2) << endl;

/// full grain quaternions
        /// angle - an additional random generation
        q0_coord = rand() / (RAND_MAX + 1.0); //        q0_coord = std::sqrt(1.0 - pow(grain_full_quaternion.at(0),2) - pow(grain_full_quaternion.at(1),2) - pow(grain_full_quaternion.at(2),2));
        // axis
        grain_full_quaternion.at(0) = std::sqrt(grain_q_quaternion.at(0)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(1) = std::sqrt(grain_q_quaternion.at(1)*(1.0 - pow(q0_coord,2)));
        grain_full_quaternion.at(2) = std::sqrt(grain_q_quaternion.at(2)*(1.0 - pow(q0_coord,2)));

// full grain quaternions vector
        grain_quaternions.push_back({q0_coord, grain_full_quaternion.at(0), grain_full_quaternion.at(1), grain_full_quaternion.at(2)});
// REPAIR        cout << "full grain quaternions " << grain_quaternions.back().at(0) << " " << grain_quaternions.back().at(1) << " " << grain_quaternions.back().at(2) << " " << grain_quaternions.back().at(3) << "  " << pow(grain_quaternions.back().at(0),2) + pow(grain_quaternions.back().at(1),2) + pow(grain_quaternions.back().at(2),2) + pow(grain_quaternions.back().at(3),2)<< endl;
    } // end for (int i = 0; i < CellNumbs.at(3 + (dim - 3)); ++i)

// REPAIR   for(auto gtq : grain_triangle_quaternions)
//        cout << "triangle Q3-quaternions " << gtq.at(0) << " " << gtq.at(1) << " " << gtq.at(2) << " "<< gtq.at(0) + gtq.at(1) + gtq.at(2) << endl;

// REPAIR    for(auto gq : grain_quaternions)
//        cout << "full grain quaternions " << gq.at(0) << " " << gq.at(1) << " " << gq.at(2) << " " << gq.at(3) << "  " << pow(gq.at(0),2) + pow(gq.at(1),2) + pow(gq.at(2),2) + pow(gq.at(3),2)<< endl;

/// Mackenzie check
    std::vector<double> disorientations_vector;
    for(unsigned int gr = 0; gr < CellNumbs.at(3 - (dim -3)); ++gr)
        for(unsigned int grn = 0; grn < gr; ++grn)
            if(AGS.coeff(gr,grn) != 0) {
                //        cout << "coeff " << AGS.coeff(gr,grn) << endl;
                disorientations_vector.push_back(Get_2grains_FCCdisorientation(grain_quaternions.at(gr), grain_quaternions.at(grn)));
                cout << disorientations_vector.back()*57.3 << endl;

            }

    return disorientations;
} // END of Mackenzie_distribution(void)

/*
int TJsAnalytics_indices(double Smax, double Smin, double Srand) {
    double j0 = 0, j1 = 0, j2 = 0, j3 = 0, j0s = 0, j1s = 0, j2s = 0, j3s = 0;
    unsigned int HGBnumber = 100;     // HGBnumber^3 is the Number of calculation points
    double dj = 1.0 / (double) HGBnumber;
    char* odir = const_cast<char*>(output_dir.c_str()); // const_cast for output directory

    ofstream Out_ip, Out_ip_rand, Out_ip_text;
    Out_ip.open(odir + "ip_analytical.txt"s, ios::trunc);
    Out_ip << "(1)index_S" << "    " << "(2)Conf_Entropy" << "    " << "(3)j0" << "    " << "(4)j1" << "    " << "(5)j2" << "    " << "(6)j3" << endl;
    Out_ip.close();
    Out_ip_rand.open(odir + "ip_random_analytical.txt"s, ios::trunc);
    Out_ip_rand << "(1)index_S" << "    " << "(2)Conf_Entropy" << "    " << "(3)j0" << "    " << "(4)j1" << "    " << "(5)j2" << "    " << "(6)j3" << endl;
    Out_ip_rand.close();
    Out_ip_text.open(odir + "ip_texture_analytical.txt"s, ios::trunc);
    Out_ip_text << "(1)index_S" << "    " << "(2)Conf_Entropy" << "    " << "(3)j0" << "    " << "(4)j1" << "    " << "(5)j2" << "    " << "(6)j3" << endl;
    Out_ip_text.close();

    Out_ip.open(odir + "ip_analytical.txt"s, ios::app);
    /// Loop over all the TJs
    for (int i = 0; i < HGBnumber; ++i) {
        for (int k = 0; k < HGBnumber; ++k) {
            for (int l = 0; l < HGBnumber; ++l) {
                /// using values with pow(10,-10) instead of 0s!
                j0 = 1.0 - j1 - j2 - j3;
                j0s = 0, j1s = 0, j2s = 0, j3s = 0;
                if (j0 != 0) j0s = j0* log2(j0); if (j1 != 0) j1s = j1* log2(j1); if (j2 != 0) j2s = j2* log2(j2); if (j3 != 0) j3s = j3* log2(j3); //Gives 0 in entropy!

                /// Configuration Entropy related with Faces
                double Conf_Entropy = - (j0s + j1s + j2s + j3s);
                double index_S = 0;
                if (Conf_Entropy >= Srand && Smax != Srand) index_S = ( Conf_Entropy - Srand)/ ( Smax - Srand );
                else if (Conf_Entropy < Srand && Smin != Srand) index_S = ( Conf_Entropy - Srand)/ ( Srand - Smin);

                if (j0 >= 0.0) Out_ip << index_S << "    " << Conf_Entropy << "    " << j0 << "    " << j1 << "    " << j2 << "    " << j3 << endl;

                /// Output
                j1 += dj;
            }
            j1 = 0.0;
            j2 += dj;
        }
        j3 += dj;
        j2 = 0.0;
    }
    j3 = 0.0;

    Out_ip.close();

    /// Random case
    Out_ip_rand.open(odir + "ip_random_analytical.txt"s, ios::app);
    Out_ip_text.open(odir + "ip_texture_analytical.txt"s, ios::app);

    double dp = 1.0/ (double) HGBnumber;
    double P00=0, P11=0, P22=0, P10=0, P20=0, P21=0, p=0, q=0, J0=0, J1=0, J2=0, J3=0, Jall=0;
    for (unsigned int k = 0; k < HGBnumber; ++k) {
        //1 случай: случайное распределение без учета кристаллографии
        P00 = p; P11 = p; P22 = p; P10 = p; P20 = p; P21 = p; P22 = p;

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;

        // Conversion from numbers to fractions
        double  jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;
        double Configurational_Face_Entropy = 0.0;
        if (j0s != 0) j0s = jj0* log2(jj0); if (j1s != 0) j1s = jj1* log2(jj1); if (j2s != 0) j2s = jj2* log2(jj2); if (j3s != 0) j3s = jj3* log2(jj3); //Gives 0 in entropy!
        /// Configuration Entropy related with Faces
        Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);
        double index_S = 0;
        if (Configurational_Face_Entropy >= Srand && Smax != Srand) index_S = ( Configurational_Face_Entropy - Srand)/ ( Smax - Srand );
        else if (Configurational_Face_Entropy < Srand && Smin != Srand) index_S = ( Configurational_Face_Entropy - Srand)/ ( Srand - Smin);

        Out_ip_rand << index_S << "    " << Configurational_Face_Entropy << "    " << jj0 << "    " << jj1 << "    " << jj2 << "    " << jj3 << endl;
        //2 случай: с учетом ограничений кристаллографии
        q = 1.0 - p;

        P00 = p;
        if (p <= 0.75) {
            P10 = (1.0 - 6.0 * pow(q, (1.0 / 2.0)) + 15.0 * q - 10.0 * pow(q, (3.0 / 2.0))) / (3.0 * q);
            P20 = (2.0 - 12.0 * pow(q, (1.0 / 2.0)) + 24.0 * q - 14.0 * pow(q, (3.0 / 2.0))) /
                  (-1.0 + 6.0 * pow(q, (1.0 / 2.0)) - 12.0 * q + 10.0 * pow(q, (3.0 / 2.0)));
            P11 = (2.0 + 8.0 * pow(q, (1.0 / 2.0)) - 10.0 * q) / (3.0 + 3.0 * pow(q, (1.0 / 2.0)));
            P21 = (1.0 - 5.0 * pow(q, (1.0 / 2.0)) + 4.0 * q) / (-1.0 + 5.0 * pow(q, (1.0 / 2.0)) - 10.0 * q);
            P22 = (3.0 + 6.0 * pow(q, (1.0 / 2.0))) / (2.0 + 10.0 * pow(q, (1.0 / 2.0)));
        }
        if (p > 0.75) {
            P10 = (3.0 - 2.0 * pow(q, (1.0 / 2.0))) / 3.0;
            P20 = 1.0;
            P11 = (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0))) / (3.0 - 3.0 * q);
            P21 = (3.0 - 4.0 * pow(q, (1.0 / 2.0))) / (3.0 - 2.0 * pow(q, (1.0 / 2.0)));
            P22 = (3.0 - 9.0 * q + 6.0 * pow(q, (3.0 / 2.0))) / (3.0 - 6.0 * q + 2.0 * pow(q, (3.0 / 2.0)));
        }

        J0 = (1.0 - P00) * (1.0 - P10) * (1.0 - P20);
        J1 = P00 * (1.0 - P11) * (1.0 - P21) + (1.0 - P00) * P10 * (1.0 - P21) + (1.0 - P00) * (1.0 - P10) * P20;
        J2 = P00 * P11 * (1.0 - P22) + P00 * (1.0 - P11) * P21 + (1.0 - P00) * P10 * P21;
        J3 = P00 * P11 * P22;
        Jall = J0 + J1 + J2 + J3;
// Conversion from numbers to fractions
        jj0 = J0/Jall, jj1 = J1/Jall, jj2 = J2/Jall, jj3 = J3/Jall, j0s = jj0, j1s = jj1, j2s = jj2, j3s = jj3;
        Configurational_Face_Entropy = 0.0;
        if (j0s != 0) j0s = jj0* log2(jj0); if (j1s != 0) j1s = jj1* log2(jj1); if (j2s != 0) j2s = jj2* log2(jj2); if (j3s != 0) j3s = jj3* log2(jj3); //Gives 0 in entropy!
        /// Configuration Entropy related with Faces
        Configurational_Face_Entropy = - (j0s + j1s + j2s + j3s);
        /// index i(S)
        index_S = 0;
        if (Configurational_Face_Entropy >= Srand && Smax != Srand) index_S = ( Configurational_Face_Entropy - Srand)/ ( Smax - Srand );
        else if (Configurational_Face_Entropy < Srand && Smin != Srand) index_S = ( Configurational_Face_Entropy - Srand)/ ( Srand - Smin);

        Out_ip_text << index_S << "    " << Configurational_Face_Entropy << "    " << jj0 << "    " << jj1 << "    " << jj2 << "    " << jj3 << endl;

        /// Increment of p fraction
        p += dp;
    } // end of for ( k < HGBnumb )

    Out_ip_rand.close();
    Out_ip_text.close();

    return 0;
} /// End of TJsAnalytics_indices() function
*/

#endif //PCC_PROCESSING_DESIGN_ANALYTICAL_SOLUTIONS_H
