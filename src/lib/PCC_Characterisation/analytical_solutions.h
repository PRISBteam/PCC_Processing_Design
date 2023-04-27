#ifndef PCC_PROCESSING_DESIGN_ANALYTICAL_SOLUTIONS_H
#define PCC_PROCESSING_DESIGN_ANALYTICAL_SOLUTIONS_H

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
