#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <set>
#include <random> // Require C++ 11 and above
#include <cmath>
#include <algorithm>

using namespace std; // standard/STL namespace

std::tuple<int,int,int,int,int,int> probe_set;
std::vector<std::tuple<double, double, double, double, double, double>> ps_vector, fs_vector;

///unsigned int number_of_fibres = 1; // set here the number of fibres
///unsigned int NewFibre_Rand(unsigned int n_fibres); // return random number of a fibre
///double NewCoord_R(double coord_range_a, double coord_range_b); //return coordinates in [a,b] interval

///Output folder
std::string output2_dir = "../"s;
//extern
char* outputdir = const_cast<char*>(output2_dir.c_str()); // 'source_path' is a path to the directory reading from the 'config/main.ini' file
std::ofstream Out_logfile2_stream, Out_numbers2_stream, Out_egaps2_stream; // 'Processing_Design.log' file output of the entire computation process as a copy of the console output

void energy_plasticity() {
/// Global log file output:
//    Out_logfile_stream.open(output_dir + "TF_modeller_new.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Out_logfile2_stream.open(output2_dir + "TF_modeller_NEW.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Out_numbers2_stream.open(output2_dir + "TF_energy_level_numbers_NEW.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Out_egaps2_stream.open(output2_dir + "TF_energy_gaps_NEW.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function

    int increment = 5;
    for(int s11 = 1; s11 < 120; s11+=increment) {
        for(int s22 = 1; s22 < 120; s22+=increment) {
            for(int s33 = 1; s33 < 120; s33+=increment) {
                for(int s12 = 1; s12 < 120; s12+=increment) {
                    for(int s13 = 1; s13 < 120; s13+=increment) {
                        for(int s23 = 1; s23 < 120; s23+=increment) {
                            get<0>(probe_set) = (double) s11; get<1>(probe_set) = (double) s22; get<2>(probe_set) = (double) s33; get<3>(probe_set) = (double) s12; get<4>(probe_set) = (double) s13; get<5>(probe_set) = (double) s23;
                            ps_vector.push_back(probe_set);
                        }}} }}}
    cout << "ps_vector size = " << ps_vector.size();

    double pressure, Mizes_stress, t_factor, tau_octo;
    std::vector<vector<double>> factors;
    double YS = 100;
    for(auto psvar : ps_vector) {
        pressure = (get<0>(psvar) + get<1>(psvar) + get<2>(psvar))/3.0;
        Mizes_stress = std::sqrt(pow((get<0>(psvar)-get<1>(psvar)),2) + pow((get<1>(psvar)-get<2>(psvar)),2) + pow((get<2>(psvar)-get<0>(psvar)),2)+6.0*(pow(get<3>(psvar),2)+pow(get<4>(psvar),2)+pow(get<5>(psvar),2)))/std::sqrt(2.0);
        t_factor = pressure/Mizes_stress;

//        tau_octo =  std::sqrt(pow((get<0>(psvar)-get<1>(psvar)),2) + pow((get<1>(psvar)-get<2>(psvar)),2) + pow((get<2>(psvar)-get<0>(psvar)),2))/3.0;
        tau_octo =  std::sqrt(pow((get<0>(psvar)-get<2>(psvar)),2));

        if (abs(Mizes_stress) > YS && abs(Mizes_stress) < 120.0) //            fs_vector.push_back(psvar);
            factors.push_back({get<0>(psvar), get<1>(psvar), get<2>(psvar), get<3>(psvar), get<4>(psvar), get<5>(psvar), pressure, Mizes_stress, t_factor});
    }

    double alpha_t = 6.49*pow(10.0,(-25.0)); //8.96//6.4904*pow(10.0,(-25.0)); // get<0>(probe_set) = (double) s11; get<1>(probe_set) = (double) s22; get<2>(probe_set) = (double) s33; get<3>(probe_set) = (double) s12; get<4>(probe_set) = (double) s13; get<5>(probe_set) = (double) s23;
    std::set<double> slip_energies, energy_gaps;
    std::vector<set<double>> energies_vector, egaps_vector;
    double A2=0, A3=0, A6=0, B2=0, B4=0, B6=0, C1=0, C3=0, C5=0, D1=0, D4=0, D6=0;

/// Out_logfile_stream << "Pressure(MPa)" << " " << "Mizes(MPa)" << " " << "t_factor" << " " << " alpha_t " << " " << " A2 " << " " << " A3 " << " "  << " A6 " << " " << " B2 " << " " << " B4 " << " "  << " B5 " << " " << " C1 " << " " << " C3 " << " "  << " C5 " << " " << " D1 " << " " << " D4 " << " " << " D6 " << endl;
    Out_logfile2_stream << "Pressure(MPa)" << " " << "Mizes(MPa)" << " " << "t_factor" << " " << " alpha_t " << " " << " energy " << " " << " id " << endl;
    Out_numbers2_stream << "Pressure(MPa)" << " " << "Mizes(MPa)" << " " << "t_factor" << " " << " alpha_t " << " " << " number_of_levels " << "e_1" << endl;
    Out_egaps2_stream << "Pressure(MPa)" << " " << "Mizes(MPa)" << " " << "t_factor" << " " << " alpha_t " << " " << " number_of_gaps " << " " << "egaps " << endl;

    double ny2 = -0.577350269, ny3 = 0.577350269, ny4 = 0.577350269; // normals
    double nz2 = 0.577350269, nz3 = 0.577350269, nz4 = 0.577350269;
    double naa2 = -0.577350269, naa3 = -0.577350269, naa4 = 0.577350269;
    double nab2 = 0.577350269, nab3 = -0.577350269, nab4 = 0.577350269;
    double ae1 = 4.42E-15, ae2 = 2.56E-10, ae3 = 4.56E-23, ae4 = 74666662.54;
    //          2-cell area   |     BurgV   |  Support Vol |  inn.prd. Coef.
    double uy7 = 0.0, uy8 = 0.707106781, uy9 = 0.707106781; // unit vectors
    double uz7 = 0.0, uz8 = -0.707106781, uz9 = 0.707106781;
    double uaa7 = 0.707106781, uaa8 = 0.0, uaa9 = 0.707106781;
    double uab7 = -0.707106781, uab8 = 0.0, uab9 = 0.707106781;
    double uac7 = -0.707106781, uac8 = 0.707106781, uac9 = 0.0;
    double uad7 = 0.707106781, uad8 = 0.707106781, uad9 = 0.0;

    for(auto f_out : factors) { // f_out.at(0) -> s11 // f_out.at(1) -> s22 // f_out.at(2) -> s33 // f_out.at(3) -> s12 // f_out.at(4) -> s13 // f_out.at(5) -> s23 //
// f_out.at(6) ->  pressure // f_out.at(7) -> Mizes_stress // // f_out.at(8) ->  t_factor //

//    A2 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(0) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(-4.08212913E-16) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(4.08212913E-16))/alpha_t;
        A2 = ((f_out.at(0)*ny2 + f_out.at(3)*ny3 + f_out.at(4)*ny4)*uz7 + (f_out.at(3)*ny2 + f_out.at(1)*ny3 + f_out.at(5)*ny4)*uz8 + (f_out.at(4)*ny2 + f_out.at(5)*ny3+f_out.at(2)*ny4)*uz9)*ae1*(pow(ae1,2)*ae2*ae4/ae3) / alpha_t;
        if (A2 > 1.E-10) {
            slip_energies.insert(A2);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << A2 << " A2 " << endl;
        }
//    get<0>(probe_set) = (double) s11; get<1>(probe_set) = (double) s22; get<2>(probe_set) = (double) s33; get<3>(probe_set) = (double) s12; get<4>(probe_set) = (double) s13; get<5>(probe_set) = (double) s23;
///    factors.push_back(s11, s22, s33, s12, s13, s23, pressure, Mizes_stress, tau_octo);
///    factors.push_back( 0,   1,   2,   3,   4,   5,      6,         7,            8);

// A3 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(4.08212913E-16) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(0) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(4.08212913E-16))/alpha_t;
        A3 = ((f_out.at(0)*ny2 + f_out.at(3)*ny3 + f_out.at(4)*ny4)*uaa7 + (f_out.at(3)*ny2 + f_out.at(1)*ny3 + f_out.at(5)*ny4)*uaa8 + (f_out.at(4)*ny2 + f_out.at(5)*ny3+f_out.at(2)*nz4)*uaa9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        if (A3 > 1.E-10) {
            slip_energies.insert(A3);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << A3 << " A3 " << endl;
        }
// A6 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(4.08212913E-16) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(0))/alpha_t;
        A6 = ((f_out.at(0)*ny2 + f_out.at(3)*ny3 + f_out.at(4)*ny4)*uad7 + (f_out.at(3)*ny2 + f_out.at(1)*ny3+f_out.at(5)*ny4)*uad8 + (f_out.at(4)*ny2 + f_out.at(5)*ny3+f_out.at(2)*ny4)*uad9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        if (A6 > 1.E-10) {
            slip_energies.insert(A6);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << A6 << " A6 " << endl;
        }
        //B2 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(0) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(-4.08212913E-16) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(4.08212913E-16))/alpha_t;
        B2 = ((f_out.at(0)*nz2 + f_out.at(3)*nz3+f_out.at(4)*nz4)*uz7 + (f_out.at(3)*nz2 + f_out.at(1)*nz3+f_out.at(5)*nz4)*uz8+(f_out.at(4)*nz2 + f_out.at(5)*nz3+f_out.at(2)*nz4)*uz9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        if (B2 > 1.E-10) {
            slip_energies.insert(B2);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << B2 << " B2 " << endl;
        }
        //B4 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(-4.08212913E-16) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(0) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(4.08212913E-16))/alpha_t;
        B4 = ((f_out.at(0)*nz2 + f_out.at(3)*nz3+f_out.at(4)*nz4)*uab7 + (f_out.at(3)*nz2 + f_out.at(1)*nz3+f_out.at(5)*nz4)*uab8+(f_out.at(4)*nz2 + f_out.at(5)*nz3+f_out.at(2)*nz4)*uab9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        if (B4 > 1.E-10) {
            slip_energies.insert(B4);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << B4 << " B4 " << endl;
        }

        B6 = ((f_out.at(0)*nz2 + f_out.at(3)*nz3+f_out.at(4)*nz4)*uac7 + (f_out.at(3)*nz2 + f_out.at(1)*nz3+f_out.at(5)*nz4)*uac8+(f_out.at(4)*nz2 + f_out.at(5)*nz3+f_out.at(2)*nz4)*uac9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
// B6 = ((f_out.at(0)*(-3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(-4.08212913E-16) + (f_out.at(3)*(-3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(-3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(0))/alpha_t;
        if (B6 > 1.E-10) {
            slip_energies.insert(B6);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << B6 << " B6 " << endl;
        }

        C1 = ((f_out.at(0)*naa2 + f_out.at(3)*naa3+f_out.at(4)*naa4)*uy7 + (f_out.at(3)*naa2 + f_out.at(1)*naa3+f_out.at(5)*naa4)*uy8+(f_out.at(4)*naa2 + f_out.at(5)*naa3+f_out.at(2)*naa4)*uy9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
// C1 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(0) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(4.08212913E-16))/alpha_t;
        if (C1 > 1.E-10) {
            slip_energies.insert(C1);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << C1 << " C1 " << endl;
        }

        C3 = ((f_out.at(0)*naa2 + f_out.at(3)*naa3+f_out.at(4)*naa4)*uaa7 + (f_out.at(3)*naa2 + f_out.at(1)*naa3+f_out.at(5)*naa4)*uaa8+(f_out.at(4)*naa2 + f_out.at(5)*naa3+f_out.at(2)*naa4)*uaa9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
// C3 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(4.08212913E-16) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(0) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(4.08212913E-16))/alpha_t;
        if (C3 > 1.E-10) {
            slip_energies.insert(C3);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << C3 << " C3 " << endl;
        }

        C5 = ((f_out.at(0)*naa2 + f_out.at(3)*naa3+f_out.at(4)*naa4)*uac7 + (f_out.at(3)*naa2 + f_out.at(1)*naa3+f_out.at(5)*naa4)*uac8+(f_out.at(4)*naa2 + f_out.at(5)*naa3+f_out.at(2)*naa4)*uac9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
// C5 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(3.74999999E-12)+f_out.at(4)*(-3.74999999E-12))*(-4.08212913E-16) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(3.74999999E-12)+f_out.at(5)*(-3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(3.74999999E-12)+f_out.at(2)*(-3.74999999E-12))*(0))/alpha_t;
        if (C5 > 1.E-10) {
            slip_energies.insert(C5);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << C5 << " C5 " << endl;
        }

        D1 = ((f_out.at(0)*nab2 + f_out.at(3)*nab3+f_out.at(4)*nab4)*uy7 + (f_out.at(3)*nab2 + f_out.at(1)*nab3+f_out.at(5)*nab4)*uy8+(f_out.at(4)*nab2 + f_out.at(5)*nab3+f_out.at(2)*nab4)*uy9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        // D1 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(0) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(4.08212913E-16))/alpha_t;
        if (D1 > 1.E-10) {
            slip_energies.insert(D1);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << D1 << " D1 " << endl;
        }

        D4 = ((f_out.at(0)*nab2 + f_out.at(3)*nab3+f_out.at(4)*nab4)*uab7 + (f_out.at(3)*nab2 + f_out.at(1)*nab3+f_out.at(5)*nab4)*uab8+(f_out.at(4)*nab2 + f_out.at(5)*nab3+f_out.at(2)*nab4)*uab9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
// D4 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(-4.08212913E-16) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(0) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(4.08212913E-16))/alpha_t;
        if (D4 > 1.E-10) {
            slip_energies.insert(D4);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << D4 << " D4 " << endl;
        }
        D6 = ((f_out.at(0)*nab2 + f_out.at(3)*nab3+f_out.at(4)*nab4)*uad7 + (f_out.at(3)*nab2 + f_out.at(1)*nab3+f_out.at(5)*nab4)*uad8+(f_out.at(4)*nab2 + f_out.at(5)*nab3+f_out.at(2)*nab4)*uad9)*ae1*(pow(ae1,2)*ae2*ae4/ae3)/alpha_t;
        // D6 = ((f_out.at(0)*(3.74999999E-12) + f_out.at(3)*(-3.74999999E-12)+f_out.at(4)*(3.74999999E-12))*(4.08212913E-16) + (f_out.at(3)*(3.74999999E-12) + f_out.at(1)*(-3.74999999E-12)+f_out.at(5)*(3.74999999E-12))*(4.08212913E-16) +(f_out.at(4)*(3.74999999E-12) + f_out.at(5)*(-3.74999999E-12)+f_out.at(2)*(3.74999999E-12))*(0))/alpha_t;
        if (D6 > 1.E-10) {
            slip_energies.insert(D6);
            Out_logfile2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << D6 << " D6 " << endl;
        }

        energies_vector.push_back(slip_energies);
        /// E-level Numbers output
        Out_numbers2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << slip_energies.size() << " " << *slip_energies.rbegin() << endl;

        for(auto it = slip_energies.begin(); it != slip_energies.end(); ++it)
            if(next(it) != slip_energies.end())
                energy_gaps.insert(*next(it) - *it);

        egaps_vector.push_back(energy_gaps);

        Out_egaps2_stream << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << alpha_t << " " << energy_gaps.size() << " ";
        for(auto sev : energy_gaps) {
            if (sev > 1.E-10)
                Out_egaps2_stream << sev << " ";
        }
        Out_egaps2_stream << endl;

        slip_energies.clear();
        energy_gaps.clear();
    }
//    Out_logfile_stream << "S11(MPa)" << " " << "S22(MPa)" << " " << "S33(MPa)" << " " << "S12(MPa)" << " " << "S13(MPa)" << " " << "S23(MPa)" << " " << "Pressure(MPa)" << " " << "Mizes(MPa)" << " " << "T_factor" << " " << " alpha_t " << " " << " A2 " << " " << " A3 " << " "  << " A6 " << " " << " B2 " << " " << " B4 " << " "  << " B5 " << " " << " C1 " << " " << " C3 " << " "  << " C5 " << " " << " D1 " << " " << " D4 " << " " << " D6 " << endl;

//for(auto f_out : factors) {
//        Out_logfile_stream << f_out.at(0) << " " << f_out.at(1) << " " << f_out.at(2) << " " << f_out.at(3) << " " << f_out.at(4) << " " << f_out.at(5) << " " << f_out.at(6) << " " << f_out.at(7) << " " << f_out.at(8) << " " << " alpha_t " << " " << " e1 " << " " << " e2 " << " " << " e3 " << " " << " " << " e4 " << " " << " e5 " << " " << " e6 " << " " << " e7 " << " " << " e8 " << " " << " e9 " << " " << " e10 " << " " << " e11 " << " " << " e12 " << " " << endl;
    //   cout << f_out.at(0) << " " << f_out.at(1) << " " << f_out.at(2) << " " << f_out.at(3) << " " << f_out.at(4) << " " << f_out.at(5) << " " << f_out.at(0) << " " << f_out.at(1) << " " << f_out.at(2) << " " << " alpha_t " << " " << " e1 " << " " << " e2 " << " " << " e3 " << " " << " .. " << endl;
    // }

    std::vector<double> t_vector;
    for(auto f_out : factors) {
        t_vector.push_back(f_out.at(5));
    }
    std::sort(t_vector.begin(), t_vector.end());

    for(auto nu : t_vector) {
        cout << nu << endl;
    }
    Out_logfile2_stream.close(); // closing off-stream 'Processing_Design.log' file
    Out_numbers2_stream.close(); // closing off-stream 'Processing_Design.log' file

}
