///================================ A part of the PCC Multiphysics module =============================================================///
///=================================================================================================================================///
/** The library contains functions setting internal energies, such as local elastic energies of stress concentrators                *
 * for all k-cells in a PCC.                                                                                                        */
///================================================================================================================================///

/// Standard C++ libraries (STL):
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// #include <execution> // Require C++ 17 and above

// external libraries
#include <Eigen/Core>
#include <Eigen/SparseCore>

// local libraries
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h"

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; // number of cells in a PCC defined globally
extern std::vector<std::string> PCCpaths; // PCCpaths to PCC files
extern int dim; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally
extern std::string output_dir;
//extern double D_coeff;

#include "multiphysics_internal_stresses.h"

///* ========================================================= PCC Internal Stress functions ======================================================= *///
///* ========================================================================================================================================= *///

/// ======# 1 #================= void Multiphysics_crack_stress_field() function ==============================================================///

using namespace std;

/// Just three overloaded templates for "sign()" function [taken from the internet]
template <typename TP> inline constexpr
int sign(TP u, std::false_type is_signed) {
    return TP(0) < u; }
template <typename TP> inline constexpr
int sign(TP u, std::true_type is_signed) {
    return (TP(0) < u) - (u < TP(0)); }
template <typename TP> inline constexpr
int sign(TP u) {
    return sign(u, std::is_signed<TP>()); }

std::vector<double> Multiphysics_crack_stress_field(Macrocrack &new_crack, Material &matrix_material, Eigen::MatrixXd &external_stress, std::tuple<double, double, double> &sample_dimensions) {

    std::ofstream Cracked_pcc_out;
    Cracked_pcc_out.open(output_dir + "Macrocrack_stress_field.txt"s, ios::trunc); // this Processing_Design.log stream will be closed at the end of the main function
    Cracked_pcc_out.close();

    Cracked_pcc_out.open(output_dir + "Macrocrack_stress_field.txt"s, ios::app); // this Processing_Design.log stream will be closed at the end of the main function


    std::vector<double> f_el_energy_densities; //function output
    double local_el_energy_density = 0.0, local_el_energy = 0.0;

    double Young_mod = matrix_material.Get_Young_modulus();
    double nu = matrix_material.Get_Poisson_ratio();
    double Shear_mod = Young_mod/ (2.0*(1.0+nu));
    double material_strenght = matrix_material.Get_Strength();

    double Len = new_crack.Get_crack_length(); //new crack length;
    Len = 0.4;  ///////////// ////////////////// !!!!!!
    Len *= 2.0; // doubled crack length
    if (new_crack.Get_crack_plane().at(0) > 0.0)
        Len *= std::get<0>(sample_dimensions);
    else
        Len *= std::get<1>(sample_dimensions);


    double Sigm = std::sqrt(0.5 * (pow((external_stress(0,0) - external_stress(1,1)), 2.0) + pow((external_stress(0,0) - external_stress(2,2)), 2.0) + pow((external_stress(1,1) - external_stress(2,2)), 2.0) ) + 3.0 * ( pow(external_stress(0,1), 2.0) + pow(external_stress(1,2), 2.0) + pow(external_stress(2,0), 2.0) )); //vonMises_stress

    double Sxx_crack = 0.0, Syy_crack = 0.0, Szz_crack = 0.0, Sxy_crack = 0.0;

   //for(auto fcv : face_coordinates_vector)
//REPAIR    cout << "fcv\t"s << face_coordinates_vector.size() << endl;

    double pcc_volume = std::get<0>(sample_dimensions)*std::get<1>(sample_dimensions)*std::get<2>(sample_dimensions);

    /// All face coordinates
    std::vector<tuple<double, double, double>> face_coordinates = Tuple3Reader(PCCpaths.at(13));

//REPAIR cout << "\t multi_dimensions\t" << get<0>(sample_dimensions) << " sd " << get<1>(sample_dimensions) << " sd " << get<2>(sample_dimensions) << endl;

/// Loop over all Faces
    unsigned int g1_n = 0, g2_n = 0; // numbers of two neighbouring grains incident to the same face in a PCC
    double g1_volume = 0.0, g2_volume = 0.0; // two grain volumes (relative: V = 1)

    for (unsigned int itr = 0; itr < face_coordinates.size(); ++itr) { // loop over all GBs and their barycentic coordinates
        double x = get<0>(face_coordinates.at(itr)) * std::get<0>(sample_dimensions);
        double y = get<1>(face_coordinates.at(itr)) * std::get<1>(sample_dimensions);
        double z = get<2>(face_coordinates.at(itr)) * std::get<2>(sample_dimensions);
//        cout << "\tLen\t" << Len << "\tx\t" << x << "\ty\t" << y << "\tz\t" << z << endl; //exit(0);

//        switch (crack_mode) {  case 1: { // crack mode I
        /// I CRACK MODE (Syy only!) ::
//        cout << " x " << x << " Len " << Len << endl;  exit(0);
        double D_val = 0.0; // pow((0.0*std::get<1>(sample_dimensions)),2.0)/4.0; // crack shift
        double r2p = pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0)/4.0 - D_val;
        double r2m = pow(x, 2.0) - pow(y, 2.0) - pow(Len, 2.0)/4.0 + D_val;
        double x2p = pow(x, 2.0) + pow(Len, 2.0)/4.0 - D_val;
        double x2m = pow(x, 2.0) - pow(Len, 2.0)/4.0 + D_val;
        double y2p = pow(y, 2.0) + pow(Len, 2.0)/4.0 - D_val;
        double y2m = pow(y, 2.0) - pow(Len, 2.0)/4.0 + D_val;

        double Pp = std::sqrt(pow(r2m, 2.0) + 4.0 * pow(x, 2.0) * pow(y,2)); //pow((y - Len/2.0), 2)); //        double Pp = std::sqrt(pow(r2m, 2.0) + 4.0 * pow(x, 2.0) * pow((y - (Len/2.0)),2)); //pow((y - Len/2.0), 2));


        double Qp = std::sqrt(Pp + r2m); // pow((x-std::sqrt(D_val)), 2.0)
        double Qm = std::sqrt(Pp - r2m);

        Sxx_crack = Sigm * (Pp * r2p * (abs(y) * Qm - abs(x) * Qp) -
                            abs(x) * Qp * (pow(y, 4.0) - 2.0 * x2p * pow(y, 2.0) - 3.0 * pow(x2m, 2.0)) +
                            abs(y) * Qm * (pow(y, 4.0) + 6.0 * x2p * pow(y, 2.0) + 5.0 * pow(x2m, 2.0))) /
                    (2.0 * std::sqrt(2.0) * pow(Pp, 3.0)); // + external_stress(0,0)

        Syy_crack = Sigm * (Pp * r2p * (abs(x) * Qp - abs(y) * Qm) +
                            abs(x) * Qp * (5.0 * pow(y, 4.0) + 6.0 * x2p * pow(y, 2.0) + pow(x2m, 2.0)) +
                            abs(y) * Qm * (3.0 * pow(y, 4.0) + 2.0 * x2p * pow(y, 2.0) - pow(x2m, 2.0))) /
                    (2.0 * std::sqrt(2.0) * pow(Pp, 3.0)); ///+ external_stress(1,1); // + external_stress(1,1)

        Szz_crack = nu * (Sxx_crack + Syy_crack);

        Sxy_crack = Sigm * sign(x * y) * (Pp * r2p * (abs(x) * Qm + abs(y) * Qp) -
                                          abs(x) * Qm * (3.0 * pow(y, 4.0) + 2.0 * x2p * pow(y, 2.0) - pow(x2m, 2.0)) +
                                          abs(y) * Qp * (pow(y, 4.0) - 2.0 * x2p * pow(y, 2.0) - 3.0 * pow(x2m, 2.0))) /
                    (2.0 * std::sqrt(2.0) * pow(Pp, 3.0)); // + external_stress(2,0)

///        local_el_energy_density = std::sqrt( 0.5*( pow((Sxx_crack - Syy_crack),2.0) + pow((Sxx_crack - Szz_crack),2.0) + pow((Syy_crack - Szz_crack),2.0) + 6.0*pow(Sxy_crack,2.0) ) );
        local_el_energy_density = pow((2.0*Young_mod),(-1.0))*(pow(Sxx_crack,2.0) + pow(Syy_crack,2.0) - 2.0*nu*Sxx_crack*Syy_crack );/// + pow((2.0*Shear_mod),(-1.0))*pow(Sxy_crack,2.0);
        f_el_energy_densities.push_back(local_el_energy_density);

/// (1.0 + nu) * (Sxx_crack + Syy_crack + Szz_crack) / 3.0
        double Vl = std::get<0>(sample_dimensions)*std::get<1>(sample_dimensions)*std::get<2>(sample_dimensions);
        Cracked_pcc_out << (pow(Sxx_crack,2) + pow(Syy_crack,2) - 2.0*nu*Sxx_crack*Syy_crack) / (2.0*Young_mod)  << "\t" << x*pow(10,6) << "\t" << y*pow(10,6) << "\t" << z*pow(10,6) << endl;

    } // for (unsigned int itr = 0; itr < face_coordinates.size(); ++itr)
    cout << "Fracture criterion:\t\t" <<  pow((2.0*Young_mod),(-1.0))*std::pow(material_strenght,2.0);

    /*
                double Pp = std::sqrt(pow(x, 2.0) - pow(y, 2.0) - pow(pow(Len, 2) / 4.0, 2.0) + 4.0 * pow(x, 2.0) * pow(y, 2.0));
                double Qp = std::sqrt(Pp + (pow(x, 2.0) - pow(y, 2.0) - pow(Len, 2.0) / 4.0));
                double Qm = std::sqrt(Pp - (pow(x, 2.0) - pow(y, 2.0) - pow(Len, 2.0) / 4.0));

                Sxx = Sigm * (
                        Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(y) * Qm - abs(x) * Qp) -
                        abs(x) * Qp * (pow(y, 4.0) - 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                       3.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                        abs(y) * Qm * (pow(y, 4.0) + 6.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) +
                                       5.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                Syy = Sigm * (
                        Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(x) * Qp - abs(y) * Qm) +
                        abs(x) * Qp * (5.0 * pow(y, 4.0) + 6.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) +
                                       pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                        abs(y) * Qm * (3.0 * pow(y, 4.0) + 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                       pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                //+ Sigm;
                Szz = nu * (Sxx + Syy);
                Sxy = Sigm * sign(x * y) *
                             (Pp * (pow(x, 2.0) + pow(y, 2.0) - pow(Len, 2.0) / 4.0) * (abs(x) * Qm + abs(y) * Qp) -
                              abs(x) * Qm *
                              (3.0 * pow(y, 4.0) + 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                               pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0)) +
                              abs(y) * Qp * (pow(y, 4.0) - 2.0 * (pow(x, 2.0) + pow(Len, 2.0) / 4.0) * pow(y, 2.0) -
                                             3.0 * pow((pow(x, 2.0) - pow(Len, 2.0) / 4.0), 2.0))
                             ) / (2.0 * std::sqrt(2.0) * pow(Pp, 3.0));
                cout << " Young_modulus: " << pow(x, 2.0) - pow(y, 2.0) - pow(pow(Len, 2) / 4.0, 2.0) + 4.0 * pow(x, 2.0) * pow(y, 2.0) << " face_energies.at(itr): " << face_energies.at(itr) << endl;
*/
    /// von Mizes energies
    /// (1) sample_dimensions.at(1) now! if (res_vonMizes_stress >  0.0) face_energies.at(itr) = pow(res_vonMizes_stress,2) * face_areas_vector.at(itr) * pow(sample_dimensions.at(1), 2.0) * GB_width / Young_modulus;

/**
            }
            case 2: { // crack mode II
                /// The second and the third (Sxy only) crack modes ::
                // S2xx := 0.5*Sigm* Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) + Abs[x]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 + 5.0*(x^2 - Len^2/4.0)^2) - Abs[y]*Qp*(3.0*y^4 + 10.0*(x^2 + Len^2/4.0)*y^2 + 7.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
                //S2yy := -0.5*Sigm* Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) - Abs[x]*Qm*(3.0*y^4 + 2.0*(x^2 + Len^2/4.0)*y^2 - 1.0*(x^2 - Len^2/4.0)^2) + Abs[y]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 - 3.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
                //S2xy := -0.5* Sigm* (Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[y]*Qm - Abs[x]*Qp) - Abs[x]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 - 3.0*(x^2 - Len^2/4.0)^2) + Abs[y]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 + 5.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
            }
            case 3: { // crack mode III
                /// The third crack mode ::
                //dl = 0.0
                //xp := x + dl
                //xm = x - dl
                //Pp3 := ((xp^2 - y^2 + y*Len)^2 + 4.0*xp^2*(y - Len/2.0)^2)^0.5
                //Qm3 := (Pp3^2 - (xp^2 - y^2 + y*Len))^0.5
                //Qp3 := (Pp3^2 + (xp^2 - y^2 + y*Len))^0.5
                //S3xy := 0.5* Sigm* (Pp3^2*(xp^2 + y^2 - y*Len)*(Qp3*(xp - y + Len/2.0) + Qm3*(xp + y - Len/2.0)) + Qm3*(xp + y - Len/2.0)*(1.0*xp^4 - 2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - 3.0*y^2*(y - Len)^2) + Qp3*(y - Len/2.0)*(3.0*xp^4 + 2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - y^2*(y - Len)^2) + Qp3*xp*(xp^4 + 6.0*xm^2*(y^2 - y*Len + Len^2/2.0) + 5.0*y^2*(y - Len)^2))/(2.0*(2.0^(0.5))*Pp3^6)
            }
        } // end of switch

        double SigmMAX = (1.0 + nu) * (Sxx + Syy + Szz) / 3.0;

        /// von Mizes stress
        double res_vonMizes_stress = 0.0;
        if (x <= Len)  res_vonMizes_stress = SigmMAX; /// assumption for the maximal stress value (!)
        else res_vonMizes_stress = std::sqrt(0.5 * (pow((Sxx - Syy),2.0) + pow((Sxx - Szz),2.0) + pow((Syy - Szz),2.0) + 6 * pow(Sxy,2.0)));

        /// von Mizes energies
        /// (1) sample_dimensions.at(1) now!
        if (res_vonMizes_stress >  0.0) face_energies.at(itr) = pow(res_vonMizes_stress,2) * face_areas_vector.at(itr) * pow(sample_dimensions.at(1), 2.0) * GB_width / Young_modulus;
//REPAIR        cout << " itr: " << itr << " face_energies.at(itr): " << face_energies.at(itr) << endl;

    } // end of for( itr )
**/
    Cracked_pcc_out << endl << endl;
    Cracked_pcc_out.close();
    cout << "\t p_crit = \t" << pow(200.0,2.0)*pow(10.0,12.0)/(2.0*Young_mod) << endl;
///    exit(0);

 return f_el_energy_densities;
} /// end of crack_modes_stress_field() function


/// HEAP
/*
std::vector<double> Energy_GBs_state() {
    std::vector<double> total_GBenergies_vector; // main function output
    std::vector<grain_boundary> all_GBs; // vector of all GBs as structures

    for (unsigned int m = 0; m < CellNumbs.at(2); m++)
        all_GBs.push_back(grain_boundary(m));

    face_areas_vector

    GB_SE_vector.at(m) = 2.0*_surface_energy*face_areas_vector.at(m);
    GB_EEE_vector.at(m) =
    //von Mizes stress
    double s_mis = sqrt(0.5*((s11 - s22)^2 + (s22 - s33)^2 + (s33 - s11)^2 + 6.0*(s12^2 + s23^2 + s31^2)));
    double ext_elastic_energy = pow(s_mis,2)*face_areas_vector.at(m)*hm/Young_modulus;

    GB_CIE_vector =
    GB_BLE_vector =
    GB_CLE_vector =

    for(unsigned int m = 0; m < CellNumbs.at(2); m++) {
        all_GBs.at(m).Set_surface_energy(GB_SE_vector);
        all_GBs.at(m).Set_external_elastic_energy(GB_EEE_vector);
        all_GBs.at(m).Set_crack_interaction_energy(GB_CIE_vector);
        all_GBs.at(m).Set_Bl_energy(GB_BLE_vector);
        all_GBs.at(m).Set_Cl_energy(GB_CLE_vector);
    } // end of for (unsigned int m = 0; m < CellNumbs.at(2); m++)

    for(unsigned int m = 0; m < CellNumbs.at(2); m++)
        total_GBenergies_vector.push_back(all_GBs.at(m).Get_surface_energy(m) + all_GBs.at(m).Get_external_elastic_energy(m) +
                                                  all_GBs.at(m).Get_crack_interaction_energy(m) + all_GBs.at(m).Get_Bl_energy(m) + all_GBs.at(m).Get_Bl_energy(m));
        return total_GBenergies_vector;
} // end of vector<double> Energy_GBs_state()

/*
void data_multiphysics_reader() {
    std::string line;
    ifstream in_data(config);

    if (in_data) { //If the file was successfully open, then
        while(getline(in_data, line, '\n'))
            for (auto it: line) {

                if (it == '?') {
            stringstream line1_stream(line);
            line1_stream >> p_max;
    } else cout << "SIMULATION_MODE() error: The file " << config << " cannot be read" << endl; // If something goes wrong


        // Number of types
        for (auto it: line) {
            if (it == '&') Processing_type = line.at(0); // simulation type
            else if (it == '`') Kinetic_type = line.at(0); // simulation Kinetic type

            else if (it == '@') res.push_back(line.at(0) - '0'); // dimension of the problem; res[1]
            else if (it == '!') res.push_back(line.at(0) - '0'); // number of i(X) designs; res[0]
                ///?? does it properly working with 10, 100 etc
                //else if (it == '!') res.push_back(line.at(0) - '0'); // number of special Face types; res[0]

            else if (it == '?') {
                stringstream line1_stream(line);
                line1_stream >> p_max;
                res.push_back(p_max);
            } // MAX fraction of Faces (calculation limit); res[2]


        }
*/