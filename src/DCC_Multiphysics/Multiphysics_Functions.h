///================================ A part of the DCC Multiphysics module =============================================================///
///=================================================================================================================================///
/** The library contains functions saetting various physical quantities for k-Cells of PCC               *
 *                                                                                              **/
///================================================================================================================================///
using namespace std;

void crack_1mode_stress_field(double Puasson_coeff = 0.3, double new_crack_length, double external_vonMizes_stress, vector<tuple<double, double, double>> const &vertex_coordinates){
    double nu = Puasson_coeff;
    double Len = new_crack_length;
    double Sigm = external_vonMizes_stress;

    for (unsigned int itr = 0; itr < vertex_coordinates.size(); ++itr) {
        double x = get<0>(vertex_coordinates.at(itr));
        double y = get<1>(vertex_coordinates.at(itr));
        double z = get<2>(vertex_coordinates.at(itr));

        /// The first crack mode ::
        double Pp = ((pow(x, 2.0) - pow(y, 2.0) - pow(pow(Len, 2) / 4.0), 2.0) + 4.0 * (pow(x, 2.0)) * (pow(y, 2.0))) ^ (0.5);
        double Qp = (Pp + (x ^ 2 - y ^ 2 - Len ^ 2 / 4.0)) ^ 0.5;
        double Qm = (Pp - (x ^ 2 - y ^ 2 - Len ^ 2 / 4.0)) ^ 0.5;
        double Sxx = Sigm * (Pp * (x ^ 2 + y ^ 2 - Len ^ 2 / 4.0) * (Abs[y] * Qm - Abs[x] * Qp) -
                             abs(x) * Qp * (y ^ 4 - 2.0 * (x ^ 2 + Len ^ 2 / 4.0) * y ^ 2 -
                                                                                        3.0 * (x ^ 2 - Len ^ 2 / 4.0) ^
                                            2) +
                             Abs[y] * Qm * (y ^ 4 + 6.0 * (x ^ 2 + Len ^ 2 / 4.0) * y ^ 2 +
                                                                                        5.0 * (x ^ 2 - Len ^ 2 / 4.0) ^
                                            2)) / (2.0 * (2.0 ^ (0.5)) * Pp ^ 3);
        double Syy = Sigm * (Pp * (x ^ 2 + y ^ 2 - Len ^ 2 / 4.0) * (Abs[x] * Qp - Abs[y] * Qm) +
                             Abs[x] * Qp * (5.0 * y ^ 4 +
                                                      6.0 * (x ^ 2 + Len ^ 2 / 4.0) * y ^ 2 + (x ^ 2 - Len ^ 2 / 4.0) ^
                                            2) +
                             Abs[y] * Qm * (3.0 * y ^ 4 +
                                                      2.0 * (x ^ 2 + Len ^ 2 / 4.0) *
                                                      y ^ 2 - (x ^ 2 - Len ^ 2 / 4.0) ^ 2)) /
                     (2.0 * (2.0 ^ (0.5)) * Pp ^ 3) + Sigm;
        double Szz = nu * (Sxx + Syy);
        double Sxy = Sigm *
                     Sign[x * y] * (Pp * (x ^ 2 + y ^ 2 - Len ^ 2 / 4.0) * (Abs[x] * Qm + Abs[y] * Qp) -
                                    Abs[x] * Qm * (3.0 * y ^ 4 +
                                                             2.0 * (x ^ 2 + Len ^ 2 / 4.0) * y ^
                                                   2 - (x ^ 2 - Len ^ 2 / 4.0) ^ 2) +
                                    Abs[y] * Qp * (y ^ 4 - 2.0 * (x ^ 2 + Len ^ 2 / 4.0) * y ^ 2 -
                                                                                               3.0 *
                                                                                               (x ^ 2 - Len ^ 2 / 4.0) ^
                                                   2)) / (2.0 * (2.0 ^ (0.5)) * Pp ^ 3);
        /// The second crack mode ::
        /*S2xx := 0.5*Sigm*
                Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) +
                              Abs[x]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 +
                                                                           5.0*(x^2 - Len^2/4.0)^2) -
                              Abs[y]*Qp*(3.0*y^4 + 10.0*(x^2 + Len^2/4.0)*y^2 +
                                                                            7.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
        S2yy := -0.5*Sigm*
                Sign[x*y]*(Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[x]*Qm + Abs[y]*Qp) -
                              Abs[x]*Qm*(3.0*y^4 + 2.0*(x^2 + Len^2/4.0)*y^2 -
                                                                           1.0*(x^2 - Len^2/4.0)^2) +
                              Abs[y]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 -
                                                                           3.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
        S2xy := -0.5*
                Sigm* (Pp^2*(x^2 + y^2 - Len^2/4.0)*(Abs[y]*Qm - Abs[x]*Qp) -
                          Abs[x]*Qp*(1.0*y^4 - 2.0*(x^2 + Len^2/4.0)*y^2 -
                                                                       3.0*(x^2 - Len^2/4.0)^2) +
                          Abs[y]*Qm*(1.0*y^4 + 6.0*(x^2 + Len^2/4.0)*y^2 +
                                                                       5.0*(x^2 - Len^2/4.0)^2))/(2.0*(2.0^(0.5))*Pp^6)
        dl = 0.0
        xp := x + dl
        xm = x - dl
       /// The third crack mode ::
        Pp3 := ((xp^2 - y^2 + y*Len)^2 + 4.0*xp^2*(y - Len/2.0)^2)^0.5
        Qm3 := (Pp3^2 - (xp^2 - y^2 + y*Len))^0.5
        Qp3 := (Pp3^2 + (xp^2 - y^2 + y*Len))^0.5
        S3xy := 0.5*
                Sigm* (Pp3^2*(xp^2 + y^2 - y*Len)*(Qp3*(xp - y + Len/2.0) +
                                                   Qm3*(xp + y - Len/2.0)) +
                           Qm3*(xp + y - Len/2.0)*(1.0*xp^4 -
                                                          2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - 3.0*y^2*(y - Len)^2) +
                           Qp3*(y - Len/2.0)*(3.0*xp^4 +
                                                     2.0*xm^2*(y^2 - y*Len + Len^2/2.0) - y^2*(y - Len)^2) +
                           Qp3*xp*(xp^4 + 6.0*xm^2*(y^2 - y*Len + Len^2/2.0) +
                                                 5.0*y^2*(y - Len)^2))/(2.0*(2.0^(0.5))*Pp3^6)
    */
        double SigmMAX = (1.0 + nu) * (Sxx + Syy + Szz) / 3.0;
        double Mizes = (0.5 * ((Sxx - Syy) ^ 2 + (Sxx - Szz) ^ 2 + (Syy - Szz) ^ 2 + 6 * Sxy ^ 2)) ^ 0.5;
    } // end of for( itr )
} /// end of crack_1mode_stress_field() function

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