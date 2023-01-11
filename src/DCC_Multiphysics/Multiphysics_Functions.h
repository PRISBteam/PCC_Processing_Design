///================================ A part of the DCC Multiphysics module =============================================================///
///=================================================================================================================================///
/** The library contains functions saetting various physical quantities for k-Cells of PCC               *
 *                                                                                              **/
///================================================================================================================================///
using namespace std;
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