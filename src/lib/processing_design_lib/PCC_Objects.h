#ifndef PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
#define PCC_PROCESSING_DESIGN_PCC_OBJECTS_H

#include "../external/Eigen-5.0/SparseCore"
#include <set>

/// ==== # I # =============== Structure for the initial configuration  ========================= ///

/*!
 * @brief This class combine PCCpaths to directories and initial variables set in the config/_.ini files with the methods of their reading like
 * @public Get_config(), Set_config()
 * @param dim, source_dir, paths, ConfigVector
 * @return Configuration_aState, Configuration_gState
 */
class Config {

public:
    // Output only one time per simulation
    bool main_reader_switch = true, subcomplex_reader_switch = true, multiphysics_reader_switch = true, processing_reader_switch = true, kinetics_reader_switch = true, characterisation_reader_switch = true, design_reader_switch = true, writer_reader_switch = true;

private:
    struct main_configuration {
        std::string pcc_source_dir;
        std::string output_dir;
        std::string pcc_standard;
        std::string main_type;
    } main_config;

    struct subcomplex_configuration {
        std::string sctype;
        std::vector<double> plane_orientation;
        double cut_length;
        unsigned int grain_neighbour_orders;
        bool is_subcomplex_log_file;
    } subcomplex_config;

    struct multiphysics_configuration {
        std::string Mid_matrix;
        std::string Mid_inclusion;
        std::tuple<double, double, double> sample_dimensions;
        double multiphysics_time_scale;
        double stress_mode;
        double inclusion_stress_intensity;
        double crack_stress_intensity;
        int cut_direction_id;
        double min_crack_lenghts, max_crack_lenghts;
        int number_of_cracks;
        int number_of_crack_sizes;
        Eigen::MatrixXd external_stress_tensor;
        double equivalent_stress;
        double pressure;
        double ambient_temperature;
        std::vector<double> macrocrack_ini;
        bool is_multiphysics_log_file;
    } multiphysics_config;

    struct processing_configuration {
        std::string pp_mode, pf_mode, pe_mode, pn_mode;
        std::string ip_mode, if_mode, ie_mode, in_mode;
        int p_multiplexity, f_multiplexity, e_multiplexity, n_multiplexity;
        std::string pp_source_path, pf_source_path, pe_source_path, pn_source_path;
        //        std::vector<std::string> p_sequence_source_paths;
        std::vector<double> pn_max_fractions, pe_max_fractions, pf_max_fractions, pp_max_fractions;
        std::vector<double> in_max_fractions, ie_max_fractions, if_max_fractions, ip_max_fractions;

        double mu; double sigma; unsigned int bins_number;
        bool is_processing_log_file;
    } processing_config;

    struct kinetics_configuration {
        // general
        std::string nk_mode, ek_mode, fk_mode, pk_mode;
        std::string material_id;
        double kinetics_time_scale;

        // corrosion
        double kinetics_corrosion_rate_scale;
        double corrosion_activation_volume;

        //irradiation
        std::tuple<double,double,double> beam_direction;
        double irradiation_damage_rate;
        double beam_energy_flux;
        double beam_current;
        double energy_dissipation_rate;
        double observation_time;

        // module output
        bool is_kinetics_log_file;

    }kinetics_config;

    struct design_configuration {
        std::string goal_function_id;
        int design_cell_type;
        std::string design_mode;
        std::string design_goal;
        unsigned int population_size;
        double mutation_rate, crossover_rate, survival_rate;
        int max_generation_number;
        int design_genes_diversity;
        bool is_design_log;

    }design_config;

    int config_dim;
//    std::string config_source_dir, config_output_dir; // Input and output directories as it is written in the 'config/main.ini' file
//    std::string pcc_standard; // PCC standard as specified in the technical documentation for the project
    std::vector<std::string> config_PCCpaths; // PCC source directory
    std::vector<int> config_ConfVector; // main module keys
//    std::string config_main_type; // 'mode' from the config/main.ini file: 'LIST' (execution one by one all the active (ON) project modules), 'TUTORIAL' as a specific education mode, 'PERFORMANCE_TEST' or the 'TASK' mode :: This define the global simulation mode: 'LIST' for the "list" of modules implementing one by one (if ON) and 'TASK' for the user-defined task scripts with modules and functions included from the project's libraries
    std::string config_sim_task; // path to the corresponding *.cpp file containing a 'simulation task' (for 'TASK' execution mode only, not 'LIST') as it is written in the 'config/main.ini' file
    bool is_log_file; // key for various config outputs

    /// The list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_aState but for 'cracked' (or induced) network of k-cells
    /* where 'n' :: "nodes", 'e' :: "edges", 'f' :: "faces", and 'p' :: "polyhedrons" */
// State_Vector in the form : [Element index] - > [Element type], like [0, 0, 2, 1, 1, 0, 2, 4, 3, 3, 2, 0,... ...,2] containing all CellNumb.at(*) element types
    std::vector<unsigned int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // Normally the State_<*>_vector of special cells can be calculated based on the corresonding special_cell_sequences
    std::vector<unsigned int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // separate vectors containing the other 'fractured' labels different from the 'special' ones. To be calculated based on the corresonding fractured_cell_sequences

/// Configuration_aState = { State_p_vector, State_f_vector, State_e_vector, State_n_vector } is a list of all 'state vectors': from (1) State_p_vector (on top, id = 0) to (4) State_n_vector (bottom, id = 3)
    std::vector<std::vector<unsigned int>> Configuration_aState;
    std::vector<std::vector<unsigned int>> Configuration_gState;
    std::vector<std::vector<unsigned int>> Configuration_iState;

public:
    /// main
    void Set_main_type(std::string &main_type_str);
    void Set_pcc_source_dir(std::string &pcc_source_directory);
    void Set_output_dir(std::string &output_directory);
    void Set_pcc_standard_id(std::string &pcc_standard_id);

    /// subcomplex
    void Set_subcomplex_mode(std::string &subcomplex_mode);
    std::string Get_subcomplex_mode(void) const;
    void Set_subcomplex_plane(std::vector<double> &new_plane_orientation);
    std::vector<double> Get_subcomplex_plane(void) const;
    void Set_subcomplex_cut_length(double &new_cut_lenth);
    double Get_subcomplex_cut_length(void) const;
    void Set_subcomplex_grain_neighbour_orders(unsigned int &new_grain_neighbour_orders);
    unsigned int Get_subcomplex_grain_neighbour_orders(void) const;
    void Set_is_subcomplex_log_file(bool new_is_log_file);
    bool Get_is_subcomplex_log_file(void) const;

    /// multiphysics
    void Set_multiphysics_matrixMaterial_id(std::string &matrix_id);
    std::string Get_multiphysics_matrixMaterial_id(void) const;
    void Set_multiphysics_inclusionMaterial_id(std::string &inclusion_id);
    std::string Get_multiphysics_inclusionMaterial_id(void) const;
    void Set_multiphysics_time_scale(double &multiphysics_time_scale);
    double Get_multiphysics_time_scale(void) const;
    void Set_multiphysics_inclusion_stress_intensity_factor(double &inclusion_stress_intensity);
    double Get_multiphysics_inclusion_stress_intensity_factor(void) const;
    void Set_multiphysics_crack_stress_intensity_factor(double &crack_stress_intensity);
    double Get_multiphysics_crack_stress_intensity_factor(void) const;
    void Set_multiphysics_crack_stress_mode(double &stress_mode);
    double Get_multiphysics_crack_stress_mode(void) const;
    void Set_multiphysics_crack_grow_direction(int &cut_direction_id);
    int Get_multiphysics_crack_grow_direction(void) const;
    void Set_multiphysics_min_crack_lenghts(double &min_crack_lenghts);
    double Get_multiphysics_min_crack_lenghts(void) const;
    void Set_multiphysics_max_crack_lenghts(double &max_crack_lenghts);
    double Get_multiphysics_max_crack_lenghts(void) const;
    void Set_multiphysics_number_of_cracks(int &number_of_cracks);
    int Get_multiphysics_number_of_cracks(void) const;
    void Set_multiphysics_number_of_crack_sizes(int &number_of_crack_sizes);
    int Get_multiphysics_number_of_crack_sizes(void) const;
    bool Get_is_multiphysics_log_file(void) const;
    void Set_multiphysics_sample_dimensions(std::tuple<double, double, double> &sample_dimensions);
    std::tuple<double, double, double>  Get_multiphysics_sample_dimensions(void);
    void Set_multiphysics_external_stress_tensor(Eigen::MatrixXd &external_stress_tensor);
    Eigen::MatrixXd Get_multiphysics_external_stress_tensor(void);
    double Get_multiphysics_vonMises_stress_in_MPa(Eigen::MatrixXd external_stress_tensor);
    double Get_multiphysics_pressure_in_MPa(Eigen::MatrixXd external_stress_tensor);
    //    void Set_macrocrack_ini(std::vector<double> &macrocrack_ini);
    void Set_is_multiphysics_log_file(bool is_log_file);
    void Set_multiphysics_temperature(double &new_temperature);
    double Get_multiphysics_temperature(void) const;

    ///processing
    void Set_processing_pp_mode(std::string &new_pp_mode);
    void Set_processing_pf_mode(std::string &new_pf_mode);
    void Set_processing_pe_mode(std::string &new_pe_mode);
    void Set_processing_pn_mode(std::string &new_pn_mode);
    std::string Get_processing_pp_mode(void);
    std::string Get_processing_pf_mode(void);
    std::string Get_processing_pe_mode(void);
    std::string Get_processing_pn_mode(void);

    void Set_processing_ip_mode(std::string &new_ip_mode);
    void Set_processing_if_mode(std::string &new_if_mode);
    void Set_processing_ie_mode(std::string &new_ie_mode);
    void Set_processing_in_mode(std::string &new_in_mode);
    std::string Get_processing_ip_mode(void);
    std::string Get_processing_if_mode(void);
    std::string Get_processing_ie_mode(void);
    std::string Get_processing_in_mode(void);

    void Set_processing_p_multiplexity(int &pp_multiplexity_value);
    void Set_processing_f_multiplexity(int &pf_multiplexity_value);
    void Set_processing_e_multiplexity(int &pe_multiplexity_value);
    void Set_processing_n_multiplexity(int &pn_multiplexity_value);
    int Get_processing_p_multiplexity(void);
    int Get_processing_f_multiplexity(void);
    int Get_processing_e_multiplexity(void);
    int Get_processing_n_multiplexity(void);

    void Set_processing_p_source_path(std::string &new_p_source_path);
    void Set_processing_f_source_path(std::string &new_f_source_path);
    void Set_processing_e_source_path(std::string &new_e_source_path);
    void Set_processing_n_source_path(std::string &new_n_source_path);
    std::string Get_processing_p_source_path(void);
    std::string Get_processing_f_source_path(void);
    std::string Get_processing_e_source_path(void);
    std::string Get_processing_n_source_path(void);

    void Set_processing_pn_max_fractions(std::vector<double> &new_pn_max_fractions);
    void Set_processing_pe_max_fractions(std::vector<double> &new_pe_max_fractions);
    void Set_processing_pf_max_fractions(std::vector<double> &new_pf_max_fractions);
    void Set_processing_pp_max_fractions(std::vector<double> &new_pp_max_fractions);
    std::vector<double> Get_processing_pn_max_fractions(void);
    std::vector<double> Get_processing_pe_max_fractions(void);
    std::vector<double> Get_processing_pf_max_fractions(void);
    std::vector<double> Get_processing_pp_max_fractions(void);

    void Set_processing_in_max_fractions(std::vector<double> &new_in_max_fractions);
    void Set_processing_ie_max_fractions(std::vector<double> &new_ie_max_fractions);
    void Set_processing_if_max_fractions(std::vector<double> &new_if_max_fractions);
    void Set_processing_ip_max_fractions(std::vector<double> &new_ip_max_fractions);
    std::vector<double> Get_processing_in_max_fractions(void);
    std::vector<double> Get_processing_ie_max_fractions(void);
    std::vector<double> Get_processing_if_max_fractions(void);
    std::vector<double> Get_processing_ip_max_fractions(void);

    void Set_processing_mu(double &new_mu);
    void Set_processing_sigma(double &new_sigma);
    void Set_processing_bins_number(unsigned int &new_bins_numb);
    double Get_processing_mu(void);
    double Get_processing_sigma(void);
    unsigned int Get_processing_bins_number(void);

//    void Set_ptype_vector(std::vector<std::string> &ptype_vector);
//    void Set_ctype_vector(std::vector<std::string> &ctype_vector);
//    void Set_pindex_vector(std::vector<double> &pindex_vector);
    void Set_is_processing_log_file(bool is_processing_log_file);
    bool Get_is_processing_log_file(void);

    void Read_config(Config &main_configuration); // Read the 'initial configuration' of the problem set in all the relevant '*.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)
    void Set_config(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<std::vector<int>> Configuration_State, std::vector<std::vector<int>> Configuration_gState); // manual setting of the configuration

    main_configuration Get_main_config() const;
    subcomplex_configuration Get_subcomplex_config() const;
    multiphysics_configuration Get_multiphysics_config() const;
    processing_configuration Get_processing_config() const;

    int Get_dim() const; //@return dim
    std::vector<int> Get_ConfVector() const; //!@return ConfVector
    std::string Get_source_dir() const; //!@return source_dir
    std::string Get_output_dir() const; //!@return output_dir
    std::string Get_pcc_standard() const; //!@return pcc_standard
    std::vector<std::string> Get_paths() const; //!@return PCC PCCpaths
    std::string Get_main_type() const; //!@return main_type
    std::string Get_sim_task() const; //!@return sim_task path to the corresponding *.cpp file containing the task code

    std::vector<std::vector<unsigned int>> Get_Configuration_aState() const; //!@return Configuration_aState
    std::vector<std::vector<unsigned int>> Get_Configuration_gState() const; //!@return Configuration_gState
    std::vector<std::vector<unsigned int>> Get_Configuration_iState() const; //!@return Configuration_iState

    ///kinetics
    void Set_kinetics_nk_mode(std::string &new_nk_mode);
    void Set_kinetics_ek_mode(std::string &new_nk_mode);
    void Set_kinetics_fk_mode(std::string &new_nk_mode);
    void Set_kinetics_pk_mode(std::string &new_nk_mode);
    std::string Get_kinetics_nk_mode(void);
    std::string Get_kinetics_ek_mode(void);
    std::string Get_kinetics_fk_mode(void);
    std::string Get_kinetics_pk_mode(void);

    void Set_kinetics_material_id(std::string &new_mat_id);
    std::string Get_kinetics_material_id(void);

    void Set_is_kinetics_log_file(bool is_kinetics_log_file);
    bool Get_is_kinetics_log_file(void);

// kinetics corrosion
    void Set_kinetics_time_scale(double &new_time_parameter);
    double Get_kinetics_time_scale(void) const;
    void Set_kinetics_corrosion_rate_scale(double &new_corrosion_rate_parameter);
    double Get_kinetics_corrosion_rate_scale(void) const;
    void Set_kinetics_corrosion_activation_volume(double &corrosion_activation_volume);
    double Get_kinetics_corrosion_activation_volume(void) const;

// kinetics irradiation
    void Set_kinetics_beam_direction(std::tuple<double,double,double> &new_beam_direction);
    std::tuple<double,double,double> Get_kinetics_beam_direction(void) const;
    void Set_kinetics_irradiation_damage_rate(double &new_irradiation_damage_rate);
    double Get_kinetics_irradiation_damage_rate(void) const;
    void Set_kinetics_beam_energy_flux(double &beam_energy_flux);
    double Get_kinetics_beam_energy_flux(void) const;
    void Set_kinetics_beam_current(double &beam_current);
    double Get_kinetics_beam_current(void) const;
    void Set_kinetics_energy_dissipation_rate(double &energy_dissipation_rate);
    double Get_kinetics_energy_dissipation_rate(void) const;
    void Set_kinetics_observation_time(double &new_observation_time);
    double Get_kinetics_observation_time(void) const;

// design module
    void Set_design_goal_function_id(std::string &new_goal_function_id);
    std::string Get_design_goal_function_id(void) const;
    void Set_design_cell_type(int &new_design_cell_type);
    int Get_design_cell_type(void) const;
    void Set_design_mode(std::string &PCCDesign_type);
    void Set_design_goal(std::string &min_max_goal);
    std::string Get_design_goal(void) const;
    std::string Get_design_mode(void) const;
    void Set_design_genes_diversity(int &genes_diversity);
    int Get_design_genes_diversity(void) const;
    void Set_design_population_size(unsigned int &population_size);
    unsigned int Get_design_population_size(void) const;
    void Set_design_mutation_rate(double &mutation_rate);
    double Get_design_mutation_rate(void) const;
    void Set_design_crossover_rate(double &crossover_rate);
    double Get_design_crossover_rate(void) const;
    void Set_design_survival_rate(double &survival_rate);
    double Get_design_survival_rate(void) const;
    void Set_design_max_generation_number(int &max_generation_number);
    int Get_design_max_generation_number(void) const;
    void Set_is_design_log_file(bool is_design_log);
    bool Get_is_design_log_file(void) const;

}; // ConfigVector (../config/main.ini) contains ALL the control variables needed for the program execution

/// ==== # I # =============== END  ========================= ///

/// ==== # II # =============== Fundamental classes representing concepts of discrete combinatorial space ========================= ///

/// ==== # II.1 # =============== PCC (polytopal cell complex) class  ========================= ///
class PCC {
    std::set<bool> internal_grains_state_vector, internal_faces_state_vector, internal_edges_state_vector, internal_nodes_state_vector; // state vectors like [0 1 1 0 0 1 ...] where '1' signifies INTERNAL element (all its (k-1)-cells on the 1-boundary have adjacent neighbours) and '0' if not

protected:
    // list of polytopes
    std::vector<unsigned int> polytope_ids;

    /// Combinatorics
    // list of polytope k-boundaries and k-co-boundaries
    // For instance, for a FACE incident EDGES and NODES are its 1-boundary and 2-boundary, while incident GRAINS are on its co-boundary, etc.
    // A polytope, by definition, is on its own 0-boundary, 0-co-boundary and 0-neighbours
    // 'Neighbours' here MUST have at least one common 1-boundary cell
    std::vector<std::vector<unsigned int>> polytope_k_boundaries_list; // list of lists for each k-polytope in a PCC
    std::vector<std::vector<unsigned int>> polytope_k_coboundaries_list; // list of lists for each k-polytope in a PCC
    std::vector<std::vector<unsigned int>> polytope_k_neighbours_list; // list of lists for each k-polytope in a PCC

    ///Geometry
    // list of the lists of triplets of node coordinates: [0] - nodes, [1] - edges, [2] - faces, [3] - grains, [4] - 4-cells, etc
    std::vector<std::vector<std::tuple<double, double, double>>> cell_barycentre_coordinates;

    /// Measures
    // list of lists of geometric measures ('volumes' for 3-cells, 'areas' for 2-cells, 'lengths' for 1-cells, 'size' for 0-cells)
    std::vector<std::vector<std::tuple<double, double, double>>> cell_measures_vector;

public:
    void Set_edge_barycentre_coordinates(void);
    void Set_face_barycentre_coordinates(void);
    std::vector<std::tuple<double, double, double>> Get_edge_barycentre_coordinates(void);
    std::vector<std::tuple<double, double, double>> Get_face_barycentre_coordinates(void);


}; // end of class PCC

/// ==== # II.2 # =============== Skeleton class  ========================= ///
class Skeleton {
private:

public:

};
/// ========== END of class Skeleton functions description

/// ==== # II.3 # =============== Polytope class  ========================= ///

class Polytope {

    std::vector<std::tuple<double, double, double>> minmax_node_coordinates; // a vector containing two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}

private:
    // list of nodes
    std::vector<unsigned int> node_ids;
    // list of faces
    std::vector<unsigned int> faces_list;
    // list of neighbours (other polytopes)
    std::vector<unsigned int> neighbours_list;

    // list of triplets of node coordinates
    std::vector<std::tuple<double, double, double>> node_coordinates;

public:
    unsigned int grain_id;

    Polytope(unsigned int grain_new_id); // constructor 1

    void Set_node_ids(Eigen::SparseMatrix<double> const &GFS, Eigen::SparseMatrix<double> const &FES, Eigen::SparseMatrix<double> const &ENS);

    void Set_faces_list(Eigen::SparseMatrix<double> const &GFS);

    std::vector<unsigned int> Get_faces_list(void) const;

    /// return - vector of all node (vertices) coordinates of a polytope
    void Set_node_coordinates(std::vector<std::tuple<double,double,double>> &vertex_coordinates_vector);

    std::vector<unsigned int> Get_node_ids(void) const;

    std::vector<std::tuple<double, double, double>> Get_node_coordinates(void) const;

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a polytope with number grain_id
    std::vector<std::tuple<double, double, double>> Get_minmax_node_coordinates(void) const;

}; // end of class Polytope
/// ========== END of class Polytope functions description

/// ==== # II # =============== END ========================= ///

/// ==== # III # =============== Classes contained descriptions of the objects obtained as output of the project Modules  ========================= ///

/// ==== # III.1 # =============== Subcomplex class  ========================= ///
/*! @breif create a PCC subcomplex of the same dimension (maximum dimension k of its k-cells) as the parent PCC
 * @protected  sub_grains_set, sub_faces_set, sub_nodes_set, internal_faces_set // sets of k_max-cells and (k_max-1)-cells of a PCC
 * @protected sub_sfaces_set, internal_sfaces_set, sub_sfaces_sequence, internal_sub_sfaces_set // sets of special 'assigned' k_max-cells and (k_max-1)-cells of a PCC
 * sub_cfaces_sequence // sets of special 'generated' k_max-cells and (k_max-1)-cells of a PCC
 * @public subcomplex_id, sub_length, a_n, b_n, c_n, D_plane //
 * @public std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane} //
 * @function
 */
class Subcomplex {

protected:
    /// 1. Combinatorics
    std::set <unsigned int> sub_grains_set;
    std::set <unsigned int> sub_faces_set;
    std::set <unsigned int> sub_nodes_set;
    std::set <unsigned int> internal_faces_set;
    std::set <unsigned int> sub_sfaces_set;
    std::set <unsigned int> internal_sfaces_set;
    std::vector <unsigned int> sub_sfaces_sequence;
    std::set <unsigned int> internal_sub_sfaces_set;
    std::vector <unsigned int> sub_cfaces_sequence;

    /// 2. Geometry
    std::vector<std::tuple<double, double, double>> sub_face_coordinates;
    std::vector<std::tuple<double, double, double>> internal_sub_face_coordinates;
    std::vector <std::tuple<double, double, double>> sub_sfaces_coord;
    std::vector <std::tuple<double, double, double>> sub_cfaces_coord;

    std::vector<std::tuple<double, double, double>> sub_grain_coordinates;

public:
    unsigned int subcomplex_id;
    double sub_length;
    double a_n, b_n, c_n, D_plane;
    std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

    Subcomplex() {} // constructor 1
    Subcomplex(std::set <unsigned int> &new_sub_grains_set); // constructor 2

    std::vector <unsigned int> Get_sub_sfaces_sequence(void) const;
    std::vector <unsigned int> Get_sub_cfaces_sequence(void) const;
    std::vector <std::tuple<double, double, double>> Get_sub_sfaces_coord(void) const;
    std::vector <std::tuple<double, double, double>> Get_sub_cfaces_coord(void) const;

    void Set_sub_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence);
    void Set_sub_cfaces_sequence(std::vector <unsigned int> const &csub_faces_sequence);
    void Set_sub_sfaces_coord(std::vector<std::tuple<double, double, double>> const &sfaces_coord);
    void Set_sub_cfaces_coord(std::vector<std::tuple<double, double, double>> const &cfaces_coord);

    /// Polytope
    // sequence
    void Set_sub_polytope_set(std::set <unsigned int> &new_sub_grains_set);
    std::set <unsigned int> Get_sub_polytope_set(void) const;
    // geometry
    void Set_sub_polytope_coordinates(std::vector<std::tuple<double, double, double>> &new_sub_grain_coordinates);
    std::vector<std::tuple<double, double, double>> Get_sub_polytope_coordinates(void) const;

    /// Faces
    // sequence
    void Set_sub_faces_set(std::set <unsigned int> &new_sub_faces_set);
    std::set <unsigned int> Get_sub_faces_set(void) const;
    void Set_internal_sub_faces_set(std::set <unsigned int> &new_internal_faces_set);
    std::set <unsigned int> Get_internal_sub_faces_set(void) const;
    void Set_sub_sfaces_set(std::set <unsigned int> &new_sfaces_set);
    std::set <unsigned int> Get_sub_sfaces_set(void) const;

    void Set_internal_sub_sfaces_set(std::set <unsigned int> &new_internal_sfaces_set);
    std::set <unsigned int> Get_internal_sub_sfaces_set(void) const;

    // special and induced [c]('cracked') fqce sequences
//    void Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence);
//    std::vector <unsigned int> Get_sfaces_sequence(void) const;
//    void Set_cfaces_sequence(std::vector <unsigned int> &sub_cfaces_sequence);
//    std::vector <unsigned int> Get_cfaces_sequence(void) const;

    /// Geometry
    void Set_sub_face_coordinates(std::vector<std::tuple<double, double, double>> &new_sub_face_coordinates);
    std::vector<std::tuple<double, double, double>> Get_sub_face_coordinates(void) const;

    void Set_sub_internal_face_coordinates(std::vector<std::tuple<double, double, double>> &new_internal_face_coordinates);
    std::vector<std::tuple<double, double, double>> Get_sub_internal_face_coordinates(void) const;

    /// Edges
    // sequence
    void Set_sub_edges_set(std::set <unsigned int> &new_sub_faces_set);
    std::set <unsigned int> Get_sub_edges_set(void) const;
    // geometry

    /// Nodes
    // sequence
    void Set_sub_nodes_set(std::set <unsigned int> &new_sub_nodes_set);
    std::set <unsigned int> Get_sub_nodes_set(void) const;
    // geometry

}; // end of class Subcomplex

/// ==== # III.1.1 # =============== PCC Section service class -- used for creation Subcomplex plains be sectioning 3D cubes ========================= ///
/*! @breif create a geometric section of a 3D cube by a plane set by 4 real coefficients saved in this object
 * @private id // a section ID
 * @public a_coef, b_coeff, c_coeff, D_coeff   // coefficients of the plane equation in a 3D space:  a_coef*X + b_coeff*Y + c_coeff*Z + D_coeff = 0
 */
class PCCSection {
private:
    double id;
public:
    double a_coef;
    double b_coeff;
    double c_coeff;
    double D_coeff;
};

/// ==== # III.2 # =============== CellEnergies class  ========================= ///
///
/// ==== # III.2.1 # =============== Material service class -- used for creation Cell Energies with tabulated material characteristics ========================= ///
/*!
 * @brief List various tabulated physical and mechanical characteristics needed for defining physical energies.
 *          The list of material IDs can be taken from 'config/CPD_material_database'
 */
class Material {
private:
    std::string material_type = "material", inclusion_type = "inclusion";
    double mass_density = 0.0;
    double melting_point = 0.0;
    double gb_cohesion_energy = 0.0;
    double gb_width;
    double Burgers_vector;
    double Young_modulus = 0.0;
    double Poisson_ratio = 0.0;
    double yield_strength = 0.0;
    double strength = 0.0;
    double fracture_toughness = 0.0;

    double gb_inclusion1_adh_energy;
    double sface_energy_agglomeration;
    double inclusion_mass_density;

    double lagbs_corrosion_current;
    double hagbs_corrosion_current;
    double sigma3_corrosion_current;

public:
    Material(std::string Mid); // constructor 1
    Material(std::string Mid, std::string Iid); // constructor 2

    // Structural
    std::string Get_material_type(void) const;
    double Get_gb_width(void) const;
    double Get_Burgers_vector(void) const;

    // Thermodynamic
    double Get_mass_density(void) const;
    double Get_melting_point(void) const;
    double Get_gb_cohesion_energy(void) const;

    // Mechanical
    double Get_Young_modulus(void) const;
    double Get_Poisson_ratio(void) const;
    double Get_Yield_strength(void) const;
    double Get_Strength(void) const;
    double Get_Fracture_toughness(void) const;

    // Inclusions
    double Get_gb_inclusion1_adh_energy(void) const;
    std::string Get_inclusion_type(void) const;
    double Get_inclusion_agglomeration_energy(void) const;
    double Get_inclusion_mass_density(void) const;

    // Corrosion
    double Get_lagbs_corrosion_current(void) const;
    double Get_hagbs_corrosion_current(void) const;
    double Get_sigma3_corrosion_current(void) const;

};

/*! @breif create a list of the energy_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
 * @private homogeneous_elastic_energy  // double average value for an entire PCC
 * @private von_Mises_elastic_stress, ambient_temperature k_elastic_energies, k_thermal_energies, k_self_energies // vector<double> for each k-cell in a PCC
 */
class CellEnergies {
private:
    std::vector<double> von_Mises_elastic_stress;
    double homogeneous_elastic_energy = 0.0;
    std::vector<double> ambient_temperature;

    /// Energies for each cell in a PCC
    std::vector<double> p_elastic_energies, f_elastic_energies, e_elastic_energies, n_elastic_energies; // elastic energies of k-cells defined at their barycentres
    std::vector<double> p_thermal_energies, f_thermal_energies, e_thermal_energies, n_thermal_energies; // thermal energies of k-cells defined at their barycentres
    std::vector<double> p_self_energies, f_self_energies, e_self_energies, n_self_energies; // any associated self-energy including the cohesion energy of grain boundaries for 'f_self_energies'

public:
    /// Set of variables
    CellEnergies() {}; // constructor
    void Set_external_von_Mises_stress(std::tuple<double, double, double, double, double, double, double, double, double> &external_stress);
    void Set_von_Mises_stress(std::vector<double> &equivalent_stress);
    std::vector<double> Get_von_Mises_stress(void) const; // [Pa]
    void Set_ambient_temperature(std::vector<double> &new_ambient_temperature); // [K]
    std::vector<double> Get_ambient_temperature(void) const; // [K]
    void Set_homogeneous_elastic_energy(std::tuple<double, double, double> &sample_dimensions, double &von_Mises_elastic_stress, Material &matrix_material); // [J]
    void Set_p_elastic_energies(std::vector<double> p_el_energies); // in [J]
    void Set_f_elastic_energies(std::vector<double> f_el_energies); // in [J]
    void Set_e_elastic_energies(std::vector<double> e_el_energies); // in [J]
    void Set_n_elastic_energies(std::vector<double> n_el_energies); // in [J]

    void Set_p_self_energies(std::vector<double> p_el_energies); // in [J]
    void Set_f_self_energies(std::vector<double> f_el_energies); // in [J]
    void Set_e_self_energies(std::vector<double> e_el_energies); // in [J]
    void Set_n_self_energies(std::vector<double> n_el_energies); // in [J]

    // Get values
    double Get_homogeneous_elastic_energy(void); // [J]

    std::vector<double> Get_p_elastic_energies(void) const;
    std::vector<double> Get_f_elastic_energies(void) const;
    std::vector<double> Get_e_elastic_energies(void) const;
    std::vector<double> Get_n_elastic_energies(void) const;

    std::vector<double> Get_p_self_energies(void) const;
    std::vector<double> Get_f_self_energies(void) const;
    std::vector<double> Get_e_self_energies(void) const;
    std::vector<double> Get_n_self_energies(void) const;

}; // END of class CellEnergies

/// ==== # III.3 # =============== CellDesign class  ========================= ///
///
/// ==== # III.3.1 # =============== Agglomeration service class -- used for creating a collection of labels at each PCC's k-cell  ========================= ///

/*!
* @brief Objects of this class store collections of labels at each PCC's k-cell
* @private (string)atype, (unsigned int)aface_number, apower, a_average_strip_length
 */
class Agglomeration {
    double adhesion_energy = 0;
    double surface_energy = 0;
private:
    std::string atype; // like "rgo"
    unsigned int aface_number = 0;
    unsigned int apower = 0;
    unsigned int a_average_strip_length = 0;
public:
    Agglomeration(unsigned int AFace); // constructor 1
    Agglomeration(unsigned int AFace, unsigned int AglPower); // constructor 2 complex

    void Set_new_agglomeration(unsigned int AFace);
    void Set_agglomeration_type(std::string type);
    void Set_agglomeration_power(std::vector<std::vector<unsigned int>> const &RW_series_vector);
    void SetAvLength(std::vector<std::vector<unsigned int>> const &RW_series_vector); // Average length of strips related to this agglomeration

    unsigned int Get_agglomeration_kcell_number() const;
    int Get_agglomeration_power() const;
    int Get_agglomeration_power(std::vector<std::vector<unsigned int>> const &RW_series_vector); /// overloaded /// BAD
    int GetAvLength() const; /// BAD
    int GetAvLength(std::vector<std::vector<unsigned int>> const &RW_series_vector); /// overloaded /// BAD

}; // end of class agglomeration

/*!
 * @brief A CellDesign object contains 'state' or 'design' vectors contained particular configuration of labels on various PCC skeletons;
 * Moreover, it contains 'sequences' of special assigned, induced and generated cell numbers in 'historical' order of their appearance.
 */
class CellDesign {
private:
    bool is_set_p_special_sequence = false, is_set_f_special_sequence = false, is_set_e_special_sequence = false, is_set_n_special_sequence = false;
    bool is_set_p_induced_sequence = false, is_set_f_induced_sequence = false, is_set_e_induced_sequence = false, is_set_n_induced_sequence = false;
    bool is_set_p_special_design = false, is_set_f_special_design = false, is_set_e_special_design = false, is_set_n_special_design = false;
    bool is_set_p_induced_design = false, is_set_f_induced_design = false, is_set_e_induced_design = false, is_set_n_induced_design = false;

    /// Configurations/Designs: special cells and induced cells
    std::vector<unsigned int> p_special_design, f_special_design, e_special_design, n_special_design; // state vectors of special k-cells
    std::vector<unsigned int> p_induced_design, f_induced_design, e_induced_design, n_induced_design; // state vectors of induced k-cells

    /// Sequences: special cells and induced cells
    std::vector<unsigned int> p_special_sequence, f_special_sequence, e_special_sequence, n_special_sequence; // sequences of UNIQUE special k-cell numbers
    std::vector<unsigned int> p_induced_sequence, f_induced_sequence, e_induced_sequence, n_induced_sequence; // sequences of UNIQUE induced k-cell numbers

    std::vector<Agglomeration> p_agglomerations_map, f_agglomerations_map, e_agglomerations_map, n_agglomerations_map; // sequences of agglomerations

    /// Series: special cells and induced cells
    std::vector<std::vector<unsigned int>> p_special_series, f_special_series, e_special_series, n_special_series; // sequences of sequences of special k-cell numbers
    std::vector<std::vector<unsigned int>> p_induced_series, f_induced_series, e_induced_series, n_induced_series; // sequences of sequences of induced k-cell numbers

public:
    /// Set of variables
    CellDesign() {}; // constructor
    void Set_special_sequences(std::vector<unsigned int> psequence, std::vector<unsigned int> fsequence, std::vector<unsigned int> esequence, std::vector<unsigned int> nsequence);
    void Set_induced_sequences(std::vector<unsigned int> p_ind_sequence, std::vector<unsigned int> f_ind_sequence, std::vector<unsigned int> e_ind_sequence, std::vector<unsigned int> n_ind_sequence);
    void Set_designes(std::vector<unsigned int> pdesign, std::vector<unsigned int> fdesign, std::vector<unsigned int> edesign, std::vector<unsigned int> ndesign);
    void Set_induced_designs(std::vector<unsigned int> p_ind_design, std::vector<unsigned int> f_ind_design, std::vector<unsigned int> e_ind_design, std::vector<unsigned int> n_ind_design);
    void Set_special_sequence(std::vector<unsigned int> sequence, int ctype);
    void Set_special_series(std::vector<std::vector<unsigned int>> special_x_series, int cell_type);
    void Set_agglomeration_sequence(std::vector<Agglomeration> &agglomeration_x_sequence, int cell_type);
    void Set_induced_sequence(std::vector<unsigned int> ind_sequence, int ctype);
    void Set_induced_series(std::vector<std::vector<unsigned int>> induced_x_series, int cell_type);
    void Set_special_configuration(std::vector<unsigned int> design, int ctype);
    void Set_induced_design(std::vector<unsigned int> ind_design, int ctype);
/// void Set_agglomerations_special_sequence

    // Get
    std::vector<unsigned int> Get_p_special_sequence(void) const;
    std::vector<unsigned int> Get_f_special_sequence(void) const;
    std::vector<unsigned int> Get_e_special_sequence(void) const;
    std::vector<unsigned int> Get_n_special_sequence(void) const;

    std::vector<Agglomeration> Get_p_agglomeration_map(void) const;
    std::vector<Agglomeration> Get_f_agglomeration_map(void) const;
    std::vector<Agglomeration> Get_e_agglomeration_map(void) const;
    std::vector<Agglomeration> Get_n_agglomeration_map(void) const;
    std::vector<std::vector<unsigned int>> Get_p_special_series(void) const;
    std::vector<std::vector<unsigned int>> Get_f_special_series(void) const;
    std::vector<std::vector<unsigned int>> Get_e_special_series(void) const;
    std::vector<std::vector<unsigned int>> Get_n_special_series(void) const;

    void Set_p_design(std::vector<unsigned int> &p_special_vector);
    std::vector<unsigned int> Get_p_design(void) const;
    std::vector<unsigned int> Get_f_design(void) const;
    std::vector<unsigned int> Get_e_design(void) const;
    std::vector<unsigned int> Get_n_design(void) const;

    std::vector<unsigned int> Get_p_induced_sequence(void) const;
    std::vector<unsigned int> Get_f_induced_sequence(void) const;
    std::vector<unsigned int> Get_e_induced_sequence(void) const;
    std::vector<unsigned int> Get_n_induced_sequence(void) const;
    std::vector<std::vector<unsigned int>> Get_p_induced_series(void) const;
    std::vector<std::vector<unsigned int>> Get_f_induced_series(void) const;
    std::vector<std::vector<unsigned int>> Get_e_induced_series(void) const;
    std::vector<std::vector<unsigned int>> Get_n_induced_series(void) const;

    // Check
    bool Check_special_sequence(int cell_type);
    bool Check_induced_sequence(int cell_type);
    bool Check_special_design(int cell_type);
    bool Check_induced_design(int cell_type);

}; // END of class CellDesign

/// ==== # III.4 # =============== Processed Complex class  ========================= ///
///
/// ==== # III.4.1 # =============== Macrocrack service class -- used in Processed Complex stored a series of macrocracks ========================= ///

/*!
 * @brief
 */
class Macrocrack {
    double total_fracture_energy = 0;
    //Subcomplex half_plane_subcomplex; // geometry part

private:
    double a_n;
    double b_n;
    double c_n;
    double D_plane;
    double crack_length;
    double crack_stress_mode;
    double real_crack_length;
    std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

public:
    Subcomplex plane_subcomplex;
    int sub_id;

    int crack_id = 0;
    double surface_energy = 0;
    double bridging_energy = 0;
    double multiple_cracking_energy = 0;
    double stress_concentrators_energy = 0;

//    Macrocrack(int crack_id_new, Subcomplex &half_plane_sub); //constructor 1
    Macrocrack(int crack_id_new, int sub_id, Subcomplex &plane_sub, double crack_new_length, double crack_mode); //constructor 2

    double Get_crack_length(void) const;

    void Set_real_crack_length(double sample_size);

    void Set_sfaces_sequence(std::vector <unsigned int> const &special_faces_sequence);
    void Set_cfaces_sequence(std::vector <unsigned int> const &induced_faces_sequence);

//    std::vector <unsigned int> Get_sfaces_sequence(void) const;
    std::set <unsigned int> Get_sfaces_set(void) const;
    std::set <unsigned int> Get_internal_sfaces_set() const;
    std::vector <unsigned int> Get_cfaces_sequence(void) const;

    void Set_sfaces_coordinates(std::vector <std::tuple<double, double, double>> const &special_faces_coord);
    void Set_cfaces_coordinates(std::vector <std::tuple<double, double, double>> const &induced_faces_coord);

    std::vector <std::tuple<double, double, double>> Get_sfaces_coordinates(void) const;
    std::vector <std::tuple<double, double, double>> Get_internal_sfaces_coordinates(void) const;
    std::vector <std::tuple<double, double, double>> Get_cfaces_coordinates(void) const;

    double Get_real_crack_length() const;

    void Set_crack_plane();

    void Set_multiple_cracking_energy(CellEnergies &cell_energies_obj, std::vector<double> &cfaces_sequence);

    void Set_bridging_energy(double adhesion_energy, std::vector<double> &sfaces_sequence);

    double Get_multiple_cracking_energy() const;

    double Get_bridging_energy() const;

    std::set <unsigned int> Get_crack_faces_set() const;

//    std::vector <unsigned int> Get_sfaces_sequence() const;

    std::vector<double> Get_crack_plane() const;

    std::vector <std::tuple<double,double,double>> Get_common_faces_coordinates(unsigned int  crack_id) const;

}; // end of class MACROCRACK
/// ========== END of class Macrocrack functions description

/*!
 * @brief
 */
class ProcessedComplex { // Essential for Characterisation module
// PCC processed with all its characteristics and design sequences

private:
    std::vector<Macrocrack> macrocrack_growth_series;
    std::vector<std::set <unsigned int>> macrocrack_sfaces_series;
    std::vector<std::vector <unsigned int>> macrocrack_sfaces;
public:
    /// Set variables
    CellDesign pcc_design;

    void Set_design(CellDesign processed_pcc_design);
    void Set_macrocrack_sfaces(std::vector<std::vector <unsigned int>> &crack_growth_sf_series);
    void Set_macrocrack_sfaces_series(std::vector<std::set <unsigned int>> &crack_growth_sf_series);
    std::vector<std::set <unsigned int>> Get_macrocrack_sfaces_series(void) const;
    std::vector<std::vector <unsigned int>> Get_macrocrack_sfaces(void) const;

    // Sequences of special k-cells
    std::vector<std::vector<unsigned int>> face_process_seq;
    std::vector<std::vector<int>> face_process_state;

    // Entropic analysis
    std::vector<double> e_entropy_mean_vector, e_entropy_skrew_vector, e_entropy_full_vector;

    std::vector<std::vector<double>> je_fractions_sface_vector, de_fractions_sface_vector;
    std::vector<double> Betti_0_sface, Betti_1_sface, Betti_2_sface, inverse_connectivity_sface;

    std::vector<std::vector<double>> je_fractions_iface_vector, de_fractions_iface_vector;
    std::vector<double> Betti_0_iface, Betti_1_iface, Betti_2_iface, inverse_connectivity_iface;

    // Agglomerations
    std::vector<std::vector<Agglomeration>> agglomerations_in_powders; // Agglomeration(unsigned int AFace, unsigned int AglPower);


    // Analytical solutions
    std::vector<std::vector<double>> j_analytical_rand_vector, d_analytical_rand_vector;
    std::vector<std::vector<double>> j_analytical_cryst_vector, d_analytical_cryst_vector;
    std::vector<std::tuple<double, double>> AnRandEntropies_vector, AnCrystEntropies_vector;

    // Laplacian lab
    std::vector<std::vector<double>> Betti_vector;
}; // end of class ProcessedComplex

#endif //PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
