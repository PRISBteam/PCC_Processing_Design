#ifndef PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
#define PCC_PROCESSING_DESIGN_PCC_OBJECTS_H

#include <Eigen/SparseCore>
#include <set>

/// ==== # 4 # =============== Agglomeration class  ========================= ///

/*!
* @brief This class combine PCCpaths to directories and initial variables set in the config/_.ini files with the methods of their reading like
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

/// ==== # 0 # =============== Structure for the initial configuration  ========================= ///

/*!
 * @brief This class combine PCCpaths to directories and initial variables set in the config/_.ini files with the methods of their reading like
 * @public Get_config(), Set_config()
 * @param dim, source_dir, paths, ConfigVector
 * @return Configuration_sState, Configuration_cState
 */
class Config {
    private:
        int config_dim;
        std::string config_source_dir, config_output_dir; // Input and output directories as it is written in the 'config/main.ini' file
        std::string pcc_standard; // PCC standard as specified in the technical documentation for the project
        //std::vector<char*> config_PCCpaths;
        std::vector<std::string> config_PCCpaths;
        std::vector<int> config_ConfVector; // main module keys
        std::string config_main_type; // 'mode' from the config/main.ini file: 'LIST' (execution one by one all the active (ON) project modules), 'TUTORIAL' as a specific education mode, 'PERFORMANCE_TEST' or the 'TASK' mode :: This define the global simulation mode: 'LIST' for the "list" of modules implementing one by one (if ON) and 'TASK' for the user-defined task scripts with modules and functions included from the project's libraries
        std::string config_sim_task; // path to the corresponding *.cpp file containing a 'simulation task' (for 'TASK' execution mode only, not 'LIST') as it is written in the 'config/main.ini' file

    /// The list of all mentioned below State_<*>_vectors and State_<*>fracture_vectors as the output of the Processing module // is the list of 'state vectors' analogous to the Configuration_sState but for 'cracked' (or induced) network of k-cells
    /* where 'n' :: "nodes", 'e' :: "edges", 'f' :: "faces", and 'p' :: "polyhedrons" */
// State_Vector in the form : [Element index] - > [Element type], like [0, 0, 2, 1, 1, 0, 2, 4, 3, 3, 2, 0,... ...,2] containing all CellNumb.at(*) element types
        std::vector<unsigned int> State_p_vector, State_f_vector, State_e_vector, State_n_vector; // Normally the State_<*>_vector of special cells can be calculated based on the corresonding special_cell_sequences
        std::vector<unsigned int> State_pfracture_vector, State_ffracture_vector, State_efracture_vector, State_nfracture_vector; // separate vectors containing the other 'fractured' labels different from the 'special' ones. To be calculated based on the corresonding fractured_cell_sequences

/// Configuration_sState = { State_p_vector, State_f_vector, State_e_vector, State_n_vector } is a list of all 'state vectors': from (1) State_p_vector (on top, id = 0) to (4) State_n_vector (bottom, id = 3)
        std::vector<std::vector<unsigned int>> Configuration_sState;
        std::vector<std::vector<unsigned int>> Configuration_cState;

public:
    void Read_config(); // Read the 'initial configuration' of the problem set in all the relevant '*.ini' files containing in the '\config' project directory using the functions from the 'ini_readers.cpp' project library (and only from there)
    void Set_config(const std::vector<int> &ConfigVector, const std::string &source_dir, int &dim, std::vector<char*> paths, std::vector<std::vector<int>> Configuration_State, std::vector<std::vector<int>> Configuration_cState); // manual setting of the configuration

    int Get_dim() const; //@return dim
    std::vector<int> Get_ConfVector() const; //!@return ConfVector
    std::string Get_source_dir() const; //!@return source_dir
    std::string Get_output_dir() const; //!@return output_dir
    std::string Get_pcc_standard() const; //!@return pcc_standard
    std::vector<std::string> Get_paths() const; //!@return PCC PCCpaths
    std::string Get_main_type() const; //!@return main_type
    std::string Get_sim_task() const; //!@return sim_task path to the corresponding *.cpp file containing the task code

    std::vector<std::vector<unsigned int>> Get_Configuration_sState() const; //!@return Configuration_sState
    std::vector<std::vector<unsigned int>> Get_Configuration_iState() const; //!@return Configuration_iState
};
// ConfigVector (../config/main.ini) contains ALL the control variables needed for the program execution

/// ==== # 1 # =============== CellDesign class  ========================= ///
class CellDesign {
private:
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

}; // END of class CellDesign

/// # 7 # The class of a MATERIAL
/*!
 *
 */
class Material {
private:
    std::string material_type = "material", inclusion_type = "inclusion";
    double mass_density = 0.0;
    double melting_point = 0.0;
    double gb_cohesion_energy = 0.0;
    double gb_width;
    double Young_modulus = 0.0;
    double Poisson_ratio = 0.0;
    double yield_strength = 0.0;
    double strength = 0.0;
    double fracture_toughness = 0.0;

    double gb_inclusion1_adh_energy;
    double sface_energy_agglomeration;
    double inclusion_mass_density;

public:
    Material(std::string Mid); // constructor 1
    Material(std::string Mid, std::string Iid); // constructor 2

// Structural
    std::string Get_material_type(void) const;
    double Get_gb_width(void) const;

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

};

/// ========== END of class Materials functions description

/// ==== # 1.2 # =============== CellEnergies class  ========================= ///
// # V # The class of a CELLS_ENERGIES :: list of the energy_vectors corresponding to different dimensions 'k' of the k-cells in a PCC
class CellEnergies {
private:
    double von_Mises_elastic_stress = 0.0;
    double homogeneous_elastic_energy = 0.0;

    /// Energies for each cell in a PCC
    std::vector<double> p_elastic_energies, f_elastic_energies, e_elastic_energies, n_elastic_energies; // elastic energies of k-cells defined at their barycentres
    std::vector<double> p_thermal_energies, f_thermal_energies, e_thermal_energies, n_thermal_energies; // thermal energies of k-cells defined at their barycentres
    std::vector<double> p_self_energies, f_self_energies, e_self_energies, n_self_energies; // any associated self-energy including the cohesion energy of grain boundaries for 'f_self_energies'

public:
    /// Set of variables
    CellEnergies() {}; // constructor
    void Set_von_Mises_stress(std::tuple<double, double, double, double, double, double, double, double, double> &external_stress); // [Pa]
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
    double Get_von_Mises_stress(void); // [Pa]
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

/// ==== # 3 # =============== Subcomplex class  ========================= ///

class Subcomplex {

protected:
    /// 1. Combinatorics
    std::set <unsigned int> sub_grains_set;
    std::set <unsigned int> sub_faces_set;
    std::set <unsigned int> sub_nodes_set;
    std::set <unsigned int> internal_faces_set;
    std::vector <unsigned int> sub_sfaces_sequence;
    std::vector <unsigned int> sub_cfaces_sequence;

    /// 2. Geometry
    std::vector<std::tuple<double, double, double>> internal_face_coordinates;
    std::vector <std::tuple<double, double, double>> sub_sfaces_coord;
    std::vector <std::tuple<double, double, double>> sub_cfaces_coord;

    std::vector<std::tuple<double, double, double>> sub_grain_coordinates;

public:
    unsigned int subcomplex_id;
    double sub_length;
    double a_n; double b_n; double c_n; double D_plane;
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
    void Set_internal_faces_set(std::set <unsigned int> &new_internal_faces_set);
    std::set <unsigned int> Get_internal_faces_set(void) const;
    
    // special and induced [c]('cracked') fqce sequences
///    void Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence);
///    std::vector <unsigned int> Get_sfaces_sequence(void) const;
///    void Set_cfaces_sequence(std::vector <unsigned int> &sub_cfaces_sequence);
///    std::vector <unsigned int> Get_cfaces_sequence(void) const;

    // geometry
    void Set_internal_face_coordinates(std::vector<std::tuple<double, double, double>> &new_internal_face_coordinates);
    std::vector<std::tuple<double, double, double>> Get_internal_face_coordinates(void) const;

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


/// ==== # x # =============== PCC Section class  ========================= ///

class PCCSection {

private:
    double id;

public:
    double a_coef;
    double b_coeff;
    double c_coeff;
    double D_coeff;

};


/// ==== # 5 # =============== PCC class  ========================= ///
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

/// ==== # X # =============== Polytope class  ========================= ///


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

/// # 6 # The class of a MACROCRACK
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

    std::vector <unsigned int> Get_sfaces_sequence(void) const;
    std::vector <unsigned int> Get_cfaces_sequence(void) const;

    void Set_sfaces_coordinates(std::vector <std::tuple<double, double, double>> const &special_faces_coord);
    void Set_cfaces_coordinates(std::vector <std::tuple<double, double, double>> const &induced_faces_coord);

    std::vector <std::tuple<double, double, double>> Get_sfaces_coordinates(void) const;
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

/// ==== # 2 # =============== Processed Complex class  ========================= ///

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
