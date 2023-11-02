#ifndef PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
#define PCC_PROCESSING_DESIGN_PCC_OBJECTS_H

#include <Eigen/SparseCore>

/// ==== # 1 # =============== CellsDesign class  ========================= ///
class CellsDesign{
private:
    std::vector<unsigned int> p_sequence, f_sequence, e_sequence, n_sequence;
    std::vector<unsigned int> p_induced_sequence, f_induced_sequence, e_induced_sequence, n_induced_sequence;

    std::vector<int> p_design, f_design, e_design, n_design;
    std::vector<int> p_induced_design, f_induced_design, e_induced_design, n_induced_design;

public:
    /// Set of variables
    CellsDesign() {}; // constructor
    void Set_sequences(std::vector<unsigned int> psequence, std::vector<unsigned int> fsequence, std::vector<unsigned int> esequence, std::vector<unsigned int> nsequence);
    void Set_induced_sequences(std::vector<unsigned int> p_ind_sequence, std::vector<unsigned int> f_ind_sequence, std::vector<unsigned int> e_ind_sequence, std::vector<unsigned int> n_ind_sequence);
    void Set_designes(std::vector<int> pdesign, std::vector<int> fdesign, std::vector<int> edesign, std::vector<int> ndesign);
    void Set_induced_designs(std::vector<int> p_ind_design, std::vector<int> f_ind_design, std::vector<int> e_ind_design, std::vector<int> n_ind_design);
    void Set_sequence(std::vector<unsigned int> sequence, int ctype);
    void Set_induced_sequence(std::vector<unsigned int> ind_sequence, int ctype);
    void Set_design(std::vector<int> design, int ctype);
    void Set_induced_design(std::vector<int> ind_design, int ctype);

    // Get
    std::vector<unsigned int> Get_p_sequence(void);
    std::vector<unsigned int> Get_f_sequence(void);
    std::vector<unsigned int> Get_e_sequence(void);
    std::vector<unsigned int> Get_n_sequence(void);
    std::vector<unsigned int> Get_p_induced_sequence(void);
    std::vector<unsigned int> Get_f_induced_sequence(void);
    std::vector<unsigned int> Get_e_induced_sequence(void);
    std::vector<unsigned int> Get_n_induced_sequence(void);

    std::vector<int> Get_p_design(void);
    std::vector<int> Get_f_design(void);
    std::vector<int> Get_e_design(void);
    std::vector<int> Get_n_design(void);
};

/// ==== # 2 # =============== Processed Complex class  ========================= ///

class ProcessedComplex { // Essential for Characterisation module
// PCC processed with all its characteristics and design sequences

private:

public:
    /// Set variables
    CellsDesign pcc_design;

    void Set_design(CellsDesign processed_pcc_design);

    // Sequences of special k-cells
    std::vector<std::vector<unsigned int>> face_process_seq;
    std::vector<std::vector<int>> face_process_state;

    // Entropic analysis
    std::vector<double> e_entropy_mean_vector, e_entropy_skrew_vector, e_entropy_full_vector;
    std::vector<std::vector<double>> je_fractions_vector, de_fractions_vector;

    // Analytical solutions
    std::vector<std::vector<double>> j_analytical_rand_vector, d_analytical_rand_vector;
    std::vector<std::vector<double>> j_analytical_cryst_vector, d_analytical_cryst_vector;
    std::vector<std::tuple<double, double>> AnRandEntropies_vector, AnCrystEntropies_vector;

    // Laplacian lab
    std::vector<std::vector<double>> Betti_vector;
}; /// end of class ProcessedComplex

/// ==== # 3 # =============== Subcomplex class  ========================= ///

class Subcomplex {

private:
    /// 1. Combinatorics
    std::vector <unsigned int> sub_grains_sequence;
    std::vector <unsigned int> sub_faces_sequence;
    std::vector <unsigned int> sub_nodes_sequence;
    std::vector <unsigned int> common_faces_sequence;
    std::vector <unsigned int> s_sub_faces_sequence;
    std::vector <unsigned int> c_sub_faces_sequence;

    /// 2. Geometry
    //vector<tuple<double, double, double>> vertex_coordinates;
    std::vector<std::tuple<double, double, double>> common_faces_coordinates;
    std::vector<std::tuple<double, double, double>> sub_grain_coordinates;

public:
    unsigned int subcomplex_id;
    double sub_length;
    double a_n; double b_n; double c_n; double D_plane;
    std::vector<double> crack_plane = {a_n, b_n, c_n, D_plane};

    Subcomplex(unsigned int subcomplex_id_new);
    Subcomplex(unsigned int subcomplex_id_new, std::vector <unsigned int> new_sub_grains_sequence);

    /// Grains
    void Set_grains_sequence(std::vector <unsigned int> new_sub_grains_sequence);
    std::vector <unsigned int> Get_grains_sequence(unsigned int subcomplex_id);
    /// geometry
    void Set_sub_grain_coordinates(std::vector<std::tuple<double, double, double>> new_sub_grain_coordinates);
    std::vector<std::tuple<double, double, double>> Get_sub_grain_coordinates(unsigned int subcomplex_id);
    /// Faces
    void Set_faces_sequence(std::vector <unsigned int> new_sub_faces_sequence);
    std::vector <unsigned int>  Get_faces_sequence(unsigned int  subcomplex_id);
    void Set_common_faces_sequence(std::vector <unsigned int> new_common_faces_sequence);

    std::vector <unsigned int> Get_common_faces_sequence(unsigned int subcomplex_id);

    void Set_sfaces_sequence(std::vector <unsigned int> const &ssub_faces_sequence);
    std::vector <unsigned int> Get_sfaces_sequence(unsigned int  subcomplex_id);
    void Set_cfaces_sequence(std::vector <unsigned int> csub_faces_sequence);
    std::vector <unsigned int> Get_cfaces_sequence(unsigned int  subcomplex_id);

    /// geometry
    void Set_common_faces_coordinates(std::vector<std::tuple<double, double, double>> new_common_faces_coordinates);
    std::vector<std::tuple<double, double, double>> Get_common_faces_coordinates(unsigned int subcomplex_id);
};

/// ==== # 4 # =============== grain3D class  ========================= ///

class Polytope {

    std::vector<std::tuple<double, double, double>> minmax_node_coordinates; // a vecor containing two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}

private:
    // list of nodes
    std::vector<unsigned int> node_ids;

    // list of Faces
    std::vector<unsigned int> Faces_list;

    // list of triplets of nodes coordinated
    std::vector<std::tuple<double, double, double>> node_coordinates;

public:
    unsigned int grain_id;

    Polytope(unsigned int grain_new_id); // constructor 1

    void Set_node_ids(unsigned int grain_id, Eigen::SparseMatrix<double> const &GFS, Eigen::SparseMatrix<double> const &FES, Eigen::SparseMatrix<double> const &ENS);

    void Set_Faces_list(unsigned int grain_id, Eigen::SparseMatrix<double> const &GFS);

    std::vector<unsigned int> Get_Faces_list();

    /// return - vector of all node (vertices) coordinates of a grain
    void Set_node_coordinates(unsigned int grain_id);

    std::vector<unsigned int> Get_node_ids(unsigned int grain_id);

    std::vector<std::tuple<double, double, double>> Get_node_coordinates(unsigned int grain_id);

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    std::vector<std::tuple<double, double, double>> Get_minmax_node_coordinates(unsigned int grain_id);

}; // end of class Polytope

/// ========== END of class Polytope functions description

#endif //PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
