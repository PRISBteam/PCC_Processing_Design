#ifndef PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
#define PCC_PROCESSING_DESIGN_PCC_OBJECTS_H

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

#endif //PCC_PROCESSING_DESIGN_PCC_OBJECTS_H
