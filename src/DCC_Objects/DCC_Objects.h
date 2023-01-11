/// Classes of special objects related to defect structures on a PCC elements

#ifndef agglomeration_H
#define agglomeration_H

/// (1) The class of agglomeration of defects in one special face element of a PCC
class agglomeration {
    double adhesion_energy = 0;
    double surface_energy = 0;

private:
    string atype; // like "rgo"
    unsigned int aface_number = 0;
    int          apower = 0;
    int          a_average_strip_length = 0;

public:
    agglomeration(unsigned int AFace) { // constructor 1 simple
        aface_number = AFace;
    }

    agglomeration(unsigned int AFace, int AglPower) { // constructor 2 complex
        aface_number = AFace;
        apower = AglPower;
    }

    void SetAFace(unsigned int AFace)
    {
        aface_number = AFace;
    }

    void SetAType(std::string type)
    {
        atype = type;
        if (atype == "rgo") {
            surface_energy = 0.2; // units [J/m^2]
            adhesion_energy = 0.4; // units [J/m^2]
        }

    }

    void SetAPower(vector<vector<int>> const &RW_series_vector) /// Agglomerations power
    {
        int AglPower = 0;
        //vector<int> aggl_vector(CellNumbs.at(2), 0); // state vector for agglomerations (with the size # cells initially filled with all 0s) : contains # of faces with agglomerations
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) AglPower += 1;

        apower = AglPower;
    }

    void SetAvLength(vector<vector<int>> const &RW_series_vector) /// Average length of strips related to this agglomeration
    {
        int ATotalLength = 0;
        for (auto RWsfv : RW_series_vector)
            for (auto RWsf: RWsfv)
                if(RWsf == aface_number) ATotalLength += RWsfv.size();

        a_average_strip_length = ATotalLength/ (double) GetAPower(RW_series_vector);
    }

    unsigned int GetAFace()
    {
        return aface_number;
    }

    int GetAPower()
    {
        return apower;
    }

    int GetAPower(vector<vector<int>> const &RW_series_vector)
    {
        if (apower != 0) {
            return apower;
        } else {
            int AglPower = 0;
            for (auto RWsfv : RW_series_vector)
                for (auto RWsf: RWsfv)
                    if(RWsf == aface_number) AglPower += 1;

            apower = AglPower;
            return apower;
        }
    }

            int GetAvLength()
    {
        return a_average_strip_length;
    }

    int GetAvLength(vector<vector<int>> const &RW_series_vector) {
        if (a_average_strip_length != 0) {
            return a_average_strip_length;
        } else {
            int ATotalLength = 0;
            for (auto RWsfv: RW_series_vector)
                for (auto RWsf: RWsfv)
                    if (RWsf == aface_number) ATotalLength += RWsfv.size();

                a_average_strip_length = ATotalLength / (double) GetAPower(RW_series_vector);
            return a_average_strip_length;
        }
    }

}; // end of class agglomeration

#endif


/// (2) The class of stress concentrators related to defects in one special face element of a PCC

class face_concentrator {
    double total_elastic_energy = 0;

private:
    //string conc_type; // bl or cl
    unsigned int conc_face_number = 0;
    unsigned int bl_index = 0;
    unsigned int cl_index = 0;
    double Bl_elastic_energy = 0;
    double Cl_elastic_energy = 0;
public:

    face_concentrator (unsigned int ConcFace) { // constructor 1 simple
        conc_face_number = ConcFace;
    }

    face_concentrator (unsigned int ConcFace, unsigned int Bl, unsigned int Cl) { // constructor 2 complex
        conc_face_number = ConcFace;
        bl_index = Bl;
        cl_index = Cl;
    }

    unsigned int GetBl_index()
    {
        return bl_index;
    }

    int GetCl_index()
    {
        return cl_index;
    }
}; // end of class

/// (3) The class of Grains in a PCC

class grain3D {
    vector<tuple<double, double, double>> minmax_node_coordinates; // a vecor containing two tuples: gmincoord{xmin,ymin,zmin},gmaxcoord{xmax,ymax,zmax}

private:
    // list of nodes
    vector<unsigned int> node_ids;

    // list of triplets of nodes coordinated
    vector<tuple<double, double, double>> node_coordinates;

public:

    unsigned int grain_id;

    grain3D(unsigned int grain_new_id) { // constructor 1 simple
        grain_id = grain_new_id;
    }

    void Set_node_ids(unsigned int grain_id, SpMat const &GFS, SpMat const &FES, SpMat const &ENS) { // set the node ids
        /// GFS -> FES -> ENS
        if(node_ids.size() == 0) {
            for(unsigned int l = 0; l < CellNumbs.at(2); l++) {// over all Faces (l)
                if (GFS.coeff(grain_id, l) == 1) {
                    for (unsigned int j = 0; j < CellNumbs.at(1); j++) // over all Edges (j)
                        if (FES.coeff(l, j) == 1) { // at the chosen Face with ID = 'l'
                            for (unsigned int i = 0; i < CellNumbs.at(0); i++) // over all Nodes
                                if (ENS.coeff(j, i) == 1) node_ids.push_back(i); // at the chosen Face with ID = 'l'
                        } // end of if (FES.coeff(l, j) == 1)
                } // end of (GFS.coeff(m, l) == 1)
            } // end of for(unsigned int l = 0; l < CellNumbs.at(2); l++) - Faces

        }/// end of if(node_ids.size() == 0)

    }

    /// return - vector of all node (vertices) coordinates of a grain
    void Set_node_coordinates(unsigned int grain_id, vector<tuple<double, double, double>> const &vertex_coordinates) { // set the node ids from Tr = triplet list
        vector<unsigned int> node_ids;
        for (auto  itr = node_ids.begin(); itr != node_ids.end(); ++itr)
            if(find(node_ids.begin(), node_ids.end(), distance(node_ids.begin(), itr)) != node_ids.end())
                node_coordinates.push_back(vertex_coordinates.at(*itr)); // vector<unsigned int> node_ids;
    }

    vector<unsigned int> Get_node_ids(unsigned int grain_id) { // set the node ids
        return node_ids;
    }

    vector<tuple<double, double, double>> Get_node_coordinates(unsigned int grain_id) { // set the node ids
        if (node_coordinates.size() != 0) {
            return node_coordinates;
        } else if (node_coordinates.size() > 0) {
            Set_node_coordinates(grain_id, vertex_coordinates);
            return node_coordinates;
        } else {
            throw std::invalid_argument(
                    "Please call Set_node_coordinates(unsigned int grain_id, vector<tuple<double, double, double>> const &vertex_coordinates) method first!");
            return node_coordinates;
        }
    }// end of  Get_node_coordinates() method

    /// return - vector with two tuples : { x_min, y_min, z_min; x_max, y_max, z_max} of a grain witn number grain_id
    vector<tuple<double, double, double>> Get_minmax_node_coordinates(unsigned int grain_id) { // min and max {x,y,z} values of vertices for a grain
        // class Grain3D function Get_node_coordinates(grain_id)
        vector<tuple<double, double, double>> minmax_tuple;
        ///Get_node_coordinates(grain_id)
        vector<tuple<double, double, double>> tup_node_coordinates = Get_node_coordinates(grain_id);
        // separating in three parts
        vector<double> x_node_coordinates, y_node_coordinates, z_node_coordinates;

        for (auto itr = tup_node_coordinates.begin(); itr != tup_node_coordinates.end(); ++itr) {
            x_node_coordinates.push_back(get<0>(*itr));
            y_node_coordinates.push_back(get<1>(*itr));
            z_node_coordinates.push_back(get<2>(*itr));
        }
//        std::for_each(tup_node_coordinates.begin(), tup_node_coordinates.end(), );

        auto itx = std::minmax_element(x_node_coordinates.begin(), x_node_coordinates.end());
        int xmin = *itx.first;
        int xmax = *itx.second;

        auto ity = std::minmax_element(y_node_coordinates.begin(), y_node_coordinates.end());
        int ymin = *ity.first;
        int ymax = *ity.second;

        auto itz = std::minmax_element(z_node_coordinates.begin(), z_node_coordinates.end());
        int zmin = *itz.first;
        int zmax = *itz.second;

        minmax_tuple = {make_tuple(xmin, ymin, zmin), make_tuple(xmax, ymax, zmax)};

        return minmax_tuple;
    }

}; // end of class grain3D

/// (4) The class of a Microcrack

class macrocrack {
    int crack_id = 0;
    double total_fracture_energy = 0;

private:
    double a_n;
    double b_n;
    double c_n;
    double D_plane;
    double crack_length;

    double stress_concentrators_energy = 0;
    double bridging_elastic_energy = 0;


public:

    macrocrack(int crack_id_new) { // constructor 1, simple
        crack_id = crack_id_new;
        a_n = 0.0;
        b_n = 0.0;
        c_n = 1.0;
        D_plane = 0.0;
        crack_length = 0.0;
    }

    macrocrack(int crack_id_new, double a_coeff, double b_coeff, double c_coeff, double D_coeff, double crack_length_new) { // constructor 2, more specific
        crack_id = crack_id_new;
        a_n = a_coeff;
        b_n = b_coeff;
        c_n = c_coeff;
        D_plane = a_coeff;
        crack_length = crack_length_new;
    }
/*
    double Get_surface_energy()
    {
        if (crack_length != 0.0)
            return ;
        else
            return 0.0;
    }
*/

}; // end of class MACROCRACK

class grain_boundary{

public:
    unsigned int GB_id; // grain boundary ID

    grain_boundary(unsigned int GB_new_id) { // class simple constructor
        GB_id = GB_new_id;
    }

/// GB state
    bool is_inclusion; // 0 - no, 1 - yes
    bool is_fractured; // 0 - no, 1 - yes
    bool is_agglomeration; // 0 - no, 1 - yes

/// Set values methods

    void Set_surface_energy(vector<double> GB_SE_vector ){
        surface_energy = GB_SE_vector.at(GB_id);
    }

    void Set_external_elastic_energy(vector<double> GB_EEE_vector) { /// EEE: external elastic energy vector
        external_elastic_energy = GB_EEE_vector.at(GB_id);
    }

    void Set_crack_interaction_energy(vector<double> GB_CIE_vector) {
        crack_interaction_energy = GB_CIE_vector.at(GB_id);
    }

    void Set_Bl_energy(vector<double> GB_BLE_vector) {
        Bl_energy = GB_BLE_vector.at(GB_id);
    }

    void Set_Cl_energy(vector<double> GB_CLE_vector){
        Cl_energy = GB_CLE_vector.at(GB_id);
    }

/// Get values methods
    double Get_surface_energy(unsigned int GB_id){
        if(surface_energy != 0)
            return surface_energy;
        else {
            Set_surface_energy(GB_SE_vector);
            return surface_energy; }
    }

    double Get_external_elastic_energy(unsigned int GB_id) { /// EEE: external elastic energy vector
        if(external_elastic_energy != 0)
            return external_elastic_energy;
        else {
            Set_external_elastic_energy(GB_EEE_vector);
            return external_elastic_energy; }
    }

    double Get_crack_interaction_energy(unsigned int GB_id) {
        if(crack_interaction_energy != 0)
            return crack_interaction_energy;
        else {
            Set_crack_interaction_energy(GB_CIE_vector);
            return crack_interaction_energy; }
    }

    double Get_Bl_energy(int GB_id){
        if(Bl_energy != 0)
            return Bl_energy;
        else {
            Set_Bl_energy(GB_BLE_vector);
            return Bl_energy; }
    }

    double Get_Cl_energy(unsigned int GB_id){
        if(Cl_energy != 0)
            return Cl_energy;
        else {
            Set_Cl_energy(GB_CLE_vector);
            return Cl_energy; }
    }

    double Get_total_energy(){

       /// total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
    }

private:

/// Combinatorial
    int GB_edges_number;  //equal to the number of neighbours
    int GB_nodes_number;

///Geometry
    double GB_area;
    double GB_perimeter;

    tuple<double,double,double> GB_barycentre_coordinates;

///Energies
    double surface_energy;
    double external_elastic_energy;
    double crack_interaction_energy;
    double Bl_energy;
    double Cl_energy;
    double total_energy; //= surface_energy + external_elastic_energy + Bl_energy + Cl_energy;
};

