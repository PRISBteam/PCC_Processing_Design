/// * Function calculates the vector<int> "EdgeTypes" of types Edges in the PCC using its FES incidence matrix and special faces sequence (special_faces_sequence) * ///

std::vector<double> j_fractions_vector(vector<int> const &TJsTypes);

std::vector<double> d_fractions_vector(vector<int> const &TJsTypes);

double Configuration_Entropy(vector<int> const &TJsTypes);

double Configuration_Entropy(std::vector<double> const &j_fractions);

std::tuple<double, double> Configuration_Entropy_tuple(std::vector<double> const &j_fractions);

/// Configuration Entropy Differences
double Configuration_Entropy_change(vector<int> const &TJsTypes, vector<int> const &TJsTypesNew);

double Configuration_Entropy_change(std::vector<double> const &j_fractions, vector<int> const &TJsTypesNew);

double Configuration_DevEntropy_change(std::vector<double> const &j_fractions, vector<int> const &TJsTypesNew);