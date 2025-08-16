///================================ PCC Design module ===================================================================================///
///=========================================================================================================================================///
///* Finds an optimal special and induced types of the State Vectors. *///
///* -----------------------------------------------------------------------------------------------------------------------------------*///
///* Created by Dr Elijah Borodin at the University of Manchester 2022-2023 years as a module of PCC Processing Design code (CPD code) *///
///* A part or the PRISB codes project (https://github.com/PRISBteam) supported by EPSRC UK via grant EP/V022687/1 in 2022-2023 years *///
/// (https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/V022687/1)                                                            *///
///==================================================================================================================================///
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

/// Attached user-defined C++ libraries:
// External
#include "../../../src/lib/external/Eigen/SparseCore"

// Internal
#include "../PCC_Support_Functions.h" // It must be here - first in this list (!)
#include "../PCC_Objects.h"
#include "../ini/ini_readers.h"

// Local
///---------------------------------------------------------
#include "functions/GeneticAlgorithm.h"
///---------------------------------------------------------

using namespace std; // standard namespace

/// External variables
extern std::vector<unsigned int> CellNumbs; //number of cells in a PCC defined globally
extern ofstream Out_logfile_stream;
extern std::string source_path;
extern std::string output_dir;
extern std::vector<std::string> paths_to_PCC_matrices; // PCCpaths to PCC files
extern int PCC_dimension; // PCC dimension: dim = 1 for graphs, dim = 2 for 2D plane polytopial complexes and dim = 3 for 3D bulk polyhedron complexes, as it is specified in the main.ini file.
extern std::vector<std::tuple<double, double, double>> node_coordinates_vector, edge_coordinates_vector, face_coordinates_vector, polytope_coordinates_vector; // coordinate vectors defined globally

#include "PCC_Design.h"
///* ========================================================= PCC SUBCOMPLEX FUNCTION ======================================================= *///
///* ========================================================================================================================================= *///
/*!
 * @details Output sequences of their evolution towards the optimum defined by a 'goal' function.
 * @param configuration
 * @return std::vector<std::vector<int>>
 */
std::vector<std::vector<int>> PCC_Design(Config &design_configuration, CellDesign &processing_cell_design){
/// Main output of the module
    std::vector<std::vector<int>> design_list_of_vectors;

    std::cout << "--- Genetic Algorithm ---" << std::endl;

    /// initial configuration from config/design.ini file
    config_reader_design(design_configuration);

    if (design_configuration.Get_design_mode() == "G"s && design_configuration.Get_design_max_generation_number() > 1) {
        int chromosome_length = CellNumbs.at(design_configuration.Get_design_cell_type());    // The length of the state vector
        int population_size = design_configuration.Get_design_population_size();
        double mutation_rate = design_configuration.Get_design_mutation_rate();
        double crossover_rate = design_configuration.Get_design_crossover_rate();
        const std::vector<int> possible_genes = {0, 1, 2, 3}; // Gene pool, e.g., {0, 1}
        int maximum_generation_number = design_configuration.Get_design_max_generation_number();

        // --- Fitness Function Selection ---
        // The desired fitness function is assigned to a std::function object.
        // This demonstrates the flexibility of the pluggable fitness function design.
        auto fitnessFunction = maximizingOnesFitness; /// (state_vector_by_sequence(processing_cell_design.Get_f_special_sequence(), design_configuration.Get_design_cell_type()));
        // auto fitnessFunction = shannonEntropyFitness;

        // --- GA Initialization ---
        GeneticAlgorithm ga(population_size, mutation_rate, crossover_rate, fitnessFunction);
        ga.initializePopulation(chromosome_length, possible_genes);

        std::cout << "\nStarting evolution..." << std::endl;
        std::cout << "Optimization Target: Maximizing the sum of genes in the state vector." << std::endl;

        // --- Evolution Loop ---
        for (int i = 0; i < maximum_generation_number; ++i) {
            ga.evolve();

            // Periodically print the progress of the optimization.
            if ((i + 1) % 10 == 0 || i == maximum_generation_number - 1) {
                Individual best = ga.getBestIndividual();
                std::cout << "Generation: " << ga.getGenerationCount()
                          << " | Best Fitness: " << best.fitness
                          << " | Preview: ";
                for (int k = 0; k < 10 && k < best.chromosome.size(); ++k) {
                    std::cout << best.chromosome[k];
                }
                std::cout << "..." << std::endl;
            }
        }

        // --- Final Results ---
        std::cout << "\nEvolution finished." << std::endl;
        Individual finalBest = ga.getBestIndividual();
        std::cout << "Final Best Fitness: " << finalBest.fitness << std::endl;
        std::cout << "Final Optimized State Vector (" << finalBest.chromosome.size() << " genes):" << std::endl;
        for (int gene: finalBest.chromosome) {
            std::cout << gene;
        }
        std::cout << std::endl;
    } // if (design_configuration.Get_design_mode() == 'G')

    return design_list_of_vectors;
}
