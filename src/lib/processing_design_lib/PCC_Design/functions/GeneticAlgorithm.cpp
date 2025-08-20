#include <iostream>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <numeric> // For std::accumulate
#include <map>
#include <set>
#include <cmath>   // For std::log2
#include <vector>

// local libraries

#include "GeneticAlgorithm.h"

#include "../../lib/processing_design_lib/PCC_Kinetics/PCC_Kinetics.h"
#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h"

/*!
 * @brief A simple fitness function that calculates the sum of genes.
 * @details The goal for the GA when using this function is to maximize the
 *          sum of values in the chromosome, effectively evolving it towards
 *          a vector of all '1's if the gene pool is {0, 1}.
 * @param chromosome The state vector (individual) to evaluate.
 * @return The fitness score, which is the sum of the genes.
 */
double maximizingOnesFitness(const std::vector<int>& chromosome) {
    return std::accumulate(chromosome.begin(), chromosome.end(), 0);
}

/*!
 * @brief Calculates the Shannon entropy for a given chromosome.
 * @details This function measures the diversity or unpredictability of the genes.
 *          A higher entropy value indicates a more diverse mix of genes. The goal
 *          for the GA would be to find a state vector that maximizes this diversity.
 * @param chromosome The state vector (individual) to evaluate.
 * @return The calculated Shannon entropy as the fitness score.
 */
double shannonEntropyFitness(const std::vector<int>& chromosome) {
    if (chromosome.empty()) {
        return 0.0;
    }

    std::map<int, int> counts;
    for (int gene : chromosome) {
        counts[gene]++;
    }

    double entropy = 0.0;
    double chromosome_size = static_cast<double>(chromosome.size());
    for (auto const& [gene, count] : counts) {
        double probability = static_cast<double>(count) / chromosome_size;
        if (probability > 0) {
            entropy -= probability * std::log2(probability);
        }
    }
    return entropy;
} // end of shannonEntropyFitness()

/*!
 *
 * @param chromosome
 * @param design_configuration
 * @param processing_cell_design
 * @return
 */
double IrradiationDamageFitness(const std::vector<int>& chromosome, Config &design_configuration, CellDesign &processing_cell_design){
    std::vector<std::vector<double>> p_cells_history;

    std::vector<unsigned int> new_p_special_vector;
    for (auto ps : chromosome)
        new_p_special_vector.push_back(ps);

    processing_cell_design.Set_p_design(new_p_special_vector);
    p_cells_history = PCC_Kinetics(design_configuration, processing_cell_design);

    std::cout << "IrradiationDamageFitness :: fraction\t\t" << p_cells_history.back().at(2) << std::endl;
    return p_cells_history.back().at(2); // area fraction of fractured GBs
}



/// --- Helper for random number generation ---
// A simple utility to get a random double between 0.0 and 1.0
double randomDouble() {
    static std::mt19937 generator(std::random_device{}());
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// A simple utility to get a random integer in a range
int randomInt(int min, int max) {
    static std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

// ==============================================================================
// CONSTRUCTOR
// ==============================================================================
GeneticAlgorithm::GeneticAlgorithm(int popSize, double mutRate, double crossRate, std::function<double(const std::vector<int>&, Config&, CellDesign&)> fitnessFunc)
    : populationSize(popSize),
      mutationRate(mutRate),
      crossoverRate(crossRate),
      fitnessFunction(std::move(fitnessFunc)),
      generationCount(0) {
    if (!fitnessFunction) {
        throw std::invalid_argument("Fitness function cannot be null.");
    }
}

// ==============================================================================
// PUBLIC METHODS
// ==============================================================================

void GeneticAlgorithm::initializePopulation(int chromosomeLength, const std::vector<int>& possibleGenes) {
    if (possibleGenes.empty()) {
        throw std::invalid_argument("Possible genes cannot be empty.");
    }
    this->genePool = possibleGenes;
    population.clear();
    population.resize(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        population[i].chromosome.resize(chromosomeLength);
        for (int j = 0; j < chromosomeLength; ++j) {
            population[i].chromosome[j] = genePool[randomInt(0, genePool.size() - 1)];
        }
    }
    // Don't evaluate fitness here, let the evolve loop do it for the first time.
}

void GeneticAlgorithm::evolve(Config &design_configuration, CellDesign &processing_cell_design) {
    // 1. Evaluate the fitness of the current population
    evaluatePopulation(design_configuration, processing_cell_design);

    // 2. Create the next generation
    std::vector<Individual> nextGeneration;
    nextGeneration.reserve(populationSize);

    // Elitism: Keep the best individual from the current generation
    Individual best = getBestIndividual();
    nextGeneration.push_back(best);

    // 3. Generate the rest of the new population through selection, crossover, and mutation
    while (nextGeneration.size() < populationSize) {
        // Select two parents
        std::vector<Individual> parents = selectParents();
        
        std::vector<int> offspring_chromosome;
        // Apply crossover if random chance passes
        if (randomDouble() < crossoverRate && parents.size() >= 2) {
            offspring_chromosome = crossover(parents[0].chromosome, parents[1].chromosome);
        } else {
            // Otherwise, just clone the first parent
            offspring_chromosome = parents[0].chromosome;
        }

        // Apply mutation
        mutate(offspring_chromosome);
        
        Individual offspring;
        offspring.chromosome = offspring_chromosome;
        nextGeneration.push_back(offspring);
    }

    // 4. Replace the old population with the new generation
    population = nextGeneration;
    generationCount++;
}

Individual GeneticAlgorithm::getBestIndividual() const {
    if (population.empty()) {
        return Individual(); // Return an empty individual
    }
    // Find the individual with the maximum fitness
    auto best = std::max_element(population.begin(), population.end(),
                                 [](const Individual& a, const Individual& b) {
                                     return a.fitness < b.fitness;
                                 });
    return *best;
}

int GeneticAlgorithm::getGenerationCount() const {
    return generationCount;
}

// ==============================================================================
// PRIVATE HELPER METHODS (GA OPERATORS)
// ==============================================================================

void GeneticAlgorithm::evaluatePopulation(Config &design_configuration, CellDesign &processing_cell_design) {
    for (auto& individual : population) {
        individual.fitness = fitnessFunction(individual.chromosome, design_configuration, processing_cell_design);
    }
}

std::vector<Individual> GeneticAlgorithm::selectParents() {
    // Implementing Tournament Selection - it's generally more robust than roulette wheel.
    // TODO: You can implement other selection methods like roulette wheel later if needed.
    
    std::vector<Individual> parents;
    parents.reserve(2);
    
    int tournament_size = 5; // A common choice for tournament size
    
    for (int i=0; i<2; ++i) { // Select 2 parents
        Individual best_in_tournament;
        for (int j = 0; j < tournament_size; ++j) {
            int randomIndex = randomInt(0, population.size() - 1);
            if (j == 0 || population[randomIndex].fitness > best_in_tournament.fitness) {
                best_in_tournament = population[randomIndex];
            }
        }
        parents.push_back(best_in_tournament);
    }
    
    return parents;
}

std::vector<int> GeneticAlgorithm::crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    // Implementing Single-Point Crossover
    // TODO: You can implement other crossover methods like uniform crossover.
    if (parent1.empty()) return {};

    int crossoverPoint = randomInt(1, parent1.size() - 1);
    
    std::vector<int> child_chromosome;
    child_chromosome.reserve(parent1.size());
    
    // Take the first part from parent1
    child_chromosome.insert(child_chromosome.end(), parent1.begin(), parent1.begin() + crossoverPoint);
    // Take the second part from parent2
    child_chromosome.insert(child_chromosome.end(), parent2.begin() + crossoverPoint, parent2.end());
    
    return child_chromosome;
}

void GeneticAlgorithm::mutate(std::vector<int>& chromosome) {
    // TODO: This is where you implement the mutation logic.
    // For each gene in the chromosome, check if it should be mutated based on mutationRate.
    // If it should, replace it with a random gene from the genePool.
    for (size_t i = 0; i < chromosome.size(); ++i) {
        if (randomDouble() < mutationRate) {
            chromosome[i] = genePool[randomInt(0, genePool.size() - 1)];
        }
    }
} 