#pragma once

#include <vector>
#include <functional>
#include <iostream>

#include "../../PCC_Objects.h"
#include "../../PCC_Support_Functions.h"


/// fitness functions
double maximizingOnesFitness(const std::vector<int>& chromosome); // T

double shannonEntropyFitness(const std::vector<int>& chromosome); // Sc

//double IrradiationDamageFitness(const std::vector<int>& chromosome); // Ic
double IrradiationDamageFitness(const std::vector<int>& chromosome, Config &design_configuration, CellDesign &processing_cell_design); // Ic

/// Helpful additions
double randomDouble();

int randomInt(int min, int max);

// 用于封装个体（染色体）及其适应度的结构体
// The 'Individual' struct bundles a chromosome (the state vector) with its fitness score.
struct Individual {
    std::vector<int> chromosome;
    double fitness;

    // 添加一个默认构造函数以简化初始化
    Individual() : fitness(0.0) {}
};

/**
 * @class GeneticAlgorithm
 * @brief A generic and reusable implementation of a Genetic Algorithm.
 *
 * This class encapsulates the entire GA process, including population management,
 * evolution through selection, crossover, and mutation. It is designed to be
 * highly configurable, especially with a pluggable fitness function.
 */
class GeneticAlgorithm {
public:
    /**
     * @brief Constructs a GeneticAlgorithm instance.
     * @param populationSize The number of individuals in the population.
     * @param mutationRate The probability of a gene mutating (0.0 to 1.0).
     * @param crossoverRate The probability of two parents creating offspring (0.0 to 1.0).
     * @param fitnessFunc A function that takes a chromosome and returns its fitness score.
     */
//    GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, std::function<double(const std::vector<int>&)> fitnessFunc);

    GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate,
                     std::function<double(const std::vector<int>&, Config&, CellDesign&)> fitnessFunc);
    /**
     * @brief Initializes the population with random individuals.
     * @param chromosomeLength The length of the chromosome (state vector) for each individual.
     * @param possibleGenes A vector of possible values for each gene (e.g., {0, 1}).
     */
    void initializePopulation(int chromosomeLength, const std::vector<int>& possibleGenes);

    /**
     * @brief Runs a single generation of the evolutionary process.
     * This involves evaluation, selection, crossover, and mutation.
     */
    void evolve(Config &design_configuration, CellDesign &processing_cell_design);

    // --- Getter Methods ---

    /**
     * @brief Retrieves the best individual from the current population.
     * @return The individual with the highest fitness score.
     */
    Individual getBestIndividual() const;

    /**
     * @brief Gets the current generation number.
     * @return The current generation count.
     */
    int getGenerationCount() const;

private:
    // --- Internal State ---
    std::vector<Individual> population;
    double mutationRate;
    double crossoverRate;
    int populationSize;
    int generationCount;
    std::function<double(const std::vector<int>&, Config&, CellDesign&)> fitnessFunction;
    std::vector<int> genePool; // Stores the possible genes for mutation

    // --- Internal GA Operators ---

    /**
     * @brief Calculates the fitness for each individual in the population.
     */
    void evaluatePopulation(Config &design_configuration, CellDesign &processing_cell_design);

    /**
     * @brief Selects parents for the next generation using tournament selection.
     * @return A vector of selected individuals (parents).
     */
    std::vector<Individual> selectParents();

    /**
     * @brief Creates an offspring from two parents using single-point crossover.
     * @param parent1 The first parent's chromosome.
     * @param parent2 The second parent's chromosome.
     * @return The child's chromosome.
     */
    std::vector<int> crossover(const std::vector<int>& parent1, const std::vector<int>& parent2);

    /**
     * @brief Applies mutation to a chromosome based on the mutation rate.
     * @param chromosome The chromosome to mutate.
     */
    void mutate(std::vector<int>& chromosome);
}; 