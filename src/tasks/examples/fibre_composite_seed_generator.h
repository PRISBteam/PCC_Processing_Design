#ifndef PCC_PROCESSING_DESIGN_FIBRE_COMPOSITE_SEED_GENERATOR_H
#define PCC_PROCESSING_DESIGN_FIBRE_COMPOSITE_SEED_GENERATOR_H

unsigned int Mersenne_Int_Rand(unsigned int OCellsNumb); // Random generation machine for a new 2-Cell number
double Mersenne_Double_Rand(double DCellsNumb); // Random generation machine for a new 2-Cell number

/*!
 * @brief
 * @param number_of_fibres
 * @param fibre_diameter_fraction
 * @return
 */
 std::set<std::tuple<double,double,double>> Fibre_composite_seed_generator(unsigned int number_of_grains, int number_of_fibres, double fibre_diameter_fraction);

#endif //PCC_PROCESSING_DESIGN_FIBRE_COMPOSITE_SEED_GENERATOR_H
