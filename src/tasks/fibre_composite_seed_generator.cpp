/// Fibre Composite Grain Seeds Generator (for Neper) task 2024
/// Dr Elijah Borodin

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <set>
#include <vector>
#include <map>
#include <random>

#include "fibre_composite_seed_generator.h"

using namespace std;
/*!
 * @details
 * @param number_of_fibres
 * @param fibre_diameter_fraction
 * @return
 */
std::set<std::tuple<double,double,double>> Fibre_composite_seed_generator(unsigned int number_of_grains, int number_of_fibres, double fibre_diameter_fraction) {
    std::set<std::tuple<double,double,double>> FC_grain_seeds; // function output - a set
    std::map<int,tuple<double,double>> fibres_to_plane_coordinates;

    /// Initial conditions output
    std::cout << "\tNumber of fibres\t:\t" << number_of_fibres << std::endl;
    std::cout << "\tNumber of grains\t:\t" << number_of_grains << std::endl;
    std::cout << "\tFibre diameter as the fraction of total length\t:\t" << fibre_diameter_fraction << std::endl;

    /// 9 possible fibres:
    //middle row
    fibres_to_plane_coordinates.insert({0,make_tuple(0.43,0.5)});
    fibres_to_plane_coordinates.insert({1,make_tuple(0.5,0.5)});
    fibres_to_plane_coordinates.insert({2,make_tuple(0.63,0.5)});
    //upper row
    fibres_to_plane_coordinates.insert({3,make_tuple(0.4,0.6)});
    fibres_to_plane_coordinates.insert({4,make_tuple(0.5,0.63)});
    fibres_to_plane_coordinates.insert({5,make_tuple(0.6,0.6)});
    //bottom row
    fibres_to_plane_coordinates.insert({6,make_tuple(0.4,0.4)});
    fibres_to_plane_coordinates.insert({7,make_tuple(0.5,0.43)});
    fibres_to_plane_coordinates.insert({8,make_tuple(0.6,0.4)});

// REPAIR    cout << " sin30 = " << sin(30.0/57.2958) << endl; exit(0);

    double x_coord =0, y_coord =0, z_coord =0;
    for (unsigned int ng = 0; ng < number_of_grains ; ++ng) {
        if (ng <= 0.01*number_of_grains) { // 1/7 of grains are distributed randomly over the space
            x_coord = Mersenne_Double_Rand(1);;
            y_coord = Mersenne_Double_Rand(1);;
            z_coord = Mersenne_Double_Rand(1);;
        }
        else {
            int f_numb = Mersenne_Int_Rand(number_of_fibres);
            double f_rad = Mersenne_Double_Rand(fibre_diameter_fraction); // from 0 to fibre radius
            double f_angle = Mersenne_Double_Rand(360); // angle in degrees
            double f_z = Mersenne_Double_Rand(1); // from 0 to fibre radius

            x_coord = std::get<0>(fibres_to_plane_coordinates.at(f_numb)) + f_rad * cos(f_angle / 57.2958);
            y_coord = std::get<1>(fibres_to_plane_coordinates.at(f_numb)) + f_rad * sin(f_angle / 57.2958);
            z_coord = f_z;
        } // end of if (ng <= 0.1*number_of_grains) .. else {}
            FC_grain_seeds.insert(make_tuple(x_coord,y_coord,z_coord));
    }
    // ALL fibres oriented along OZ-axis (!)
    // Let us assume 3D cube 1.0 x 1.0 x 1.0

// RAND    std::cout << "\tGrain seeds\t:\t" << std::endl; for (auto gs : FC_grain_seeds) std::cout << std::get<0>(gs) << "\t\t" << std::get<1>(gs) << "\t\t" << std::get<2>(gs) << std::endl;

    return FC_grain_seeds;

} // end of Fibre_composite_seed_generator()

unsigned int Mersenne_Int_Rand(unsigned int OCellsNumb){ // Random generation machine for a new 2-Cell number
    std::random_device rd; // seed for a device generating unsigned random integers
    std::mt19937 mt(rd()); // advanced random engine based on the Mersenne Twister 19937 algorithm proposed in [M. Matsumoto and T. Nishimura, ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, Pages 3–30, https://dl.acm.org/doi/pdf/10.1145/272991.272995]
    uniform_int_distribution<unsigned int> uni_rand (0, OCellsNumb - 1); // uniformly distributed from 0 to OCellsNumb-1 inclusive

    return uni_rand(mt); // random generation of the boundary number in the range from 0 to OrdinaryCellNumbs.size()-1
} // END of Mersenne_Int_Rand(unsigned int OCellsNumb)

double Mersenne_Double_Rand(double DCellsNumb){ // Random generation machine for a new 2-Cell number
    std::random_device drd; // seed for a device generating unsigned random integers
    std::mt19937 mtd(drd()); // advanced random engine based on the Mersenne Twister 19937 algorithm proposed in [M. Matsumoto and T. Nishimura, ACM Transactions on Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, Pages 3–30, https://dl.acm.org/doi/pdf/10.1145/272991.272995]
    std::uniform_real_distribution<double> dis (0, DCellsNumb); // uniformly distributed double from 0 to DCellsNumb inclusive

    return dis (mtd); // random generation of the boundary number in the range from 0 to OrdinaryCellNumbs.size()-1
} // END of Mersenne_Double_Rand(double DCellsNumb)

/*
// Pattern for main.cpp TASK mode::

int number_of_fibres = 9;
unsigned int number_of_grains = 100000;
double fibre_diameter_fraction = 0.05;
#include "tasks/fibre_composite_seed_generator.h"
std::set<std::tuple<double,double,double>> grain_seeds_set = Fibre_composite_seed_generator(number_of_grains, number_of_fibres, fibre_diameter_fraction);

std::ofstream Out_task_data; // performance_test output stream to the output file
Out_task_data.open(output_dir + "Out_task_data_100k_9f_b5.txt"s, ios::trunc); // creates 'performance_test.txt' file in the 'output_dir'
for (auto gss : grain_seeds_set)
Out_task_data << std::get<0>(gss) << "\t" << std::get<1>(gss) << "\t" << std::get<2>(gss) << endl;
Out_task_data.close();
*/