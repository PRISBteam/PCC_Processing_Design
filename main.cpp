///******************************************************************************************************************************///
///*                                                                                                                            *///
///******************************************************************************************************************************///
///**********************************   Elijah Borodin, Manchester, Spring 2022   ***********************************************///
///******************************************************************************************************************************///
///****************************************   DCCAnalyser(c) utility   ********************************************************///
///******************************************************************************************************************************///

/// Standard (STL) C++ libraries:
///------------------------------
#include <iostream>
#include <ctime>
///-------------------------------

/// Attached user defined C++ libraries:
///-------------------------------------
#include "src/DCCProcessing.h" // Change element types of the DCC itself
#include "src/DCCKinetic.h" // Generate a process on the elements of the DCC without any changes in the complex

///*.....................................................................    Main    .............................................*///
int main() {
/** Two functions (dependent on the problem's spatial dimension - 2D or 3D) are launching here with the arguments of the paths for
 * adjacency AN, AE, AF, (AG) and incidence (boundary operators) BEN, BFE, (BGF) matrices. **/

    using namespace std; ///Standard namespace

    cout << "=============================================" << endl << "Start of the DCC Processing" << endl << "-------------------------------------------------------------" << endl;

/// Space dimension of the problem
    int dim = 3;
  /**  do { ///Manual user input of the space dimension value (dim)
        cout << " Please, input the dimension of the problem: 3 for 3D and 2 for 2D "s << endl;
        cin >> dim;
        cout << "The problem dimension is\t" << dim << endl;
        if (dim != 2 && dim != 3) cout << "Input Error: Please retype 2 or 3" << endl;
    }while (dim != 2 && dim != 3); **/

/// Please, input below the correct source directory and simulation type:
    char simulation_type = 'R'; // 'R', 'S', 'F' or 'I', 'E' :: This char define the process type: 'R' for Random, 'S' for maximum configuration Entropy production, 'F' for the 3D one-layer film, 'I' for the Ising-like model, 'E' for experimental data obtained by EBSD ///
    /**
     do { ///Manual user input of the simulation type
          cout << " Please, input the symbol of particular simulation type of the problem: E (experiment), R (rndom), S (entropy maximisation) or I (Ising model):"s << endl;
          cin >> simulation_type;
          cout <<"The simulation type is\t"<< simulation_type << endl;
          if (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I') cout << "Input Error: Please retype 'E', 'R', 'S' or 'I' for the specific simulation type" << endl;
      }while (simulation_type != 'E' && simulation_type != 'R' && simulation_type != 'S' && simulation_type != 'I');
   **/

/// Here the names of the directories with the source files "srcdir" and output files "outdir" must be defined
    string srcdir = "/Users/user/Dropbox/OFFICE/NEPER/resultsApril2022/"s, problem_folder_path = "1k3cells/"s;
    if (simulation_type == 'F') string srcdir = "/Users/user/Dropbox/OFFICE/NEPER/resultsApril2022/2D_film_seeds/"s, problem_folder_path = "2D_10k/"s;
/**    /// Manual user input of the DCC files folder path
        cout << " Please, input the name of the folder where the DCC source files are (like 1k3cells/ ):"s << endl;
        cin >> problem_folder_path;
**/
    srcdir = srcdir + problem_folder_path;
    cout << srcdir << endl;

    string outdir = "/Users/user/Dropbox/OFFICE/CProjects/Voro3D/"s, results_output_folder_path = "test/"s;
    if (simulation_type == 'R') string outdir = "/Users/user/Dropbox/OFFICE/CProjects/Voro3D_film/"s, results_output_folder_path = "test/"s;
    if (dim == 2) string outdir = "/Users/user/Dropbox/OFFICE/CProjects/Voro2D/", results_output_folder_path = "test/"s;
/**    /// Manual user input of the simulation results output folder path
        cout << " Please, input the simulation results output folder path (like test/ ):"s << endl;
        cin >> results_output_folder_path;
**/

    outdir = outdir + results_output_folder_path;
    cout << outdir << endl;

    char* output_folder = const_cast<char*>(outdir.c_str()); // From string to char for the passing folder path to a function

/// Below the file names with the sparse DCC matrices must be defined
    string ssd0 = srcdir + "A0.out"s, ssd1 = srcdir + "A1.out"s, ssd2 = srcdir + "A2.out"s, ssd3 = srcdir + "A3.out"s, ssd4 = srcdir + "B10.out"s, ssd5 = srcdir + "B21.out"s, ssd6 = srcdir + "B32.out"s;
    //The next line just a technical procedure string to char arrays transformation needed to use them as the function arguments
    char* sd0 = const_cast<char*>(ssd0.c_str()); char *sd1 = const_cast<char*>(ssd1.c_str()); char *sd2 = const_cast<char*>(ssd2.c_str()); char *sd3 = const_cast<char*>(ssd3.c_str()); char *sd4 = const_cast<char*>(ssd4.c_str()); char *sd5 = const_cast<char*>(ssd5.c_str()); char *sd6 = const_cast<char*>(ssd6.c_str());
/// File path with the amounts of the different cells  (1st line for Nodes, 2nd line for Edges, 3rd line for Faces and (in 3D) 4th line for grains)
    string  ncells = srcdir + "number_of_cells.out"s; char* number_of_cells = const_cast<char*>(ncells.c_str());
///File path with the configuration profile (types of special faces and calculating parameters)
    string  config = srcdir + "config.txt"s; char* configuration = const_cast<char*>(config.c_str());

    if (dim == 3)
        if (simulation_type == 'F')
            HAGBsKinetic3D(sd0, sd1, sd2, sd3, sd4, sd5, sd6, number_of_cells, configuration, output_folder, simulation_type);
        else
            HAGBsProbability3D(sd0, sd1, sd2, sd3, sd4, sd5, sd6, number_of_cells, configuration, output_folder, simulation_type);
    else if (dim == 2)
        HAGBsProbability2D(sd0, sd1, sd2, sd4, sd5, number_of_cells, configuration, output_folder, simulation_type);

/// Elapsing time
unsigned int end_time = clock();
double fulltime = (double) end_time/ 1000.0;
cout << "HAGBsProbability " << dim << "D " << "runtime is equal to  " << fulltime <<  "  seconds" << endl;

cout << "-------------------------------------------------------------" << endl << "The end of the VoroCAnalyser program" << endl << "=============================================" << endl ;

return 0;
}/// The end of Main function
