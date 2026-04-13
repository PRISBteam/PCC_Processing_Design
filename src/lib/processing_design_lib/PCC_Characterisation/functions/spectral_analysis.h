#ifndef PCC_PROCESSING_DESIGN_SPECTRAL_ANALYSIS_H
#define PCC_PROCESSING_DESIGN_SPECTRAL_ANALYSIS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace std; // standard namespace
using namespace Eigen; // Eigen library namespace
//using namespace Spectra; // Spectra library namespace

typedef Eigen::SparseMatrix<double> SpMat; // <Eigen> library class, which declares a column-major sparse matrix type of doubles with the nickname 'SpMat'


std::vector<Eigen::SparseMatrix<double>> ReducedFaceLaplacians(std::vector<unsigned int> &face_sequence_vector, SpMat &ENS, SpMat &FES, SpMat &GFS);

std::vector<Eigen::SparseMatrix<double>> ReducedFaceLaplacians(std::vector<unsigned int> &face_sequence_vector, SpMat &ENS, SpMat &FES, SpMat &GFS);

std::vector<double> OperatorSpectrum(Eigen::SparseMatrix<double> const &OSM);

Eigen::MatrixXcd OperatorEigvectors(Eigen::SparseMatrix<double> const &OSM);

std::vector<double> OperatorsBetti(std::vector<unsigned int> &face_sequence_vector, SpMat &ENS, SpMat &FES, SpMat &GFS);

#endif //PCC_PROCESSING_DESIGN_SPECTRAL_ANALYSIS_H
