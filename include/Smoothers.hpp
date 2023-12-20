#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include <vector>
#include "Lattice.hpp"


void residual             (Lattice &mesh, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r);
void gseidel              (Lattice &mesh,       std::vector<double> &u, const std::vector<double> &b);
void jacobi               (Lattice &mesh,       std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b);
void jacobi_parallel_naive(Lattice &mesh,       std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b);
void jacobi_parallel      (Lattice &mesh,       std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b);


#endif // !__SMOOTHERS_HPP__
