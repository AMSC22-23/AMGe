#ifndef __SMOOTHERS_HPP__
#define __SMOOTHERS_HPP__


#include <vector>
#include "Lattice.hpp"


void residual(Lattice &mesh, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r) {
	// @TODO: add also parallelization here
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		r[i] = 4.0 * u[i] - u[nord] - u[sud] - u[ovest] - u[est] - b[i];
	}
}


void gseidel(Lattice &mesh, std::vector<double> &u, const std::vector<double> &b) {
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + u[nord] + u[sud] + u[ovest] + u[est]);
	}
}


void jacobi_parallel(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	// sembra che openmp non capisca bene i range based loop di c++11, quindi bisogna adattare il for a qualcosa di tradizionale (for(int i = 0; ...))
	#pragma omp parallel for

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}
	#pragma omp barrier
}


void jacobi(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}
}


#endif // !__SMOOTHERS_HPP__
