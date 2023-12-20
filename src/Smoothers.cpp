#include "Smoothers.hpp"


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


void jacobi(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}
}


void jacobi_parallel_naive(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	// implementazione naive, forse openmp non riesce a parallelizzare i for di c++11
	#pragma omp parallel for

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}

	#pragma omp barrier
}


void jacobi_parallel(Lattice &mesh, std::vector<double> &u, const std::vector<double> &old, const std::vector<double> &b) {
	// adattato il for loop alla sua forma tradizionale
	const auto inner_nodes = mesh.get_inner_nodes();

	#pragma omp parallel for

	for (int j = 0; j < inner_nodes.size() ; j ++) {
		auto i = inner_nodes[j];

		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (b[i] + old[nord] + old[sud] + old[ovest] + old[est]);
	}

	#pragma omp barrier
}
