#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include "LatticeMesh.hpp"
#include "Utils.hpp"


using namespace std;


double compute_error(const std::vector<double> &exact_solution, const std::vector<double> &computed_solution, int dim){
	std::vector<double> error(dim);
	for(Index i = 0; i < dim; i++){
		error[i]=exact_solution[i]-computed_solution[i];
	}

	return norm(error);
	//return *std::max_element(error.begin(), error.end());
}


double compute_residual(LatticeMesh mesh, vector<double> &u, vector<double> &r, vector<double> &b){    //compute residual vecotr
	const double h = mesh.hx;

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		r[i] = h*h*b[i] + u[nord] + u[sud] + u[ovest] + u[est] - 4.0 * u[i];
	}

	return norm(r);
}


void gauss_seidel(LatticeMesh &mesh, std::vector<double> &u, std::vector<double> &b) {
	const double h = mesh.hx;

	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);

		u[i] = 0.25 * (h*h*b[i] + u[nord] + u[sud] + u[ovest] + u[est]);
	}
}


double g(double x,double y){
	return std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
}


double f(double x,double y){
	return 8.0 * M_PI * M_PI * std::cos(2.0 * M_PI * x) * std::cos(2.0 * M_PI * y);
}

int main(int argc, char *argv[]){
	for (int lvl = 0; lvl < 8; ++lvl) {
		int N = (2 << lvl) + 1;
		LatticeMesh mesh(0.0, 0.0, 1.0, 1.0, N, N);

		std::vector<double> U_exact(N*N);
		std::vector<double> residual(N*N);
		std::vector<double> U_computed(N*N);
		std::vector<double> F(N*N);


		mesh.evaluate_function(U_exact, g);
		mesh.evaluate_function(U_computed, g);
		mesh.evaluate_function(F, f);


		for (Index i : mesh.get_inner_nodes()) {
			U_computed[i] = 0.0;
		}


		int it = 0;
		do {
			gauss_seidel(mesh, U_computed, F);
			++it;
		} while (compute_residual(mesh, U_computed, residual, F) > 1e-14);


		// here error should be diminish with N but it is not
		std::cout << "iterazioni:" << it << std::endl;
		std::cout << "N:         " << N << std::endl;
		std::cout << "errore:    " << compute_error(U_exact, U_computed, N*N) << std::endl;
		std::cout << "================" << std::endl;
	}


	return 0;
}
