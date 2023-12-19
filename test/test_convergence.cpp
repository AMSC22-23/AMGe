#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"


using namespace std;


double compute_error(const std::vector<double> &exact_solution, const std::vector<double> &computed_solution, int dim){
	std::vector<double> error(dim);
	for(Index i = 0; i < dim; i++){
		error[i]=exact_solution[i]-computed_solution[i];
	}

	return norm(error);
	//return *std::max_element(error.begin(), error.end());
}

void set_initial_guess(Lattice &mesh, std::vector<double> &u, double (*g)(double x, double y)) {
	mesh.evaluate_function(u, g);

	for (Index i : mesh.get_inner_nodes()) {
		u[i] = 0.0;
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
		Lattice mesh(0.0, 0.0, 1.0, 1.0, N);

		std::vector<double> U_exact(mesh.numel());
		std::vector<double> r(mesh.numel());
		std::vector<double> U_computed(mesh.numel());
		std::vector<double> F(mesh.numel());


		mesh.evaluate_function(U_exact, g);
		set_initial_guess(mesh, U_computed, g);
		mesh.evaluate_forcing_term(F, f);


		int it = 0;
		do {
			gseidel(mesh, U_computed, F);
			++it;
			residual(mesh, U_computed, F, r);
		} while (norm(r) > 1e-8);


		// here error should be diminish with N but it is not
		std::cout << "iterazioni:" << it << std::endl;
		std::cout << "N:         " << N << std::endl;
		std::cout << "errore:    " << compute_error(U_exact, U_computed, N*N) << std::endl;
		std::cout << "================" << std::endl;
	}


	return 0;
}
