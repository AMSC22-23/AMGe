#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"


using namespace std;


void residual(Lattice &mesh, const std::vector<double> &u, const std::vector<double> &b, std::vector<double> &r) {
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


void set_initial_guess(Lattice &mesh, std::vector<double> &u, double (*g)(double x, double y)) {
	mesh.evaluate_function(u, g);

	for (Index i : mesh.get_inner_nodes()) {
		u[i] = 0.0;
	}
}

void two_level(Lattice &fine, Lattice &coarse, std::vector<double> &u, std::vector<double> &b){
	std::vector<double> r_fine(fine.numel());
	std::vector<double> e_fine(fine.numel());
	std::vector<double> r_coarse(coarse.numel());
	std::vector<double> e_coarse(coarse.numel());

	for (int pre = 0; pre < 10; ++pre) {
			gseidel(fine, u, b);
		}
	
		residual(fine, u, b, r_fine);

		fine.project_on_coarse(coarse, r_fine, r_coarse);


		for (auto &x : e_coarse) {
			x = 0.0;
		}

		for (int coarse_it = 0; coarse_it < 300; ++coarse_it) {
			gseidel(coarse, e_coarse, r_coarse);
		}

		fine.interpolate_on_fine(coarse, e_fine, e_coarse);

		for (Index i : fine.get_inner_nodes()) {
			u[i] -= e_fine[i];
		}

		for (int post = 0; post < 10; ++post) {
			gseidel(fine, u, b);
		}
}


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


int main (int argc, char *argv[]) {
	const int N = 65;
	
	
	Lattice fine(0.0, 0.0, 1.0, 1.0, N);
	Lattice coarse = fine.build_coarse();


	std::vector<double> u_fine(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r_fine(fine.numel());
	std::vector<double> e_fine(fine.numel());


	std::vector<double> r_coarse(coarse.numel());
	std::vector<double> e_coarse(coarse.numel());


	set_initial_guess(fine, u_fine, g);
	fine.evaluate_forcing_term(b, f);
	
	for (Index i : fine.get_inner_nodes()) {
		u_fine[i] = 0.0;
	}



	for (int i = 0; i < 1000; ++i) {
		two_level(fine, coarse, u_fine, b);
		residual(fine, u_fine, b, r_fine);
		std::cout << norm(r_fine) << std::endl;
	}



	return 0;
}
