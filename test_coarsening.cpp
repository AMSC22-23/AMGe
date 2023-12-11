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


void set_initial_guess(Lattice &mesh, std::vector<double> &u, double (*g)(double x, double y)) {
	mesh.evaluate_function(u, g);

	for (Index i : mesh.get_inner_nodes()) {
		u[i] = 0.0;
	}
}


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}



int main (int argc, char *argv[]) {
	const int N = 5;


	Lattice fine(0.0, 0.0, 1.0, 1.0, N);
	Lattice coarse = fine.build_coarse();


	std::vector<double> u_fine(fine.numel());
	std::vector<double> r(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r_interpolated(fine.numel());
	std::vector<double> r_coarse(coarse.numel());


	set_initial_guess(fine, u_fine, g);
	fine.evaluate_forcing_term(b, f);
	for (Index i : fine.get_inner_nodes()) {
		u_fine[i] = 0.0;
	}


	residual(fine,u_fine, b, r);
	fine.print_vector(r, "residual");

	fine.project_on_coarse(coarse, r, r_coarse);
	coarse.print_vector(r_coarse, "residual_coarse");

	fine.interpolate_on_fine(coarse, r_interpolated, r_coarse);
	fine.print_vector(r_interpolated, "r_interpolated");


	return 0;
}
