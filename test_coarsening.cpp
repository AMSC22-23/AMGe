#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"
#include "Smoothers.hpp"


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


	std::vector<double> u(fine.numel());
	std::vector<double> b(fine.numel());
	std::vector<double> r(fine.numel());


	std::vector<double> r_interpolated(fine.numel());
	std::vector<double> r_coarse(coarse.numel());


	fine.evaluate_boundary_conditions(u, g);
	fine.evaluate_forcing_term(b, f);


	residual(fine, u, b, r);
	fine.print_vector(r, "residual");

	fine.project_on_coarse(coarse, r, r_coarse);
	coarse.print_vector(r_coarse, "residual_coarse");

	fine.interpolate_on_fine(coarse, r_interpolated, r_coarse);
	fine.print_vector(r_interpolated, "r_interpolated");


	return 0;
}
