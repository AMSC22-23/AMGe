#include <cmath>
#include <iostream>
#include <vector>
#include "Lattice.hpp"
#include "Multigrid.hpp"
#include "Smoothers.hpp"


double g(double x,double y) {
	return x * x * x * x * x / 20.0 + y * y * y * y * y / 20.0;
}


double f(double x,double y) {
	return -(x * x * x + y * y * y);
}


int main() {
	Lattice mesh(0.0, 0.0, 1.0, 1.0, 129);

	std::vector<double> u(mesh.numel());
	std::vector<double> r(mesh.numel());
	std::vector<double> b(mesh.numel());


	mesh.evaluate_boundary_conditions(u, g);
	mesh.evaluate_forcing_term(b, f);


	Multigrid solver(mesh, 3, 3, 2);


	for (int i = 0; i < 100; ++i) {
		solver.step(u, b);
		residual(mesh, u, b, r);

		std::cout << norm(r) << std::endl;
	}


	return 0;
}
