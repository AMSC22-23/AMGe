#include <iostream>
#include "Lattice.hpp"
#include "Multigrid.hpp"


int main() {
	Lattice mesh(0.0, 0.0, 1.0, 1.0, 129);
	Multigrid solver(mesh, 3, 3, 2);


	solver.test_allocation();


	return 0;
}
