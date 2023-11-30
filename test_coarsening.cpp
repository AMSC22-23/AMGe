#include <iostream>
#include <cmath>
#include <vector>
#include "LatticeMesh.hpp"
#include "Utils.hpp"


using namespace std;



int main (int argc, char *argv[]) {
	const int Nx = 7;
	const int Ny = 7;


	LatticeMesh mesh(-1.0, 0.0, 2.0, 2.0, Nx,Ny);
	LatticeMesh coarse = mesh.build_coarse();


	std::vector<double> u_fine(Nx*Ny);
	std::vector<double> u_coarse((1+Nx/2) * (1+Ny/2));


	for (int i = 0; i < u_fine.size(); ++i) {
		u_fine[i] = i;
	}


	mesh.project_on_coarse(u_fine, u_coarse);


	for (auto x : u_coarse) {
		std::cout << x << std::endl;
	}


	return 0;
}
