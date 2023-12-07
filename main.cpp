#include <iostream>
#include <cmath>
#include <vector>
#include "Lattice.hpp"
#include "Utils.hpp"


using namespace std;

/*

double g(double x,double y){  //function boundary conditions
	return 1.0;
}


double f(double x,double y){  //main function
	return 1.0;
}




void gauss_seidel(LatticeMesh mesh, vector<double> &U, vector<double> &b){  //1 cycle of Gauss Siedel
	const double hx_square = mesh.hx * mesh.hx;
	const double hy_square = mesh.hy * mesh.hy;
	const double den       = 2.0 * ((1.0 / hx_square) + (1.0 / hy_square));
	

	for(Index i : mesh.get_inner_nodes()){
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);
		
		
		U[i] = (
			  b[i]
			+ ((U[ovest] + U[est]) / hx_square)
			+ ((U[nord]  + U[sud]) / hy_square)
		) / den;
	}

}


double compute_residual(LatticeMesh mesh, vector<double> &U, vector<double> &residual, vector<double> &b){    //compute residual vecotr
	double hx_square= mesh.hx*mesh.hx;
	double hy_square= mesh.hy*mesh.hy;
	double factor =-2.0*(hx_square+hy_square);
	

	for(Index i : mesh.get_inner_nodes()){
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);
		
		
		residual[i] = (
						hx_square*hy_square*b[i]
					   +((U[nord]+U[sud])*hy_square)
					   +((U[ovest]+U[est])*hx_square)
					   +factor*U[i]
					   );
	}

	return norm(residual);
}



void jacobi(LatticeMesh &mesh, std::vector<double> &U, std::vector<double> &Old, std::vector<double> &b){
	const double hx_square = mesh.hx * mesh.hx;
	const double hy_square = mesh.hy * mesh.hy;
	const double den       = 2.0 * ((1.0 / hx_square) + (1.0 / hy_square));


	for (Index i : mesh.get_inner_nodes()) {
		const auto [nord, sud, ovest, est] = mesh.get_cardinal_neighbours(i);


		U[i] = (
			  b[i]
			+ ((Old[ovest] + Old[est]) / hx_square)
			+ ((Old[nord]  + Old[sud]) / hy_square)
		) / den;
	}
}
*/


int main (int argc, char *argv[]) {
	Lattice mesh(0.0, 0.0, 1.0, 1.0, 5);

	/*
	const int Nx = 50;
	const int Ny = 50;


	LatticeMesh mesh(-1.0, 0.0, 2.0, 2.0, Nx,Ny);
	vector<double> Old(Nx*Ny);
	vector<double> U(Nx*Ny);
	vector<double> F(Nx*Ny);


	mesh.evaluate_function(U, g);
	mesh.evaluate_function(Old, g);
	mesh.evaluate_function(F, f);

	vector<double> residual(Nx*Ny);
	float residual_norm=100;
	int count=0;


	for (int it = 0; it < 1000; ++it) {
		gauss_seidel(mesh, U, F);
		//jacobi(mesh, U, Old, F);
		compute_residual(mesh,U,residual,F);
		count++;
		//std::swap(U, Old);
		residual_norm=norm(residual);
	}


	export_to_matlab("U", U);
	*/


	return 0;
}
