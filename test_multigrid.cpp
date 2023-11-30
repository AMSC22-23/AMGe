#include <iostream>
#include <cmath>
#include <vector>
#include "LatticeMesh.hpp"
#include "Utils.hpp"


using namespace std;



double g(double x,double y){  //function boundary conditions
	return ((x*x*x)/6.0)+((y*y*y)/6.0);
}


double f(double x,double y){  //main function
	return -x-y;
}



double compute_error(const std::vector<double> &exact_solution, const std::vector<double> &computed_solution, int dim){
	std::vector<double> error(dim);
	for(Index i = 0; i < dim; i++){
		error[i]=exact_solution[i]-computed_solution[i];
	}
	return norm(error);
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


void gauss_seidel(LatticeMesh mesh, std::vector<double> &U, std::vector<double> &b){  //1 cycle of Gauss Siedel
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


int main(int argc, char *argv[]){
	const int Nx = 129;
	const int Ny = 129;


	LatticeMesh fine(0.0, 0.0, 1.0, 1.0, Nx,Ny);
	LatticeMesh coarse = fine.build_coarse();


	std::vector<double> U_exact(Nx*Ny);


	std::vector<double> U_fine(Nx*Ny);
	std::vector<double> residual_fine(Nx*Ny);
	std::vector<double> err_fine(Nx*Ny);
	std::vector<double> err_coarse((Nx / 2 + 1) * (Ny / 2 + 1));
	std::vector<double> residual_coarse((Nx / 2 + 1) * (Ny / 2 + 1));
	std::vector<double> F(Nx*Ny);


	fine.evaluate_function(U_exact, g);
	fine.evaluate_function(U_fine, g);
	fine.evaluate_function(F, f);


	for (Index i : fine.get_inner_nodes()) {
		U_fine[i] = 0.0;
	}


	/* iniziamo con il multigrid */
	for (int it = 0; it < 10000; ++it) {
		gauss_seidel(fine, U_fine, F);
		gauss_seidel(fine, U_fine, F);
		gauss_seidel(fine, U_fine, F);
		compute_residual(fine, U_fine, residual_fine, F);


		fine.project_on_coarse(residual_fine, residual_coarse);


		// bisogna inizializzare con una buona stima iniziale U_coarse, chiedere
		fine.project_on_coarse(residual_fine, err_coarse);
		gauss_seidel(coarse, err_coarse, residual_coarse);
		gauss_seidel(coarse, err_coarse, residual_coarse);
		gauss_seidel(coarse, err_coarse, residual_coarse);
		fine.interpolate_on_fine(coarse, err_fine, err_coarse);


		for (Index i : fine.get_inner_nodes()) {
			U_fine[i] += err_fine[i];
		}


		gauss_seidel(fine, U_fine, F);
		gauss_seidel(fine, U_fine, F);
		gauss_seidel(fine, U_fine, F);
	}


	std::cout << compute_error(U_fine, U_exact, Nx*Ny);
}

