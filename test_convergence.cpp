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



double compute_error(vector<double> exact_solution, vector<double> computed_solution,int dim){
	vector<double> error(dim);
	for(Index i = 0; i < dim; i++){
		error[i]=exact_solution[i]-computed_solution[i];
	}
	return norm(error);
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


int main(int argc, char *argv[]){
    const int Nx = 5;
    const int Ny = 5;


    LatticeMesh mesh(0.0, 0.0, 1.0, 1.0, Nx,Ny);
    
    vector<double> U(Nx*Ny);
	vector<double> F(Nx*Ny);


    mesh.evaluate_function(U, g);
	mesh.evaluate_function(F, f);

    vector<double> solution(Nx*Ny);

    mesh.evaluate_function(solution, g);
    for (int it = 0; it < 1000; ++it) {
		//gauss_seidel(mesh, U, F);
		gauss_seidel(mesh, U, F);
	}

    cout<<compute_error(solution,U,Nx*Ny);


}