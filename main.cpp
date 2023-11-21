#include <iostream>
#include <cmath>
#include <vector>
#include "ComputeFunctionNode.hpp"
#include "Utils.hpp"


using namespace std;

#define Nx 50
#define Ny 50


int index(int i, int j) {
	return i * Ny + j;
}


double g(double x,double y){  //function boundary conditions
	return 1.0;
}


double f(double x,double y){  //main function
	return 1.0;
}


void setup_solution(Mesh &M, vector<double> &U, ComputeFunctionNode &b) {  //initialize function U to border value
	for (int i = 0; i < M.getDimensionY(); ++i) {
		for (int j = 0; j < M.getDimensionX(); ++j) {
			U[i * M.getDimensionY() + j] = b.getValue(i, j);
		}
	}
}

void jacobi_iteration(double h_square, vector<double> &U_old, vector<double> &U, ComputeFunctionNode &function,Mesh m){  //1 cycle of Jacobi
	double b;

	for (int i = 1; i < m.getDimensionY()-1; ++i) {
		for (int j = 1; j < m.getDimensionX()-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				
			b = function.getValue(i, j) * h_square;

			U[index(i,j)] =
				( b
				+ U_old[index(i+1,j)]
				+ U_old[index(i-1,j)]
				+ U_old[index(i,j+1)]
				+ U_old[index(i,j-1)] ) / 4.0;
		}
	}
}


double compute_residual(double h_square, vector<double> &U, vector<double> &residual, ComputeFunctionNode &function){  
	for (int i = 1; i < Ny-1; ++i) {
		for (int j = 1; j < Nx-1; ++j) { 
			residual[index(i,j)] =
				 h_square*function.getValue(i,j)
				+U[index(i,j-1)]
				+U[index(i-1,j)]
				+U[index(i+1,j)]
				+U[index(i,j+1)]
				-4*U[index(i,j)];
		}
	}

	return norm(residual);
}


int main (int argc, char *argv[]) {
	Mesh mesh(-1.0, 0.0, 2.0, 2.0, Nx,Ny);

	ComputeFunctionNode bordo(&mesh,g);
	ComputeFunctionNode funzione(&mesh,f);


	vector<double> U_old(mesh.getDimensionX()*mesh.getDimensionY());
	vector<double> U(mesh.getDimensionX()*mesh.getDimensionY());
	
	vector<double> residual(mesh.getDimensionX()*mesh.getDimensionY());
	double residual_norm=100;

	int count =0;

	setup_solution(mesh, U_old, bordo);
	setup_solution(mesh, U    , bordo);

	double h_2dimension = mesh.getDiscretizationStepX()*mesh.getDiscretizationStepY();


	while (residual_norm > 1e-6) {
		// Jacobi iteration
		jacobi_iteration(h_2dimension,U_old,U,funzione,mesh);
		residual_norm=compute_residual(h_2dimension,U,residual,funzione);
		swap(U, U_old);
	}


	export_to_matlab("U", U);


	return 0;
}
