#include <iostream>
#include <cmath>
#include <vector>
#include "ComputeFunctionNode.hpp"
#include "Utils.hpp"


using namespace std;

#define Nx 50
#define Ny 50
#define N_INTERNAL (N-2)


/*void print_matrix(double A[N_INTERNAL*N_INTERNAL][N_INTERNAL*N_INTERNAL]){
	for(int i = 0; i < N_INTERNAL*N_INTERNAL; ++i){
		for(int j = 0; j < N_INTERNAL*N_INTERNAL; ++j){
			std::cout << A[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


void print_vector(double b[N_INTERNAL*N_INTERNAL]) {
		for(int i = 0; i < N_INTERNAL*N_INTERNAL; ++i){
			std::cout << b[i] << std::endl;
		}
}*/


int index(int i, int j) {
	return i * Ny + j;
}


/*int local_index(int i, int j) {			//index matrix in a vector
	return (i-1) * N_INTERNAL + j - 1;
}*/



double g(double x,double y){  //function buondary conditions
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

void jacobi_iteration(vector<double> &U_old, vector<double> &U, ComputeFunctionNode &function){  //1 cycle of Jacobi
	double b;
	double hx_square= function.getMesh()->getDiscretizationStepX()*function.getMesh()->getDiscretizationStepX();
	double hy_square= function.getMesh()->getDiscretizationStepY()*function.getMesh()->getDiscretizationStepY();
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	

	for (int i = 1; i < function.getMesh()->getDimensionY()-1; ++i) {
			for (int j = 1; j < function.getMesh()->getDimensionX()-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function.getValue(i, j);
				U[index(i,j)] =
					( b
					+ ((U_old[index(i+1,j)]+U_old[index(i-1,j)])/hx_square)
					+ ((U_old[index(i,j+1)]+U_old[index(i,j-1)])/hy_square) ) / den;		
			}
	}
}


double compute_residual(vector<double> &U, vector<double> &residual, ComputeFunctionNode &function){    //compute residual vecotr
	double h_squarex= function.getMesh()->getDiscretizationStepX()*function.getMesh()->getDiscretizationStepX();
	double h_squarey= function.getMesh()->getDiscretizationStepY()*function.getMesh()->getDiscretizationStepY();
	double factor =-2.0*(h_squarex+h_squarey);
	
	for (int i = 1; i < function.getMesh()->getDimensionY()-1; ++i) {
		for (int j = 1; j < function.getMesh()->getDimensionX()-1; ++j) { 
			residual[index(i,j)] =
				 (h_squarex*h_squarey*function.getValue(i,j))
				+((U[index(i,j-1)]+U[index(i,j+1)])*h_squarex)
				+((U[index(i-1,j)]+U[index(i+1,j)])*h_squarey)
				+factor*U[index(i,j)];
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



	while(residual_norm>1.e-8){
		// Jacobi iteration
		jacobi_iteration(U_old,U,funzione);
		residual_norm=compute_residual(U,residual,funzione);
		swap(U, U_old);
		cout<<residual_norm<<endl;
	}

	/*for (auto x : U) {
		cout << x << endl;
	}*/

	cout<<residual_norm<<endl;


	return 0;
}
