#include <iostream>
#include <cmath>
#include <vector>
#include "ComputeFunctionNode.hpp"


#define N 50
#define N_INTERNAL (N-2)


void print_matrix(double A[N_INTERNAL*N_INTERNAL][N_INTERNAL*N_INTERNAL]){
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
}


int index(int i, int j) {
	return i * N + j;
}


int local_index(int i, int j) {			//index matrix in a vector
	return (i-1) * N_INTERNAL + j - 1;
}


/*double bordo(int i, int j) {
	return 1.0;
}*/


double g(double x,double y){  //function buondary conditions
	return 1.0;
}


double f(double x,double y){  //main function
	return 0.0;
}


/*double f(int i, int j) {
	return 0.0;

}*/


void setup_solution(Mesh &M, std::vector<double> &U, ComputeFunctionNode &b) {
	for (int i = 0; i < M.getDimension(); ++i) {
		for (int j = 0; j < M.getDimension(); ++j) {
			U[i * M.getDimension() + j] = b.getValue(i, j);
		}
	}
}




int main (int argc, char *argv[]) {
	Mesh mesh(-1.0, 0.0, 2.0, 2.0, N);
	ComputeFunctionNode bordo(&mesh,g);
	ComputeFunctionNode funzione(&mesh,f);


	std::vector<double> U_old(N*N);
	std::vector<double> U(N*N);

	setup_solution(mesh, U_old, bordo);
	setup_solution(mesh, U    , bordo);



	for (int k = 0; k < 1000; ++k) {
		// Jacobi iteration
		double h = mesh.getDiscretizationStep();
		for (int i = 1; i < N-1; ++i) {
			for (int j = 1; j < N-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				double b = funzione.getValue(i, j) * h * h;

				U[index(i,j)] =
					( b
					- U_old[index(i+1,j)]
					- U_old[index(i-1,j)]
					- U_old[index(i,j+1)]
					- U_old[index(i,j-1)] ) / (-4.0);
			}
		}


		std::swap(U, U_old);
	}

	for (auto x : U) {
		std::cout << x << std::endl;
	}


	//print_matrix(A);
	//print_vector(b);

	return 0;
}
