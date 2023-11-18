#include <iostream>
#include <cmath>
#include "ComputeFunctionNode.hpp"


#define N 10
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


/*int index(int i, int j) {
	return i * (N_INTERNAL) + j;
}*/


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




int main (int argc, char *argv[]) {
	double A[(N-2)*(N-2)][(N-2)*(N-2)] = {0.0};
	double b[(N-2)*(N-2)]              = {0.0};
	double h = 1.0 / static_cast<double>(N-1);


	Mesh mesh(0.0, 0.0, 1.0, 1.0, N);
	ComputeFunctionNode bordo(&mesh,g);
	ComputeFunctionNode funzione(&mesh,f);

	int row = 0;

	for (int i = 1; i < N-1; ++i) {
		for (int j = 1; j < N-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
			A[row][local_index(i,  j)] = -4.0;
			b[row]=funzione.getValue(i,j);

			if (i != 1) {
				A[row][local_index(i-1,j)] =  1.0;
			}
			else {
				b[row] -= bordo.getValue(i-1,j);
			}

			if (i != N-2) {
				A[row][local_index(i+1,j)] =  1.0;
			}
			else {
				b[row] -= bordo.getValue(i+1,j);
			}

			if (j != 1) {
				A[row][local_index(i,j-1)] =  1.0;
			}
			else {
				b[row] -= bordo.getValue(i,j-1);
			}

			if (j != N-2) {
				A[row][local_index(i,j+1)] =  1.0;
			}
			else {
				b[row] -= bordo.getValue(i,j+1);
			}

			++row;
		}
	}

	print_matrix(A);    
	//print_vector(b);

	return 0;
}
