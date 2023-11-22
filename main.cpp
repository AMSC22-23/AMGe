#include <iostream>
#include <cmath>
#include <vector>
#include "ComputeFunctionNode.hpp"
#include "Utils.hpp"


using namespace std;


#define Mx 100
#define My 500


double g(double x,double y){  //function boundary conditions
	return 1.0;
}


double f(double x,double y){  //main function
	return 1.0;
}


void setup_solution(vector<double> &U, ComputeFunctionNode &b,Mesh m) {  //initialize function U to border value
	for (int i = 0; i < m.Nx; ++i) {
		for (int j = 0; j < m.Ny; ++j) {
			U[i * m.Ny + j] = b.getValue(i, j);
		}
	}
}


void jacobi_iteration(vector<double> &U_old, vector<double> &U, ComputeFunctionNode &function,Mesh m){  //1 cycle of Jacobi
	double b;
	double hx_square= m.hx*m.hx;
	double hy_square= m.hy*m.hy;
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	

	for (int i = 1; i < m.Nx-1; ++i) {
			for (int j = 1; j < m.Ny-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function.getValue(i, j);
				U[m.index(i,j)] =
					( b
					+ ((U_old[m.index(i+1,j)]+U_old[m.index(i-1,j)])/hx_square)
					+ ((U_old[m.index(i,j+1)]+U_old[m.index(i,j-1)])/hy_square) ) / den;
			}
	}
}


void Gauss_Siedel_iteration(vector<double> &U, ComputeFunctionNode &function,Mesh m){  //1 cycle of Gauss Siedel
	double b;
	double hx_square= m.hx*m.hx;
	double hy_square= m.hy*m.hy;
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	

	for (int i = 1; i < m.Nx-1; ++i) {
			for (int j = 1; j < m.Ny-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function.getValue(i, j);
				U[m.index(i,j)] =
					( b
					+ ((U[m.index(i+1,j)]+U[m.index(i-1,j)])/hx_square)
					+ ((U[m.index(i,j+1)]+U[m.index(i,j-1)])/hy_square) ) / den;
			}
	}
}


double compute_residual(vector<double> &U, vector<double> &residual, ComputeFunctionNode &function,Mesh m){    //compute residual vecotr
	double h_squarex= m.hx*m.hx;
	double h_squarey= m.hy*m.hy;
	double factor =-2.0*(h_squarex+h_squarey);
	
	for (int i = 1; i < m.Nx-1; ++i) {
		for (int j = 1; j < m.Ny-1; ++j) {
			residual[m.index(i,j)] =
				 (h_squarex*h_squarey*function.getValue(i,j))
				+((U[m.index(i,j-1)]+U[m.index(i,j+1)])*h_squarex)
				+((U[m.index(i-1,j)]+U[m.index(i+1,j)])*h_squarey)
				+factor*U[m.index(i,j)];
		}
	}

	return norm(residual);
}


int main (int argc, char *argv[]) {
	Mesh mesh(-1.0, 0.0, 2.0, 2.0, Mx,My);

	ComputeFunctionNode bordo(&mesh,g);
	ComputeFunctionNode funzione(&mesh,f);


	vector<double> U_old(mesh.Nx*mesh.Ny);
	vector<double> U(mesh.Nx*mesh.Ny);
	
	vector<double> residual(mesh.Nx*mesh.Ny);
	double residual_norm=100;

	int count =0;

	setup_solution( U_old, bordo,mesh);
	setup_solution( U, bordo,mesh);

	//call to jacobi

	/*while (residual_norm > 1e-6) {
		// Jacobi iteration
		jacobi_iteration(U_old,U,funzione,mesh);
		residual_norm=compute_residual(U,residual,funzione);
		swap(U, U_old);
		count++;
	}
	cout<<residual_norm<<endl;
	cout<<count<<endl;*/



	//call to Gauss-Siedel

	while (residual_norm > 1e-6) {
		// Jacobi iteration
		Gauss_Siedel_iteration(U,funzione,mesh);
		residual_norm=compute_residual(U,residual,funzione,mesh);
		count++;
	}
	cout<<residual_norm<<endl;
	cout<<count<<endl;



	//export_to_matlab("U", U);


	return 0;
}
