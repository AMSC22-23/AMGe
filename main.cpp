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


void setup_solution(vector<double> &U, ComputeFunctionNode &b) {  //initialize function U to border value
	for (int i = 0; i < b.getMesh()->getDimensionY(); ++i) {
		for (int j = 0; j < b.getMesh()->getDimensionX(); ++j) {
			U[i * b.getMesh()->getDimensionY() + j] = b.getValue(i, j);
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


void Gauss_Siedel_iteration(vector<double> &U, ComputeFunctionNode &function){  //1 cycle of Gauss Siedel
	double b;
	double hx_square= function.getMesh()->getDiscretizationStepX()*function.getMesh()->getDiscretizationStepX();
	double hy_square= function.getMesh()->getDiscretizationStepY()*function.getMesh()->getDiscretizationStepY();
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	

	for (int i = 1; i < function.getMesh()->getDimensionY()-1; ++i) {
			for (int j = 1; j < function.getMesh()->getDimensionX()-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function.getValue(i, j);
				U[index(i,j)] =
					( b
					+ ((U[index(i+1,j)]+U[index(i-1,j)])/hx_square)
					+ ((U[index(i,j+1)]+U[index(i,j-1)])/hy_square) ) / den;		
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

	setup_solution( U_old, bordo);
	setup_solution( U, bordo);

	//call to jacobi

	/*while (residual_norm > 1e-6) {
		// Jacobi iteration
		jacobi_iteration(U_old,U,funzione);
		residual_norm=compute_residual(U,residual,funzione);
		swap(U, U_old);
		count++;
	}
	cout<<residual_norm<<endl;
	cout<<count<<endl;*/



	//call to Gauss-Siedel

	while (residual_norm > 1e-6) {
		// Jacobi iteration
		Gauss_Siedel_iteration(U,funzione);
		residual_norm=compute_residual(U,residual,funzione);
		count++;
	}
	cout<<residual_norm<<endl;
	cout<<count<<endl;



	//export_to_matlab("U", U);


	return 0;
}
