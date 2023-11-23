#include <iostream>
#include <cmath>
#include <vector>
#include "LatticeMesh.hpp"
#include "ComputeFunctionNode.hpp"
#include "Utils.hpp"


using namespace std;


double g(double x,double y){  //function boundary conditions
	return 1.0;
}


double f(double x,double y){  //main function
	return 1.0;
}


void setup_solution(vector<double> &U, ComputeFunctionNode &b,LatticeMesh m) {  //initialize function U to border value
	for (int i = 0; i < m.Nx; ++i) {
		for (int j = 0; j < m.Ny; ++j) {
			U[i * m.Ny + j] = b.getValue(i, j);
		}
	}
}


void jacobi_iteration(vector<double> &U_old, vector<double> &U, LatticeMesh m,unsigned int level, vector<double> &function){  //1 cycle of Jacobi
	double b;
	double hx_square= m.hx*pow(2,level)*m.hx*pow(2,level);
	double hy_square= m.hy*pow(2,level)*m.hy*pow(2,level);
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	
	if(m.Nx/pow(2,level)>3&&m.Ny>3){

	for (int i = 1; i < (m.Nx/pow(2,level))-1; ++i) {
			for (int j = 1; j < m.Ny-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function[m.index(i,j)];
				U[m.index(i,j)] =
					( b
					+ ((U_old[m.index(i+1,j)]+U_old[m.index(i-1,j)])/hx_square)
					+ ((U_old[m.index(i,j+1)]+U_old[m.index(i,j-1)])/hy_square) ) / den;
			}
	}
	}
}


/*void Jacobi2(vector<double> &U_old, vector<double> &U, vector<double> &function,LatticeMesh m){
	double b;
	double hx_square= m.hx*m.hx*4.0;  //m.hx*m.hx*pow(2,level)
	double hy_square= m.hy*m.hy*4.0;
	double den      = ((2.0)*((1/hx_square)+(1/hy_square))); 
	for(int i=2;i<m.Nx-1;i+=2){  //i=pow(2,level), i+=pow(2,level)
		for(int j=2;j<m.Nx-1;j+=2){
				b = function[m.index(i,j)];
				U[m.index(i,j)] =
					( b
					+ ((U_old[m.index(i+1,j)]+U_old[m.index(i-1,j)])/hx_square)
					+ ((U_old[m.index(i,j+1)]+U_old[m.index(i,j-1)])/hy_square) ) / den;
		}
	}
}*/



void Gauss_Siedel_iteration(vector<double> &U, vector<double> &function,LatticeMesh m){  //1 cycle of Gauss Siedel
	double b;
	double hx_square= m.hx*m.hx;
	double hy_square= m.hy*m.hy;
	double den      = ((2.0)*((1/hx_square)+(1/hy_square)));
	

	for (int i = 1; i < m.Nx-1; ++i) {
			for (int j = 1; j < m.Ny-1; ++j) {     //iterating over nodes and computing for each one a row of matrix A
				b = function[m.index(i,j)];
				U[m.index(i,j)] =
					( b
					+ ((U[m.index(i+1,j)]+U[m.index(i-1,j)])/hx_square)
					+ ((U[m.index(i,j+1)]+U[m.index(i,j-1)])/hy_square) ) / den;
			}
	}
}


double compute_residual(vector<double> &U, vector<double> &residual, ComputeFunctionNode &function,LatticeMesh m){    //compute residual vecotr
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


int main (int argc, char *argv[]) {
	const int Nx = 50;
	const int Ny = 50;


	LatticeMesh mesh(-1.0, 0.0, 2.0, 2.0, Nx,Ny);
	vector<double> Old(Nx*Ny);
	vector<double> U(Nx*Ny);
	vector<double> F(Nx*Ny);


	mesh.evaluate_function(U, g);
	mesh.evaluate_function(Old, g);
	mesh.evaluate_function(F, f);


	for (int it = 0; it < 1000; ++it) {
		jacobi(mesh, U, Old, F);
		std::swap(U, Old);
	}


	export_to_matlab("U", U);


	return 0;
}
