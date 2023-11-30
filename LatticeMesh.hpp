#ifndef __LATTICE_MESH_HPP__
#define __LATTICE_MESH_HPP__


#include "Point.hpp"
#include "Mesh.hpp"
#include <iostream>
#include <string>
#include <tuple>
#include "Utils.hpp"


class LatticeMesh : Mesh {
public:
	LatticeMesh(
		  double _x
		, double _y
		, double _width
		, double _height
		, unsigned int _Nx
		, unsigned int _Ny
		) :
		  x(_x)
		, y(_y)
		, width(_width)
		, height(_height)
		, Nx(_Nx)
		, Ny(_Ny)
		, hx(width / static_cast<double>(Nx-1))
		, hy(width / static_cast<double>(Ny-1)) {

		for (int j = 1; j < Ny-1; ++j) {
			for (int i = 1; i < Nx-1; ++i) {
				inner_nodes.push_back(i + j * Nx);
			}
		}


		// se una LatticeMesh fosse 3x65, potrebbe ammettere ancora suddivisioni sulla seconda dimensione, sistemare
		is_minimal = 
			   (Ny <= 3)
			or (Nx <= 3)
			or (!is_2nplusone(Nx))
			or (!is_2nplusone(Ny));
	}


	Index index(int i, int j) {
		return i + j * Nx;
	}


	std::pair<int,int> inverse_index(Index i){
		return std::pair<int,int>(i%Nx, i/Nx);
	} 


	void save(std::ostream &);


	virtual Point get_point_coordinate(Index i) override {
		return Point(
			  x + (i % Nx) * hx
			, y + (i / Nx) * hy
		);
	}


	virtual void evaluate_function(std::vector<double> &U, double (*f)(double, double)) override {
		for (Index i = 0; i < Nx*Ny; ++i) {
			Point point = get_point_coordinate(i);
			U[i] = f(point.x, point.y);
		}
	}


	virtual std::tuple<Index, Index, Index, Index> get_cardinal_neighbours(Index i) override {
		// nessun controllo fatto ai bordi, usare con cautela
		return std::make_tuple(
			/* NORD  */   i + Nx
			/* SUD   */ , i - Nx
			/* OVEST */ , i + 1
			/* EST   */ , i - 1
		);
	}


	virtual std::vector<Index> get_inner_nodes() override {
		return inner_nodes;
	}


	LatticeMesh build_coarse() {
		return LatticeMesh(
			  x
			, y
			, width
			, height
			, Nx / 2 + 1
			, Ny / 2 + 1
		);
	}


	void project_on_coarse(const std::vector<double> &u_fine, std::vector<double> &u_coarse) {
		std::pair<int,int> p;
		int j = 0;

		for (Index i = 0;i < Nx*Ny; i++) {
			auto [k, l] = inverse_index(i);

			if(k % 2 ==0 and l % 2 == 0) {
				u_coarse[j] = u_fine[i];
				j++;
			}
		}
	}


	void interpolate_on_fine(LatticeMesh &coarse, std::vector<double> &u_fine, const std::vector<double> &u_coarse) {
		// devo iterare solo sui nodi interni
		for (int i = 1; i < Nx-1; ++i) {
			for (int j = 1; j < Ny-1; ++j) {

				if (i % 2 == 0 and j % 2 == 0) {
					u_fine[index(i,j)] = u_coarse[coarse.index(i/2, j/2)];
				}
				else if (i % 2 == 1 and j % 2 == 0) {
					u_fine[index(i,j)] =
						  0.5 * u_coarse[coarse.index(    i/2, j/2)]
						+ 0.5 * u_coarse[coarse.index(1 + i/2, j/2)];
				}
				else if(i % 2 == 0 and j % 2 == 1) {
					u_fine[index(i,j)] =
						  0.5 * u_coarse[coarse.index(i/2,     j/2)]
						+ 0.5 * u_coarse[coarse.index(i/2, 1 + j/2)];
				}
				else {
					u_fine[index(i,j)] =
						  0.25 * u_coarse[coarse.index(    i/2,     j/2)]
						+ 0.25 * u_coarse[coarse.index(    i/2, 1 + j/2)]
						+ 0.25 * u_coarse[coarse.index(1 + i/2,     j/2)]
						+ 0.25 * u_coarse[coarse.index(1 + i/2, 1 + j/2)];
				}
			}
		}
	}


	~LatticeMesh() {};


	const unsigned int Nx, Ny;
	const double x,y,width,height;
	const double hx,hy;
	bool is_minimal;
};


#endif
