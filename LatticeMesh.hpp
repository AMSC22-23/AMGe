#ifndef __LATTICE_MESH_HPP__
#define __LATTICE_MESH_HPP__


#include "Point.hpp"
#include "FunctionNode.hpp"
#include "Mesh.hpp"
#include <iostream>
#include <string>
#include <tuple>


class LatticeMesh : Mesh {   //map from mesh to omega
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
		, hy(width / static_cast<double>(Ny-1))
	{
		for (int i = 1; i < Nx; ++i) {
			for (int j = 1; j < Ny; ++j) {
				inner_nodes.push_back(i + j * Nx);
			}
		}
	}


	Point getValue(int i, int j) {
		return Point(
			x + i * hx
			,y + j * hy
		);
	}


	int index(int i, int j) {
		return i * Ny + j;
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


	~LatticeMesh() {};


	const unsigned int Nx, Ny;
	const double x,y,width,height;
	const double hx,hy;
};


#endif
