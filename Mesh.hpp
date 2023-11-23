#ifndef __MESH_HPP__
#define __MESH_HPP__


#include <vector>
#include "Point.hpp"


using Index = int;


class Mesh {
public:
	Mesh(){};
	virtual Point get_point_coordinate(Index i) = 0;
	virtual void evaluate_function(std::vector<double> &U, double (*f)(double, double)) = 0;
	virtual std::tuple<Index, Index, Index, Index> get_cardinal_neighbours(Index i) = 0;
	~Mesh(){};


	// mi piacerebbe aggiungere un iteratore ai nodi interni, ma poi non sarebbe più così facile da parallelizzare
private:
	std::vector<Index> inner_nodes;
};


#endif
