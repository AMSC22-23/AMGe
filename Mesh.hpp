#ifndef __MESH_HPP__
#define __MESH_HPP__


#include <vector>
#include "Point.hpp"


class Mesh {
public:
	Mesh(){};
	virtual Point get_point_coordinate(int index) = 0;
	virtual void evaluate_function(std::vector<double> &U, double (*f)(double, double)) = 0;
	~Mesh(){};
};


#endif
