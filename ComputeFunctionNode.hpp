#ifndef __COMPUTE_FUNCTION_NODE_HPP__
#define __COMPUTE_FUNCTION_NODE_HPP__


#include "LatticeMesh.hpp"
#include "Point.hpp"


class ComputeFunctionNode : public FunctionNode<double> {  //compute mapping of a function from mesh index to R
public:
    ComputeFunctionNode(LatticeMesh *_M, double (*_f)(double, double)) : FunctionNode() {
	M = _M;
	f = _f;
    };


    double getValue(int i, int j) override {
	Point p = M->getValue(i, j);
	return f(p.x, p.y);
    }



    virtual ~ComputeFunctionNode() override = default;


private:
    LatticeMesh *M;
    double (*f)(double, double);
};


#endif
