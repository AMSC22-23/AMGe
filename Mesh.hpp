#include<iostream>

#define N 10
#define L 10
#define H 10
#include "FunctionNode.hpp"

class  Mesh: public FunctionNode<std::pair<double,double>> {   //map from mesh to omega


public:
    Mesh() : FunctionNode() {
        hx=static_cast<double>(L)/static_cast<double>(N);
        hy=static_cast<double>(H)/static_cast<double>(N);
    };

    std::pair<double, double> getValue(int i,int j) override{
		std::pair<double,double> a;
		a.first=hx*static_cast<double>(i);
		a.second=hy*static_cast<double>(j);
		return a;
	}

    ~Mesh() override = default;


    private:
        double hx;
        double hy;

};