#include "Point.hpp"
#include "FunctionNode.hpp"


class Mesh : public FunctionNode<Point> {   //map from mesh to omega
public:
    Mesh(double _x, double _y, double _width, double _height, unsigned int _Nx, unsigned int _Ny) : FunctionNode() {
	x      = _x;
	y      = _y;
	width  = _width;
	height = _height;
	Nx     = _Nx;
	Ny     = _Ny;
    };


    Point getValue(int i, int j) override {
	return Point(
	     x + i * getDiscretizationStepX()
	    ,y + j * getDiscretizationStepY()
	);
    }


    unsigned int getDimensionX() {
	return Nx;
    }

    unsigned int getDimensionY() {
	return Ny;
    }


    double getDiscretizationStepX() {
	return width / static_cast<double>(Nx);
    }

    double getDiscretizationStepY() {
	return height / static_cast<double>(Ny);
    }

    double getWidth(){
        return width;
    }
    
    double getHeight(){
        return height;
    }

    double getXCenter(){
        return x;
    }

    double getYCenter(){
        return y;
    }


    ~Mesh() override = default;


private:
    unsigned int Nx, Ny;
    double x,y,width,height;
};
