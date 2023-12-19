#ifndef __POINT_HPP__
#define __POINT_HPP__

//@note: if everything is public just go with a `struct`
//@note: (minor) is you used an std::array could easily generalize to 3D
class Point {
public:
    Point(double _x, double _y) {
	x = _x;
	y = _y;
    }


    double x,y;
};


#endif
