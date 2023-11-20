#ifndef __UTILS_HPP__
#define __UTILS_HPP__


double norm(std::vector<double> &v){
	double result = 0.0;

	for (auto x : v) {
		result += x * x;
	}

	return std::sqrt(result);
}


#endif
