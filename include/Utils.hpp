#ifndef __UTILS_HPP__
#define __UTILS_HPP__


#include <bitset>
#include <cmath>
#include <string>


double norm(const std::vector<double> &v){
	double result = 0.0;

	for (double x : v) {
		result += x * x;
	}

	return std::sqrt(result);
}


void export_to_matlab(const std::string &name, const std::vector<double> &v, std::ostream &out = std::cout){
	out
		<< "# name: "      << name     << std::endl
		<< "# type: matrix"            << std::endl
		<< "# rows: 1"                 << std::endl
		<< "# columns: "   << v.size() << std::endl;


	for (auto x : v) {
		out << x << ' ';
	}


	out << '\n' << '\n';
}


bool is_2nplusone(int n){
	int count=0;
	std::bitset<64> bits(n);
	
	if(bits.test(0)==0){
		return false;
	}

	for (int i=1;i<64;i++){
		
		if(bits.test(i)==1){
			count++;
			
			if(count==2){
				return false;
			}

		}
	
	}

	if(count==0){
		return false;
	}
	
	return true;
}


#endif
