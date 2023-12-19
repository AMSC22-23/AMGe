#ifndef __UTILS_HPP__
#define __UTILS_HPP__


#include <iostream>
#include <vector>
#include <string>


double norm(const std::vector<double> &v);
void export_to_matlab(const std::string &name, const std::vector<double> &v, std::ostream &out = std::cout);
bool is_2nplusone(int n);


#endif
