#ifndef DOMAINPARSER_HPP
#define DOMAINPARSER_HPP

#include <flags.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


std::vector<std::string> split(const std::string& s);
void DomainParser(
    std::string file_name, std::vector <FL_TYPE> &element_vector,
    std::vector<FL_TYPE> &domain
    );

#endif