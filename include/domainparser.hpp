#ifndef DOMAINPARSER_HPP
#define DOMAINPARSER_HPP

#include <flags.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

std::vector<std::string> split(const std::string& s);
void DomainParser(std::string file_name,std::vector <double> &element_vector,unsigned int num_of_nodes,unsigned int num_of_elements);

#endif