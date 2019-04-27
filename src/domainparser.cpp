#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <domainparser.hpp>


std::vector<std::string> split(const std::string& s)
{
   std::vector<std::string> tokens;
   std::string token;
   tokens.reserve(9);
   std::istringstream tokenStream(s);
   while (getline(tokenStream, token,' '))
   {
      tokens.emplace_back(token);
   }
   return tokens;
}

void DomainParser(std::string file_name,std::vector <double> &element_vector,unsigned int num_of_nodes,unsigned int num_of_elements)
{
    unsigned int counter = 1;
    unsigned int counter1 = 0;
    unsigned int n_count;
    unsigned int e_count;
    std::vector <double> node_vector;
    int ia1,ia2,ia3;
    std::string line;
    std::ifstream infile; 
    infile.open(file_name);
    std::vector <std::string> data;
    while (getline(infile,line)) 
    {  
        if (counter == 5)
        {
            data = split(line);
            num_of_nodes = std::stoi(data[0]);
            node_vector.resize(3*num_of_nodes);
        }

        if ( counter > 5 && counter < num_of_nodes+6)
        {
            data = split(line);
            n_count = (counter - 6)*3;
            node_vector[n_count] = std::stod(data[1]);
            node_vector[n_count+1] = std::stod(data[2]);
            node_vector[n_count+2] = std::stod(data[3]);
        }

        if (counter == num_of_nodes+8)
        {
            data = split(line);
            num_of_elements = std::stoi(data[0]);
            element_vector.resize(9*num_of_elements);
            
        }

        if ( counter > num_of_nodes + 8 && counter < num_of_nodes+9+num_of_elements)
        {
            data = split(line);
            if (std::stoi(data[1]) == 2)
            {
                ia1 = (std::stoi(data[5])-1)*3;
                ia2 = (std::stoi(data[6])-1)*3;
                ia3 = (std::stoi(data[7])-1)*3;
                e_count = counter1*9;
                element_vector[e_count] = node_vector[ia1];
                element_vector[e_count+1] = node_vector[ia1+1];
                element_vector[e_count+2] = node_vector[ia1+2];
                element_vector[e_count+3] = node_vector[ia2];
                element_vector[e_count+4] = node_vector[ia2+1];
                element_vector[e_count+5] = node_vector[ia2+2];
                element_vector[e_count+6] = node_vector[ia3];
                element_vector[e_count+7] = node_vector[ia3+1];
                element_vector[e_count+8] = node_vector[ia3+2];
                counter1++;
            }    
        }
        counter++;
    }
    num_of_elements = counter1;
    element_vector.resize(num_of_elements*9);
    infile.close();
}