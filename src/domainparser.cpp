#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <domainparser.hpp>
#include <flags.hpp>


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
    FL_TYPE ax, ay, az, bx, by, bz, cx, cy, cz; // node coords
    FL_TYPE e01x, e01y, e01z, e02x, e02y, e02z, nx, ny, nz; // edges and surface normals
    FL_TYPE normal_length;
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
            element_vector.resize(element_size*num_of_elements);
            
        }

        if ( counter > num_of_nodes + 8 && counter < num_of_nodes+element_size+num_of_elements)
        {
            data = split(line);
            if (std::stoi(data[1]) == 2)
            {
                ia1 = (std::stoi(data[5])-1)*3;
                ia2 = (std::stoi(data[6])-1)*3;
                ia3 = (std::stoi(data[7])-1)*3;
                e_count = counter1*element_size;
                element_vector[e_count] = ax = node_vector[ia1];
                element_vector[e_count+1] = ay = node_vector[ia1+1];
                element_vector[e_count+2] = az = node_vector[ia1+2];
                element_vector[e_count+3] = bx = node_vector[ia2];
                element_vector[e_count+4] = by = node_vector[ia2+1];
                element_vector[e_count+5] = bz = node_vector[ia2+2];
                element_vector[e_count+6] = cx = node_vector[ia3];
                element_vector[e_count+7] = cy = node_vector[ia3+1];
                element_vector[e_count+8] = cz = node_vector[ia3+2];
                
                // get the edges
                e01x = bx - ax;
                e01y = by - ay;
                e01z = bz - az;
                
                e02x = cx - ax;
                e02y = cy - ay;
                e02z = cz - az;
                
                // calculate the surface normal
                nx = e01y * e02z - e01z * e02y;
                ny = e01z * e02x - e01x * e02z;
                nz = e01x * e02y - e01y * e02x;

                normal_length = sqrt(nx * nx + ny * ny + nz * nz);

                element_vector[e_count+9] = nx/normal_length;
                element_vector[e_count+10] = ny/normal_length;
                element_vector[e_count+11] = nz/normal_length;

                element_vector[e_count+12] = e01x;
                element_vector[e_count+13] = e01y;
                element_vector[e_count+14] = e01z;
                element_vector[e_count+15] = e02x;
                element_vector[e_count+16] = e02y;
                element_vector[e_count+17] = e02z;

                element_vector[e_count+18] = normal_length;

                counter1++;
            }    
        }
        counter++;
    }
    num_of_elements = counter1;
    element_vector.resize(num_of_elements*element_size);
    infile.close();
}