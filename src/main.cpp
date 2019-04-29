#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>

#include <chrono>

#ifdef COLOR
    const size_t channels = 3; // rgb
#else
    const size_t channels = 1; // monochrome
#endif

std::vector<FL_TYPE> camera(3, 0);


int main(int argc, char** argv)
{
    camera[0] = 30;
    camera[1] = 30;
    camera[2] = 60;

    std::vector <FL_TYPE> element_vector;
    unsigned int num_of_nodes;
    unsigned int num_of_elements;
    std::string file = "shadow";
    std::cout << "Parsing GMSH domain \"" << file << "\"\n";
    DomainParser(file,element_vector,num_of_nodes,num_of_elements);
    std::cout << "Parsing finished\n" << "Elements: " << element_vector.size() / element_size << "\n";
    
    // create image plane
    // [r g b  r g b ...]
    FL_TYPE *image_plane = \
        (FL_TYPE *)calloc(w * h * channels, sizeof(FL_TYPE));
    if(!image_plane)
    {
        throw std::runtime_error("can't allocate memory for image plane");
    }
    // create the rays
    // will be updated by void RayTrace::updateRays
    // [origin direction   origin direction...]
    FL_TYPE *rays = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h * 2 * 3);
    if(!rays)
    {
        throw std::runtime_error("can't allocate memory for rays");
    }
    RayTrace::updateRays(camera, rays);

    FL_TYPE *nodes = &element_vector[0];
  
    FL_TYPE lights[] = {lx, ly, lz};
    auto start_time = std::chrono::high_resolution_clock::now();
    render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane);
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";

    RayTrace::writeImage(image_plane, "a.ppm");
}


// TODO: normalize rays