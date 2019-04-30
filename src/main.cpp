#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <algorithm>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>
#include <tree.hpp>

#include <chrono>

#ifdef COLOR
    const size_t channels = 3; // rgb
#else
    const size_t channels = 1; // monochrome
#endif

std::vector<FL_TYPE> camera(3, 0);
std::vector<FL_TYPE> domain_limits(6, 0);

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
    DomainParser(file,element_vector, domain_limits);

    std::cout << "Parsing finished\n" << "Elements: " << element_vector.size() / element_size << "\n";
    std::cout << "Domain limits: [" << domain_limits[0] << ", "\
        << domain_limits[1] << ", "\
        << domain_limits[2] << "], ["\
        << domain_limits[3] << ", "\
        << domain_limits[4] << ", "\
        << domain_limits[5] << "]\n";
    
    FL_TYPE max_length = std::max(
        fabs(domain_limits[0] - domain_limits[3]),
        std::max(
            fabs(domain_limits[1] - domain_limits[4]),
            fabs(domain_limits[2] - domain_limits[5])
        )
    );
    
    std::cout << "Domain max length: " << max_length << "\n";

    Node::node_count = 0;
    Node root(
        domain_limits[0] - fabs(max_length*(1-tree_minus_tol)),
        domain_limits[4] + fabs(max_length*(1-tree_plus_tol)),
        domain_limits[2] - fabs(max_length*(1-tree_minus_tol)),
        max_length*tree_plus_tol/tree_minus_tol,
        2);
        
    std::cout << "Octree nodes: " << Node::node_count << " \n";


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
    
    fillTree(nodes, element_vector.size() / element_size);
    root.numElementsInside();
    std::cout<< "tree has: " << root.numElementsInside() << " elements\n";
  
    FL_TYPE lights[] = {lx, ly, lz};
    auto start_time = std::chrono::high_resolution_clock::now();
    render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane, root);
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";

    RayTrace::writeImage(image_plane, "a.ppm");
    free(image_plane);
    free(rays);
}


// TODO: normalize rays