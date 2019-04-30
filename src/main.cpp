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
    DomainParser(file,element_vector,num_of_nodes,num_of_elements, domain_limits);
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
    Node root(domain_limits[0], domain_limits[4], domain_limits[2], max_length, 2);
    // Node root(10, 10, 10, 10, 4);
    // size_t *nodes_hit = (size_t *)calloc(Node::node_count, sizeof(size_t));
    // int idx = 0;
    std::cout << "Octree nodes: " << Node::node_count << "\n";
    // root.rayIntersection(16, 6, 25, 0, 0, -1, nodes_hit, idx);
    // for(size_t _k = 0; _k < Node::node_count; _k++)
    // {
    //     std::cout << nodes_hit[_k] << "\n";
    //     // std::cout << Node::all_nodes[_k]->numElementsInside() << "\n";
    // }
    

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
    render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane, root);
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";

    RayTrace::writeImage(image_plane, "a.ppm");
}


// TODO: normalize rays