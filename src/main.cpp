#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>

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
    camera[2] = 25;

    lx = -lx;
    ly = -ly;
    lz = -lz;
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
        (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h * channels);
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
    // FL_TYPE nodes[] = {-3.5, -2, -2, 0.5, -2, -3, 0.5, 2, -3.1, 
                        // -1, -1.5, -2.5, 3, -1.5, -3.2, 3, 2.5, -2.4 };  
    
    // FL_TYPE nodes[] = {-2, -2, 2, 2, -2, 2, 2, 2, 2};

    FL_TYPE lights[] = {lx, ly, lz};
    render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane);
    // render(rays, nodes, 1, lights, 1, image_plane);
    
    // for(size_t _i = 0; _i < h; _i++)
    // {
    //     for(size_t _j = 0; _j < w; _j++)
    //     {
    //         std::cout << image_plane[_i * w + _j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    RayTrace::writeImage(image_plane, "a.ppm");
}


// TODO: normalize rays and normals