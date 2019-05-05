#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <chrono>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>
#include <tree.hpp>


#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/System.hpp>


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
        max_length*tree_plus_tol/(tree_minus_tol*powf(2, max_depth))
    );
        
    std::cout << "Octree nodes: " << Node::node_count << " \n";


    // create image plane
    // [r g b  r g b ...]
    // create the rays
    // will be updated by void RayTrace::updateRays
    // [origin direction   origin direction...]
    FL_TYPE *rays = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h * 2 * 3);
    if(!rays)
    {
        throw std::runtime_error("can't allocate memory for rays");
    }

    FL_TYPE *nodes = &element_vector[0];
    
    fillTree(nodes, element_vector.size() / element_size);
    root.numElementsInside();
    std::cout<< "tree has: " << root.numElementsInside() << " elements\n";
  
    FL_TYPE lights[] = {lx, ly, lz};
    
    std::vector<FL_TYPE> vec_tree(0), vec_tree_ids(0);
    flattenTree(vec_tree, vec_tree_ids);

    sf::RenderWindow window(sf::VideoMode(w, h), "Render");    
    
    sf::Uint8 *sf_pixbuf = new sf::Uint8[w*h*4]; // rgba

    while(window.isOpen())
    {
        FL_TYPE *image_plane = \
        (FL_TYPE *)calloc(w * h * channels, sizeof(FL_TYPE));
        if(!image_plane)
        {
            throw std::runtime_error("can't allocate memory for image plane");
        }

        auto start_time = std::chrono::high_resolution_clock::now();
        sf::Texture texture;
        texture.create(w, h);
        sf::Sprite sprite(texture);

        sf::Event event;
        while(window.pollEvent(event))
        {
            if(event.type == sf::Event::Closed)
            {
                window.close();
            }
            if(event.type == sf::Event::KeyReleased)
            {
                if(event.key.code == sf::Keyboard::Up)
                {
                    camera[2] += 5;
                }
                else if(event.key.code == sf::Keyboard::Down)
                {
                    camera[2] -= 5;
                }
                if(event.key.code == sf::Keyboard::Left)
                {
                    camera[0] -= 0.5;
                }
                else if(event.key.code == sf::Keyboard::Right)
                {
                    camera[0] += 0.5;
                }

            }
            // sf::Image image()
        }

        RayTrace::updateRays(camera, rays);
        render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane, root);
        
        for(size_t _i = 0; _i < w*h; _i++)
        {
            FL_TYPE pp = image_plane[_i];
            pp = (pp > 1.0)?1.0:pp;
            pp = (pp < 0.0)?0.0:pp;
            sf_pixbuf[_i * 4] = sf::Uint8(pp * 255);
            sf_pixbuf[_i * 4 + 1] = sf::Uint8(pp * 255);
            sf_pixbuf[_i * 4 + 2] = sf::Uint8(pp * 255);
            sf_pixbuf[_i * 4 + 3] = (uint8_t) 255;
        }
        texture.update(sf_pixbuf);
        window.draw(sprite);
        window.display();

        auto end_time = std::chrono::high_resolution_clock::now();
        double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
        std::cout << "Render time taken: " << time_spent << "ms\n";
        free(image_plane);
    }

    // RayTrace::writeImage(image_plane, "a.ppm");
    free(rays);
}


// TODO: normalize rays