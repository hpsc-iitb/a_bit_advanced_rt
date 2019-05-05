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

__device__
void d_normalize(
    float &x, float &y, float &z,
    float &nx, float &ny, float &nz
)
{
    float l = sqrt(x * x + y * y + z * z);
    nx = x / l;
    ny = y / l;
    nz = z / l;
}

__device__
bool d_rayTreeIntersection(
    float rox, float roy, float roz,
    float rdx, float rdy, float rdz,
    float *tree, float *tree_idx, unsigned int *ids,
    unsigned int &idx, unsigned int current_node,
    bool normalized
)
{
    unsigned int ctpos = (unsigned int)tree_idx[current_node];
    if(tree[ctpos + 26] == 0)
    {
        return false;
    }
    rdx = (!rdx)?1e-8:rdx;
    rdy = (!rdy)?1e-8:rdy;
    rdz = (!rdz)?1e-8:rdz;
    
    if(!normalized)
    {
        d_normalize(rdx, rdy, rdz, rdx, rdy, rdz);
    }

    float swap_tmp;

    float txmin = (tree[ctpos+1+0] - rox) / rdx; // vertex a
    float txmax = (tree[ctpos+1+3] - rox) / rdx; // vertex b
    
    float tymin = (tree[ctpos+1+7] - roy) / rdy; // vertex c
    float tymax = (tree[ctpos+1+1] - roy) / rdy; // vertex a

    float tzmin = (tree[ctpos+1+2] - roz) / rdz; // vertex a
    float tzmax = (tree[ctpos+1+14] - roz) / rdz; // vertex e


    // account for negatives
    if(rdx < 0)
    {
        swap_tmp = txmax;
        txmax = txmin;
        txmin = swap_tmp;
    }
    if(rdy < 0)
    {
        swap_tmp = tymax;
        tymax = tymin;
        tymin = swap_tmp;
    }
    if(rdz < 0)
    {
        swap_tmp = tzmax;
        tzmax = tzmin;
        tzmin = swap_tmp;
    }

    float tmin = fmax(
        txmin, fmax(tymin, tzmin)
    );

    float tmax = fmin(
        txmax, fmin(tymax, tzmax)
    );

    if(tmin <= tmax)
    {
        // there exists a parameter t for which ray intersects nodes
        if(tree[ctpos+25])
        {
            // last node, add own id to intersecting nodes
            ids[idx++] = tree[ctpos];
            return true;
        }
        else
        {
            bool retval = false;
            for(unsigned int _k = 0; _k < 8; _k++)
            {
                retval = retval | d_rayTreeIntersection(
                    rox, roy, roz, rdx, rdy, rdz,
                    tree, tree_idx, ids, idx, tree[ctpos+27+_k],
                    true
                );
            }
            return retval;
        }
    }
    else
    {
        return false;
    }
}

__global__
void render_gpu(
    float *rays, float *elems, float *tree, float *tree_idx,
    float *image_plane, float *lights, int num_lights,
    int num_tree_sz, int elem_sz
)
{
    float rox, roy, roz, rdx, rdy, rdz, nrdx, nrdy, nrdz;
    unsigned int pixel_num = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int ray_start = pixel_num * 6;
    image_plane[pixel_num] = 0.0;
    unsigned int ids[20];
    unsigned int idx = 0;

    // unsigned int ct = 0;
    
    unsigned long _j;
    
    float t; // parametric eqn
    // float illum;
    float u, v; // barycentric coords
    float ax, ay, az;//, bx, by, bz, cx, cy, cz; // nodes
    // float px, py, pz; // intersection point
    float nx, ny, nz, nl; // surface normal
    float e01x, e01y, e01z, e02x, e02y, e02z;
    // float ray_length_1; // 1/shadow_ray_length
    float tx, ty, tz, D, Dy, Dz;

    unsigned int hit_elem = -1;

    rox = rays[ray_start];
    roy = rays[ray_start+1];
    roz = rays[ray_start+2];

    rdx = rays[ray_start+3];
    rdy = rays[ray_start+4];
    rdz = rays[ray_start+5];
    // image_plane[pixel_num] = 0.5;
    float min_dis = 1e20;
    unsigned int leaf_pos;
    if(d_rayTreeIntersection(
        rox, roy, roz, rdx, rdy, rdz,
        tree, tree_idx, ids, idx, 0, false
    ))
    {
        for(unsigned int _i = 0; _i < idx; _i++)
        {
            leaf_pos = (unsigned int)tree_idx[ids[_i]];
            for(unsigned int _m = 0; _m < tree[leaf_pos+26]; _m++)
            {
                _j = tree[leaf_pos + 27 + _m] * 19;
                ax = elems[_j];
                ay = elems[_j+1];
                az = elems[_j+2];

                nx = elems[_j + 9];
                ny = elems[_j + 10];
                nz = elems[_j + 11];

                e01x = elems[_j + 12];
                e01y = elems[_j + 13];
                e01z = elems[_j + 14];
                
                e02x = elems[_j + 15];
                e02y = elems[_j + 16];
                e02z = elems[_j + 17];

                nl = elems[_j + 18];

                D = - nl*(nx*rdx + ny*rdy + nz*rdz);  // |-d e1 e2| = -n.d

                if(fabs(D) < 1e-6)
                {
                    continue;
                }

                tx = rox - ax;
                ty = roy - ay;
                tz = roz - az;
                // FL_TYPE px, py, pz;
                // cross(rdx, rdy, rdz, e02x, e02y, e02z, px, py, pz);

                Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);
                u = Dy / D;

                if(u < 0 || u > 1)
                {
                    continue;
                }

                Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);

                v = Dz / D;

                if(v < 0 || v + u > 1)
                {
                    continue;
                }

                t = (tx*nx + ty*ny + tz*nz)*nl/D;

                if(t < 0)
                {
                    continue;
                }
                if(t < min_dis)
                {
                    // near_ele_num = _j+_m;
                    hit_elem = _j;
                    min_dis = t;
                }
            }
        }
    }

    if(min_dis < 1e19)
    {
        rox = rox + rdx*min_dis;
        roy = roy + rdy*min_dis;
        roz = roz + rdz*min_dis;

        rdx = lights[0] - rox;
        rdy = lights[1] - roy;
        rdz = lights[2] - roz;

        d_normalize(rdx, rdy, rdz, nrdx, nrdy, nrdz);

        image_plane[pixel_num] = fabs(elems[hit_elem+9]*nrdx + elems[hit_elem+10]*nrdy + elems[hit_elem+11]*nrdz);
    }
    
    idx = 0;
    // return;

    if(d_rayTreeIntersection(
        rox, roy, roz, rdx, rdy, rdz,
        tree, tree_idx, ids, idx, 0, false
    ))
    {
        // return;
    // image_plane[pixel_num] = ((float)idx)/10.0;
        for(unsigned int _i = 0; _i < idx; _i++)
        {
            leaf_pos = (unsigned int)tree_idx[ids[_i]];
            for(unsigned int _m = 0; _m < tree[leaf_pos+26]; _m++)
            {
                _j = tree[leaf_pos + 27 + _m] * 19;
                if(_j == hit_elem)
                {
                    continue;
                }
                ax = elems[_j];
                ay = elems[_j+1];
                az = elems[_j+2];

                // bx = elems[_j+3];
                // by = elems[_j+4];
                // bz = elems[_j+5];

                // cx = elems[_j+6];
                // cy = elems[_j+7];
                // cz = elems[_j+8];

                nx = elems[_j + 9];
                ny = elems[_j + 10];
                nz = elems[_j + 11];

                e01x = elems[_j + 12];
                e01y = elems[_j + 13];
                e01z = elems[_j + 14];
                
                e02x = elems[_j + 15];
                e02y = elems[_j + 16];
                e02z = elems[_j + 17];

                nl = elems[_j + 18];

                D = - nl*(nx * rdx + ny*rdy + nz*rdz);  // |-d e1 e2| = -n.d

                if(fabs(D) < 1e-6)
                {
                    continue;
                }

                tx = rox - ax;
                ty = roy - ay;
                tz = roz - az;
                // FL_TYPE px, py, pz;
                // cross(rdx, rdy, rdz, e02x, e02y, e02z, px, py, pz);

                Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);
                u = Dy / D;

                if(u < 0 || u > 1)
                {
                    continue;
                }

                Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);

                v = Dz / D;

                if(v < 0 || v + u > 1)
                {
                    continue;
                }

                t = (tx*nx + ty*ny + tz*nz)*nl/D;

                if(t < 0)
                {
                    continue;
                }
                image_plane[pixel_num] = 0.0;
                return;
            }
        }
    }

}

int main(int argc, char** argv)
{
    camera[0] = 30;
    camera[1] = 30;
    camera[2] = 60;

    std::vector <FL_TYPE> element_vector;
    // unsigned int num_of_nodes;
    // unsigned int num_of_elements;
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

    FL_TYPE *image_plane = \
    (FL_TYPE *)calloc(w * h * channels, sizeof(FL_TYPE));
    if(!image_plane)
    {
        throw std::runtime_error("can't allocate memory for image plane");
    }


    float *d_tree, *d_tree_idx, *d_elems, *d_rays, *d_image_plane, *d_lights;
    
    cudaMalloc((void **)&d_rays, w*h*6*sizeof(FL_TYPE));
    cudaMalloc((void **)&d_image_plane, w*h*sizeof(FL_TYPE));
    cudaMalloc((void **)&d_lights, 3*sizeof(FL_TYPE));
    cudaMalloc((void **)&d_elems, element_vector.size()*sizeof(FL_TYPE));
    cudaMalloc((void **)&d_tree, vec_tree.size()*sizeof(FL_TYPE));
    cudaMalloc((void **)&d_tree_idx, vec_tree_ids.size()*sizeof(FL_TYPE));
    
    cudaMemcpy(d_image_plane, image_plane, w*h*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(d_lights, lights, 3*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(d_elems, nodes, element_vector.size()*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tree, &vec_tree[0], vec_tree.size()*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tree_idx, &vec_tree_ids[0], vec_tree_ids.size()*sizeof(FL_TYPE), cudaMemcpyHostToDevice);

    
    while(window.isOpen())
    {
        // image_plane = \
        // (FL_TYPE *)calloc(w * h * channels, sizeof(FL_TYPE));
        // if(!image_plane)
        // {
        //     throw std::runtime_error("can't allocate memory for image plane");
        // }

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

                else if(event.key.code == sf::Keyboard::Key::S)
                {
                    RayTrace::writeImage(image_plane, "a.ppm");
                }

            }
            // sf::Image image()
        }

        RayTrace::updateRays(camera, rays);

        // render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane, root);

        cudaMemcpy(d_rays, rays, w*h*6*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
        render_gpu<<<w, h>>>(
            d_rays, d_elems, d_tree, d_tree_idx, d_image_plane,
            d_lights, 1, vec_tree_ids.size(), element_vector.size()/element_size
        );
   
        cudaMemcpy(image_plane, d_image_plane, w*h*sizeof(FL_TYPE), cudaMemcpyDeviceToHost);

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
    }
    free(image_plane);
    free(rays);
    delete sf_pixbuf;
    cudaFree(d_rays);
    cudaFree(d_image_plane);
    cudaFree(d_lights);
    cudaFree(d_elems);
    cudaFree(d_tree);
    cudaFree(d_tree_idx);
}


// TODO: normalize rays