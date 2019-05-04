#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>
#include <cuda.h>
#include <cuda_runtime_api.h>

#include <chrono>

#ifdef COLOR
    const size_t channels = 3; // rgb
#else
    const size_t channels = 1; // monochrome
#endif

std::vector<FL_TYPE> camera(3, 0);



__global__ 
void checkintersection(float* nodes,float* rays, float* ts, size_t elem_size)
{
    int _i = blockIdx.x * blockDim.x + threadIdx.x;
    // int _j = get_global_id(1);
    float ax,ay,az,bx,by,bz,cx,cy,cz,e01x,e01y,e01z,e02x,e02y,e02z,\
    tx,ty,tz,u,v,D, Dy, Dz, t, nx, ny, nz, rox, roy, roz, rdx, rdy, rdz;

    float min_dis = 1e20;

    unsigned int sz, eidx;
    
    __shared__ float elems_shared[100*9];
    

    rox = rays[_i*6];
    roy = rays[_i*6+1];
    roz = rays[_i*6+2];
        
    rdx = rays[_i*6+3];
    rdy = rays[_i*6+4];
    rdz = rays[_i*6+5];
    
    for(unsigned int _j = 0; _j < elem_size*9; _j+=100*9)
    {
        sz = ((elem_size*9 - _j)>=100*9)?100*9:(elem_size*9 - _j);
        if(threadIdx.x == 0)
        {
            for(unsigned int _k = 0; _k < sz; _k++)
            {
                elems_shared[_k] = nodes[_j + _k];
            }
        }
        __syncthreads();
        for(unsigned int _m = 0; _m < sz; _m+= 9)
        {
            eidx = _j + _m;
            ax = elems_shared[_m];
            ay = elems_shared[_m+1];
            az = elems_shared[_m+2];

            bx = elems_shared[_m+3];
            by = elems_shared[_m+4];
            bz = elems_shared[_m+5];

            cx = elems_shared[_m+6];
            cy = elems_shared[_m+7];
            cz = elems_shared[_m+8];

            // nx = nodes[_j + 9];
            // ny = nodes[_j + 10];
            // nz = nodes[_j + 11];

            e01x = bx - ax;
            e01y = by - ay;
            e01z = bz - az;
            
            e02x = cx - ax;
            e02y = cy - ay;
            e02z = cz - az;

            nx = e01y * e02z - e01z * e02y;
            ny = e01z * e02x - e01x * e02z;
            nz = e01x * e02y - e01y * e02x;

    //         // e01x = nodes[_j + 12];
    //         // e01y = nodes[_j + 13];
    //         // e01z = nodes[_j + 14];
                    
    //         // e02x = nodes[_j + 15];
    //         // e02y = nodes[_j + 16];
    //         // e02z = nodes[_j + 17];

    //         //nl = nodes[_j + 18];

            D = - (nx * rdx + ny*rdy + nz*rdz);  // |-d e1 e2| = -n.d

            if(fabs(D) < 1e-6)
            {
                continue;
            }

            tx = rox - ax;
            ty = roy - ay;
            tz = roz - az;

            Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);
            u = Dy / D;

            if(u < 0 || u > 1)
            {
                continue;
            }
            // ts[_i] = elems_shared[_m];

            Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);

            v = Dz / D;

            if(v < 0 || v + u > 1)
            {
                continue;
            }

            t = (tx*nx + ty*ny + tz*nz)/D;

            if(t < 0)
            {
                continue;
            }
            if(t < min_dis)
            {
                min_dis = t;
            }
        }
        __syncthreads();
    }
    if(min_dis < 1e19)
    {
        ts[_i] = min_dis;
    }
}



int main(int argc, char** argv)
{
    camera[0] = 30;
    camera[1] = 30;
    camera[2] = 60;

    std::vector <FL_TYPE> element_vector;
    unsigned int num_of_nodes = 0;
    unsigned int num_of_elements = 0;
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
    FL_TYPE *ts = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h);
    for(size_t _i = 0; _i < w*h; _i++)
    {
        ts[_i] = 1e20;
    }
    

    FL_TYPE *drays, *dnodes, *dts;

    cudaMalloc((void **)&drays, w*h*6*sizeof(FL_TYPE));
    cudaMalloc((void **)&dts, w*h*sizeof(FL_TYPE));
    cudaMalloc((void **)&dnodes, element_vector.size()*sizeof(FL_TYPE));

    cudaMemcpy(dts, ts, w*h*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes, nodes, element_vector.size()*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
  
    FL_TYPE lights[] = {lx, ly, lz};
    // render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane);
    std::cout << "zzz\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    cudaMemcpy(drays, rays, w*h*6*sizeof(FL_TYPE), cudaMemcpyHostToDevice);

    checkintersection<<<h, w>>>(dnodes, drays, dts, (element_vector.size()/element_size));

    cudaMemcpy(ts, dts, w*h*sizeof(FL_TYPE), cudaMemcpyDeviceToHost);

    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";

    // for(size_t _i = 0; _i < h; _i++)
    // {
    //     for(size_t _j = 0; _j < w; _j++)
    //     {
    //         std::cout << ts[_i*w + _j] << " ";
    //     }
    //     std::cout << std::endl;
        
    // }
    

    // RayTrace::writeImage(image_plane, "a.ppm");
    cudaFree(drays);
    cudaFree(dts);
    cudaFree(dnodes);
}


// TODO: normalize rays