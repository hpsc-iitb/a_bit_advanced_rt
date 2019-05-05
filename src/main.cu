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
void primary_rays(float* nodes,float* rays, float* image_plane, float* lights, size_t elem_size)
{
    int _i = blockIdx.x * blockDim.x + threadIdx.x;
    // int _j = get_global_id(1);
    float ax,ay,az,bx,by,bz,cx,cy,cz,e01x,e01y,e01z,e02x,e02y,e02z,\
    tx,ty,tz,u,v,D, Dy, Dz, t, nx, ny, nz, rox, roy, roz, rdx, rdy, rdz, ray_length_1,\
    srox, sroy, sroz, srdx, srdy, srdz, nl_1;

    float min_dis = 1e20;
    int near_ele_num = -1, eid;

    unsigned int sz;//, eidx;

    FL_TYPE lx = lights[0];
    FL_TYPE ly = lights[1];
    FL_TYPE lz = lights[2];

    
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
            //eidx = _j + _m;
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
                near_ele_num = _j+_m;
                min_dis = t;
            }
        }
        __syncthreads();
    }
    
    // if(min_dis<1e19)
    // {
    //     is_hit[_i] = near_ele_num;
    //     ts[_i] = min_dis;
    // }

    //for (size_t _j = 0; _j < w * h * 6; _j+=6)
    //{

        //for (size_t _j = 0; _j < w * h * 6; _j+=6)
    //{
    if(near_ele_num == -1)
    {
        image_plane[_i] = 0.0;
    }
    else
    {
        eid = near_ele_num;
        ax = nodes[eid];
        ay = nodes[eid+1];
        az = nodes[eid+2];

        bx = nodes[eid+3];
        by = nodes[eid+4];
        bz = nodes[eid+5];

        cx = nodes[eid+6];
        cy = nodes[eid+7];
        cz = nodes[eid+8];

        e01x = bx - ax;
        e01y = by - ay;
        e01z = bz - az;
        
        e02x = cx - ax;
        e02y = cy - ay;
        e02z = cz - az;

        nx = e01y * e02z - e01z * e02y;
        ny = e01z * e02x - e01x * e02z;
        nz = e01x * e02y - e01y * e02x;

        srox = rox + rdx * min_dis;
        sroy = roy + rdy * min_dis;
        sroz = roz + rdz * min_dis;

        srdx = lx - srox;
        srdy = ly - sroy;
        srdz = lz - sroz;
        
        ray_length_1 = 1/sqrt(srdx * srdx + srdy * srdy + srdz * srdz);
        nl_1 = 1/sqrt(nx * nx + ny * ny + nz * nz);
        // srdx *= ray_length_1;
        // srdy *= ray_length_1;
        // srdz *= ray_length_1;

        image_plane[_i] = fabs(srdx * nx + srdy * ny + srdz * nz)*ray_length_1*nl_1;
    }
    // return;  
    __syncthreads();

    for(unsigned int _p = 0; _p < elem_size*9; _p+=100*9)
    {
        sz = ((elem_size*9 - _p)>=100*9)?100*9:(elem_size*9 - _p);
        if(threadIdx.x == 0)
        {
            for(unsigned int _q = 0; _q < sz; _q++)
            {
                elems_shared[_q] = nodes[_p + _q];
            }
        }
        __syncthreads();

        // srox = rox + rdx * min_dis;
        // sroy = roy + rdy * min_dis;
        // sroz = roz + rdz * min_dis;
        
        // rdx = lx - srox;
        // rdy = ly - sroy;
        // rdz = lz - sroz;

        for(unsigned int _r = 0; _r < sz; _r+= 9)
        {
            //eidx = _j + _r;
            if(_p + _r == near_ele_num)
            {
                // image_plane[0] = 0;
                continue;
            }
            ax = elems_shared[_r];
            ay = elems_shared[_r+1];
            az = elems_shared[_r+2];

            bx = elems_shared[_r+3];
            by = elems_shared[_r+4];
            bz = elems_shared[_r+5];

            cx = elems_shared[_r+6];
            cy = elems_shared[_r+7];
            cz = elems_shared[_r+8];

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

            D = - (nx * srdx + ny*srdy + nz*srdz);  // |-d e1 e2| = -n.d

            if(fabs(D) < 1e-6)
            {
                continue;
            }

            tx = srox - ax;
            ty = sroy - ay;
            tz = sroz - az;

            Dy = srdx*(tz*e02y - ty*e02z) + tx*(srdy*e02z - e02y*srdz) + e02x*(ty*srdz - srdy*tz);
            u = Dy / D;

            if(u < 0 || u > 1)
            {
                continue;
            }
            // ts[_i] = elems_shared[_m];

            Dz = srdx*(e01z*ty - e01y*tz) + e01x*(tz*srdy - ty*srdz) + tx*(e01y*srdz - srdy*e01z);

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
                image_plane[_i] = 0.0;
                return;
            }
        }
        __syncthreads(); 
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
    FL_TYPE *shadow_rays = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h * 2 * 3);
    if(!rays)
    {
        throw std::runtime_error("can't allocate memory for rays");
    }
    RayTrace::updateRays(camera, rays);

    FL_TYPE *nodes = &element_vector[0];
    FL_TYPE *ts = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h);
    FL_TYPE *is_hit = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h);
    for(size_t _i = 0; _i < w*h; _i++)
    {
        ts[_i] = 1e20;
        is_hit[_i] = -1;
    }
    
    FL_TYPE *drays, *dnodes, *dimage_plane, *dlights;//, *dts;
    FL_TYPE lights[] = {lx, ly, lz};

    cudaMalloc((void **)&drays, w*h*6*sizeof(FL_TYPE));
    cudaMalloc((void **)&dimage_plane, w*h*sizeof(FL_TYPE));
    cudaMalloc((void **)&dlights, 3*sizeof(FL_TYPE));
    cudaMalloc((void **)&dnodes, element_vector.size()*sizeof(FL_TYPE));

    cudaMemcpy(dimage_plane, image_plane, w*h*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(dlights, lights, 3*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(dnodes, nodes, element_vector.size()*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(drays, rays, w*h*6*sizeof(FL_TYPE), cudaMemcpyHostToDevice);

    // render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane);
    //std::cout << "zzz\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    primary_rays<<<h, w>>>(dnodes, drays, dimage_plane, dlights, (element_vector.size()/element_size));

    cudaMemcpy(image_plane, dimage_plane, w*h*sizeof(FL_TYPE), cudaMemcpyDeviceToHost);

    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";
    //cudaMemcpy(drays, shadow_rays, w*h*6*sizeof(FL_TYPE), cudaMemcpyHostToDevice);
    

    // for(size_t _i = 0; _i < h; _i++)
    // {
    //     for(size_t _j = 0; _j < w; _j++)
    //     {
    //         std::cout << image_plane[_i*w + _j] << " ";
    //     }
    //     std::cout << std::endl;
        
    // }
    

    RayTrace::writeImage(image_plane, "a.ppm");
    cudaFree(drays);
    cudaFree(dlights);
    cudaFree(dimage_plane);
    cudaFree(dnodes);
}
