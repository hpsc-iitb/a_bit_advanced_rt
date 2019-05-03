#include <cstdint>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include <flags.hpp>
#include <utils.hpp>
#include <raytrace.hpp>
#include <domainparser.hpp>

#include <CL/cl.hpp>
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
    FL_TYPE *shadow_rays = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h * 2 * 3);
    if(!rays)
    {
        throw std::runtime_error("can't allocate memory for rays");
    }
    RayTrace::updateRays(camera, rays);

    FL_TYPE *nodes = &element_vector[0];
  
    FL_TYPE lights[] = {lx, ly, lz};
    //render(rays, nodes, element_vector.size() / element_size, lights, 1, image_plane);
    num_of_elements =  element_vector.size() / element_size;
    //#######################
    //default device of the default platform
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);
    cl::Platform default_platform=all_platforms[0];
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    cl::Device default_device=all_devices[0];
    // std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n"; // for priting the default device to be used

    cl::Context context({default_device});
 
    cl::Program::Sources sources;

    // kernel calculates for trapozidal
    std::string kernel_code= "__kernel void checkintersection(__global double* nodes,__global double* rays,__global double* ts, __global int *elem_size)"
        "{"
        "   int _i = get_global_id(0);"
        "   "
        "   double nx, ny, nz, rox, roy, roz, rdx, rdy, rdz, ax,ay,az,bx,by,bz, cx, cy, cz, e01x,e01y,e01z,e02x,e02y,e02z,nl,px,py,pz,d,d_1,tx,ty,tz,u,v,w,qx,qy,qz,modr,nnx,nny,nnz,rdnx,rdny,rdnz, D, Dy, Dz, t;"
        "   "
        "   double min_dis = 1e20;"
        "   "
        "   rox = rays[_i*6];"
        "   roy = rays[_i*6+1];"
        "   roz = rays[_i*6+2];"
        "   "    
        "   rdx = rays[_i*6+3];"
        "   rdy = rays[_i*6+4];"
        "   rdz = rays[_i*6+5];"
        "   double local_nodes[19];"
        "   "
        "   for(unsigned int _j = 0; _j < (*elem_size)*19; _j+=19)"
            "{"
            "   for(unsigned short _k = 0; _k < 19; _k++)"
            "   {"
            "       local_nodes[_k] = nodes[_j+_k];"
            "   }"
            "   ax = local_nodes[0];"
            "   ay = local_nodes[1];"
            "   az = local_nodes[2];"
            ""
            "   bx = local_nodes[3];"
            "   by = local_nodes[4];"
            "   bz = local_nodes[5];"
            "" 
            "   cx = local_nodes[6];"
            "   cy = local_nodes[7];"
            "   cz = local_nodes[8];"
            ""
            "   nx = local_nodes[  9];"
            "   ny = local_nodes[  10];"
            "   nz = local_nodes[  11];"
            ""
            "   e01x = local_nodes[  12];"
            "   e01y = local_nodes[  13];"
            "    e01z = local_nodes[  14];"
            "     "       
            "    e02x = local_nodes[  15];"
            "    e02y = local_nodes[  16];"
            "    e02z = local_nodes[  17];"
            ""
            "    nl = local_nodes[  18];"
            ""
            "    D = - nl*(nx * rdx + ny*rdy + nz*rdz);" 
            ""
            "    if(fabs(D) < 1e-6)"
            "    {"
            "        continue;"
            "    }"
            ""  
            "    tx = rox - ax;"
            "    ty = roy - ay;"
            "    tz = roz - az;"
            ""
            ""
            "    Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);"
            "    u = Dy / D;"
            ""
            "    if(u < 0 || u > 1)"
            "    {"
            "        continue;"
            "    }"
            ""
            "    Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);"
            ""
            "    v = Dz / D;"
            ""
            "    if(v < 0 || v + u > 1)"
            "    {"
            "        continue;"
            "    }"
            ""
            "    t = (tx*nx + ty*ny + tz*nz)*nl/D;"
            ""
            "    if(t < 0)"
            "    {"
            "        continue;"
            "    }"
            "    if(t < min_dis)"
            "    {"
            "        min_dis = t;"
            "    }"
            "}"
            "if(min_dis < 1e19)"
            "{"
            "    ts[_i] = min_dis;"
            "}"
        "}";

    sources.push_back({kernel_code.c_str(),kernel_code.length()});

    cl::Program program(context,sources);


    // Below is the command to check if the kernal is having any error while executing
    if(program.build({default_device})!=CL_SUCCESS){
        std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
        exit(1);
    }


    //passing a buffer form host to device
    cl::Buffer buffer_node(context,CL_MEM_READ_WRITE,sizeof(FL_TYPE)*element_vector.size());
    cl::Buffer buffer_ray(context,CL_MEM_READ_WRITE,sizeof(FL_TYPE)*w * h * 6);
    cl::Buffer buffer_ts(context,CL_MEM_READ_WRITE,sizeof(FL_TYPE)*w*h);
    cl::Buffer buffer_elem_size(context,CL_MEM_READ_WRITE,sizeof(int)*1);

    //create queue to which we will push commands for the device.
    std::cout << "check\n";
    cl::CommandQueue queue(context,default_device);
 
    FL_TYPE *t_s = (FL_TYPE *)calloc(w * h, sizeof(FL_TYPE));
    for (size_t _i = 0; _i < w*h; _i++)
    {
        t_s[_i] = 1e20;
    }
    queue.enqueueWriteBuffer(buffer_node,CL_TRUE,0,sizeof(FL_TYPE)*element_vector.size(),nodes);
    queue.enqueueWriteBuffer(buffer_elem_size,CL_TRUE,0,sizeof(int)*1,&num_of_elements);
    

    cl::Kernel kernel_rt=cl::Kernel(program,"checkintersection");
    //write arrays A and B to the device
    auto start_time = std::chrono::high_resolution_clock::now();
    queue.enqueueWriteBuffer(buffer_ray,CL_TRUE,0,sizeof(FL_TYPE)*w * h * 6,rays);
    queue.enqueueWriteBuffer(buffer_ts, CL_TRUE, 0, w*h*sizeof(FL_TYPE), t_s);

    //run the kernel
    kernel_rt.setArg(0,buffer_node);
    kernel_rt.setArg(1,buffer_ray);
    kernel_rt.setArg(2, buffer_ts);
    kernel_rt.setArg(3,buffer_elem_size);
    queue.enqueueNDRangeKernel(kernel_rt,cl::NullRange,cl::NDRange(w * h),cl::NullRange);
    queue.finish();
    queue.enqueueReadBuffer(buffer_ts,CL_TRUE,0,sizeof(FL_TYPE)*w*h,t_s);
    // for (size_t _i = 0; _i < w*h; _i++)
    // {
    //     std::cout << ((t_s[_i] > 51)?(t_s[_i]):0.0);
    // }
    //##############################
    FL_TYPE t;
    FL_TYPE rox, roy, roz;
    for (size_t _i = 0; _i < w * h * 6; _i+=6)
    {
        t = t_s[_i/6];
        
        rox = rays[_i] + rays[_i+3] * t;
        roy = rays[_i+1] + rays[_i+4] * t;
        roz = rays[_i+2] + rays[_i+5] * t;
        shadow_rays[_i] = rox;
        shadow_rays[_i+1] = roy;
        shadow_rays[_i+2] = roz;
        shadow_rays[_i+3] = lights[0] - rox;
        shadow_rays[_i+4] = lights[1] - roy;
        shadow_rays[_i+5] = lights[2] - roz;
    }

    //sources.push_back({kernel_code.c_str(),kernel_code.length()});

    //cl::Program program(context,sources);

    


    // Below is the command to check if the kernal is having any error while executing
    // if(program.build({default_device})!=CL_SUCCESS){
    //     std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
    //     exit(1);
    // }

    cl::Buffer buffer_s_ray(context,CL_MEM_READ_WRITE,sizeof(FL_TYPE)*w * h * 6);

    // cl::CommandQueue queue(context,default_device);
    for (size_t _i = 0; _i < w*h; _i++)
    {
        t_s[_i] = 1e20;
    }

    queue.enqueueWriteBuffer(buffer_s_ray,CL_TRUE,0,sizeof(FL_TYPE)*w * h * 6,shadow_rays);

    //run the kernel
    // cl::Kernel kernel_rt=cl::Kernel(program,"checkintersection");
    kernel_rt.setArg(0,buffer_node);
    kernel_rt.setArg(1,buffer_s_ray);
    kernel_rt.setArg(2, buffer_ts);
    kernel_rt.setArg(3,buffer_elem_size);
    queue.enqueueNDRangeKernel(kernel_rt,cl::NullRange,cl::NDRange(w * h),cl::NullRange);
    queue.finish();

    // FL_TYPE t_s_1[w*h];
    queue.enqueueReadBuffer(buffer_ts,CL_TRUE,0,sizeof(FL_TYPE)*w*h,t_s);
    auto end_time = std::chrono::high_resolution_clock::now();
    double time_spent = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    std::cout << "Render time taken: " << time_spent << "ms\n";

    RayTrace::writeImage(image_plane, "a.ppm");
}


// TODO: normalize rays