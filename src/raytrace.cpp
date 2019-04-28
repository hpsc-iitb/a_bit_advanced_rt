#include <flags.hpp>
#include <raytrace.hpp>
#include <iostream>
#include <cmath>


void render(
    FL_TYPE *rays, FL_TYPE *nodes, size_t num_elements,
    FL_TYPE *lights, size_t num_lights, FL_TYPE *image_plane
)
{
    size_t *is_hit = (size_t  *)malloc(sizeof(size_t) * w * h); // ray hit
    for (size_t _i = 0; _i < w * h; _i++)
    {
        is_hit[_i] = -1;
    }
    
    FL_TYPE *ts = (FL_TYPE *)malloc(sizeof(FL_TYPE) * w * h); // ray intersections
    FL_TYPE *illums = (FL_TYPE *)calloc(w * h * num_lights, sizeof(FL_TYPE)); // all illums

    FL_TYPE t; // parametric eqn
    FL_TYPE illum;
    FL_TYPE u, v; // barycentric coords
    FL_TYPE rox, roy, roz, rdx, rdy, rdz; // ray vector
    FL_TYPE ax, ay, az, bx, by, bz, cx, cy, cz; // nodes
    FL_TYPE px, py, pz; // intersection point
    FL_TYPE nx, ny, nz; // surface normal
    FL_TYPE e01x, e01y, e01z, e02x, e02y, e02z;
    FL_TYPE ray_length_1; // 1/shadow_ray_length

    for (size_t _i = 0; _i < w * h * 6; _i+=6)
    {
        rox = rays[_i];
        roy = rays[_i+1];
        roz = rays[_i+2];
        
        rdx = rays[_i+3];
        rdy = rays[_i+4];
        rdz = rays[_i+5];
    // std::cout << rdx << " " << rdy << " " << rdz << " " << "\n";

        FL_TYPE min_dis = 1e10;
        for (size_t _j = 0; _j < num_elements * element_size; _j+=element_size)
        {    
            ax = nodes[_j];
            ay = nodes[_j+1];
            az = nodes[_j+2];

            bx = nodes[_j+3];
            by = nodes[_j+4];
            bz = nodes[_j+5];

            cx = nodes[_j+6];
            cy = nodes[_j+7];
            cz = nodes[_j+8];

            nx = nodes[_j + 9];
            ny = nodes[_j + 10];
            nz = nodes[_j + 11];

            e01x = nodes[_j + 12];
            e01y = nodes[_j + 13];
            e01z = nodes[_j + 14];
            
            e02x = nodes[_j + 15];
            e02y = nodes[_j + 16];
            e02z = nodes[_j + 17];
            
            if(checkIntersection(
                rox, roy, roz, rdx, rdy, rdz, ax, ay, az,
                bx, by, bz, cx, cy, cz, nx, ny, nz,
                e01x, e01y, e01z, e02x, e02y, e02z,
                t, illum, u, v
            ))
            {
                px = rox + t*rdx;
                py = roy + t*rdy;
                pz = roz + t*rdz;
                // optimisation, use norm
                FL_TYPE dis = distance(
                    px, py, pz, rox, roy, roz
                );
                if(dis < min_dis)
                {
                    is_hit[_i/6] = _j/element_size;
                    ts[_i/6] = t;
                    min_dis = dis;
                    image_plane[_i/6] = fabs(rdx * nx + rdy * ny + rdz * nz)/20;
                }
            }
        }
    }

    std::cout << "Primary ray hit complete\n";

    return;
    // calculate illumination
    for (size_t _j = 0; _j < w * h * 6; _j+=6)
    {
        if(is_hit[_j/6] == -1)
        {
            image_plane[_j] = 0.0;
            continue;
        }
        
        t = ts[_j/6];
        int eid = is_hit[_j/6];
        nx = nodes[eid * element_size + 9];
        ny = nodes[eid * element_size + 10];
        nz = nodes[eid * element_size + 11];

        rox = rays[_j] + rays[_j+3] * t;
        roy = rays[_j+1] + rays[_j+4] * t;
        roz = rays[_j+2] + rays[_j+5] * t;
        for (size_t _i = 0; _i < num_lights; _i++)
        {
            FL_TYPE lx = lights[_i * 3];
            FL_TYPE ly = lights[_i * 3 + 1];
            FL_TYPE lz = lights[_i * 3 + 2];

            
            rdx = lx - rox;
            rdy = ly - roy;
            rdz = lz - roz;

            ray_length_1 = 1/sqrt(rdx * rdx + rdy * rdy + rdz * rdz);
            rdx *= ray_length_1;
            rdy *= ray_length_1;
            rdz *= ray_length_1;

            illum = fabs(rdx * nx + rdy * ny + rdz * nz);
            // illums[_i * w * h + _j/6] = 0.5;

        }        
    }
    
    std::cout << "Done calculating illumination\n";


    // now do the shadow rays
    // for (size_t _i = 0; _i < h * w; _i++)
    // {
    //     std::cout << is_hit[_i] << " ";
    // }
    
    // for(size_t _i = 0; _i < h; _i++)
    // {
    //     for(size_t _j = 0; _j < w; _j++)
    //     {
    //         std::cout << image_plane[_i * w + _j] << " ";
    //     }
    //     std::cout << "\n";
    // }


    for (size_t _i = 0; _i < w * h * 6; _i+=6)
    {
        if (is_hit[_i/6] == -1)
        {
            // no planes hit
            image_plane[_i/6] = 0.0;
            continue;
        }
        t = ts[_i/6];
        
        rox = rays[_i] + rays[_i+3] * t;
        roy = rays[_i+1] + rays[_i+4] * t;
        roz = rays[_i+2] + rays[_i+5] * t;
        
        // TODO: multiple lights: done
        for (size_t _k = 0; _k < num_lights; _k++)
        {

            rdx = lights[_k * 3 + 0] - rox;
            rdy = lights[_k * 3 + 1] - roy;
            rdz = lights[_k * 3 + 2] - roz;

            // image_plane[_i/6] = 1.0;
            
            for (size_t _j = 0; _j < num_elements * element_size; _j+=element_size)
            {   
                if(_j / element_size == is_hit[_i / 6]) // stop self hits
                {
                    continue;
                }
                ax = nodes[_j];
                ay = nodes[_j+1];
                az = nodes[_j+2];

                bx = nodes[_j+3];
                by = nodes[_j+4];
                bz = nodes[_j+5];

                cx = nodes[_j+6];
                cy = nodes[_j+7];
                cz = nodes[_j+8];

                nx = nodes[_j + 9];
                ny = nodes[_j + 10];
                nz = nodes[_j + 11];

                e01x = nodes[_j + 12];
                e01y = nodes[_j + 13];
                e01z = nodes[_j + 14];
                
                e02x = nodes[_j + 15];
                e02y = nodes[_j + 16];
                e02z = nodes[_j + 17];

                if(checkIntersection(
                    rox, roy, roz, rdx, rdy, rdz, ax, ay, az,
                    bx, by, bz, cx, cy, cz, nx, ny, nz,
                    e01x, e01y, e01z, e02x, e02y, e02z,
                    t, illum, u, v
                ))
                {
                    if(t < 0)
                    {
                        // image_plane[_i/6] = 0.0;
                        illums[_k * w * h + _i/6] = 0; // no illuminaiton due to this light
                        std::cout << "shadowed";
                    }
                    
                    // std::cout << roz << " " << rdz << " " << ax << " " << ay << " " << az << " " << t << " " << u << " " << v <<  "\n";
                    break;
                }
            }
        }
    }

    for (size_t _i = 0; _i < w * h * num_lights; _i++)
    {
        image_plane[_i % (w*h)] += illums[_i];
    }
    
}

inline bool checkIntersection(
     FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
     FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
     FL_TYPE ax, FL_TYPE ay, FL_TYPE az,
     FL_TYPE bx, FL_TYPE by, FL_TYPE bz,
     FL_TYPE cx, FL_TYPE cy, FL_TYPE cz,
     FL_TYPE nx, FL_TYPE ny, FL_TYPE nz,
     FL_TYPE e01x, FL_TYPE e01y, FL_TYPE e01z,
     FL_TYPE e02x, FL_TYPE e02y, FL_TYPE e02z,
     FL_TYPE &t, FL_TYPE &illum, FL_TYPE &u,
     FL_TYPE &v
)
{   // calc the edges
#ifndef USE_PRECOMPUTED_NORMALS
    // FL_TYPE e01x = bx - ax;
    // FL_TYPE e01y = by - ay;
    // FL_TYPE e01z = bz - az;

    // FL_TYPE e02x = cx - ax;
    // FL_TYPE e02y = cy - ay;
    // FL_TYPE e02z = cz - az;

    FL_TYPE px, py, pz; // _nx, _ny, _nz;
    cross(rdx, rdy, rdz, e02x, e02y, e02z, px, py, pz);
    
    // calculate the triple product
    FL_TYPE d = e01x * px + e01y * py + e01z * pz;

    // check for paralle surface
    if(fabs(d) < 1e-7)
    {
        return false;
    }

    // illum = fabs(d/20.0);

    FL_TYPE d_1 = 1.0 / d;

    FL_TYPE tx = rox - ax;
    FL_TYPE ty = roy - ay;
    FL_TYPE tz = roz - az;

    u = tx * px + ty * py + tz * pz;
    u *= d_1;
    if(u < 0 || u > 1)
    {
        return false;
    }

    FL_TYPE qx, qy, qz;
    cross(
        tx, ty, tz, e01x, e01y, e01z, qx, qy, qz
    );

    v = qx * rdx + qy * rdy + qz * rdz;
    v *= d_1;
    if(v < 0 || v + u > 1)
    {
        return false;
    }

    t = e02x * qx + e02y * qy + e02z * qz;
    t *= d_1;
    
    // cross(e01x, e01y, e01z, e02x, e02y, e02z, _nx, _ny, _nz);
    FL_TYPE rdnx, rdny, rdnz, nnx, nny, nnz;
    normalize(rdx, rdy, rdz, rdnx, rdny, rdnz);
    normalize(nx, ny, nz, nnx, nny, nnz);

    illum = -rdnx * nnx - rdny * nny - rdnz * nnz;
    // std::cout <<"ray dir: " <<  rdx << " " << rdy << " " << rdz << " " << "\n";
    // std::cout <<"normal: " <<  nx << " " << ny << " " << nz << " " << "\n";
    // std::cout << "illum: " << illum << "\n";
    illum = fabs(illum);
#else
    FL_TYPE D = - (nx * rdx + ny*rdy + nz*rdz);  // |-d e1 e2| = -n.d

    if(fabs(D) < 1e-6)
    {
        return false;
    }

    FL_TYPE tx = rox - ax;
    FL_TYPE ty = roy - ay;
    FL_TYPE tz = roz - az;
    // FL_TYPE px, py, pz;
    // cross(rdx, rdy, rdz, e02x, e02y, e02z, px, py, pz);

    FL_TYPE Dy = rdx*(tz*e02y - ty*e02z) + tx*(rdy*e02z - e02y*rdz) + e02x*(ty*rdz - rdy*tz);
    u = Dy / D;

    if(u < 0 || u > 1)
    {
        return false;
    }

    FL_TYPE Dz = rdx*(e01z*ty - e01y*tz) + e01x*(tz*rdy - ty*rdz) + tx*(e01y*rdz - rdy*e01z);

    v = Dz / D;

    if(v < 0 || v + u > 1)
    {
        return false;
    }

    t = (tx*nx + ty*ny + tz*nz)/D;
#endif
    return true;
}


inline void cross(
    FL_TYPE &ax, FL_TYPE &ay, FL_TYPE &az,
    FL_TYPE &bx, FL_TYPE &by, FL_TYPE &bz,
    FL_TYPE &cx, FL_TYPE &cy, FL_TYPE &cz
)
{
    cx = ay * bz - az * by;
    cy = az * bx - ax * bz;
    cz = ax * by - ay * bx;
}


inline FL_TYPE distance(
    FL_TYPE ax, FL_TYPE ay, FL_TYPE az,
    FL_TYPE bx, FL_TYPE by, FL_TYPE bz
)
{
    return sqrt(
        (bx - ax) * (bx - ax) + (by - ay) * (by - ay) + (bz - az) * (bz - az)
    );
}

inline void normalize(
    FL_TYPE &x, FL_TYPE &y, FL_TYPE &z,
    FL_TYPE &nx, FL_TYPE &ny, FL_TYPE &nz
)
{
    FL_TYPE l = sqrt(x * x + y * y + z * z);
    nx = x / l;
    ny = y / l;
    nz = z / l;
}