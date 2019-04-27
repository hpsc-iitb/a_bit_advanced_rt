#include <flags.hpp>
#include <raytrace.hpp>
#include <iostream>


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
    FL_TYPE t; // parametric eqn
    FL_TYPE illum;
    FL_TYPE u, v; // barycentric coords
    FL_TYPE rox, roy, roz, rdx, rdy, rdz;
    FL_TYPE ax, ay, az, bx, by, bz, cx, cy, cz;
    FL_TYPE px, py, pz;
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
        for (size_t _j = 0; _j < num_elements * 9; _j+=9)
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
            // if(_i/6 == 29)
            // {
            //     std::cout << "29\n"; 
            // }
            if(checkIntersection(
                rox, roy, roz, rdx, rdy, rdz, ax, ay, az,
                bx, by, bz, cx, cy, cz, t, illum, u, v
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
                    is_hit[_i/6] = _j/9;
                    ts[_i/6] = t;
                    min_dis = dis;
                    image_plane[_i/6] = illum;
                }
            }
        }
    }
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
        
        // TODO: multiple lights
        rdx = lights[0] - rox;
        rdy = lights[1] - roy;
        rdz = lights[2] - roz;

        // image_plane[_i/6] = 1.0;
        
        for (size_t _j = 0; _j < num_elements * 9; _j+=9)
        {   
            if(_j / 9 == is_hit[_i / 6]) // stop self hits
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
            if(checkIntersection(
                rox, roy, roz, rdx, rdy, rdz, ax, ay, az,
                bx, by, bz, cx, cy, cz, t, illum, u, v
            ))
            {
                if(t < 0)
                {
                    image_plane[_i/6] = 0.0;
                }
                
                // std::cout << roz << " " << rdz << " " << ax << " " << ay << " " << az << " " << t << " " << u << " " << v <<  "\n";
                break;
            }
        }
    }
}

inline bool checkIntersection(
     FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
     FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
     FL_TYPE ax, FL_TYPE ay, FL_TYPE az,
     FL_TYPE bx, FL_TYPE by, FL_TYPE bz,
     FL_TYPE cx, FL_TYPE cy, FL_TYPE cz,
     FL_TYPE &t, FL_TYPE &illum, FL_TYPE &u,
     FL_TYPE &v
)
{   // calc the edges
    FL_TYPE e01x = bx - ax;
    FL_TYPE e01y = by - ay;
    FL_TYPE e01z = bz - az;

    FL_TYPE e02x = cx - ax;
    FL_TYPE e02y = cy - ay;
    FL_TYPE e02z = cz - az;

    FL_TYPE px, py, pz, nx, ny, nz;
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
    
    cross(e01x, e01y, e01z, e02x, e02y, e02z, nx, ny, nz);
    FL_TYPE rdnx, rdny, rdnz, nnx, nny, nnz;
    normalize(rdx, rdy, rdz, rdnx, rdny, rdnz);
    normalize(nx, ny, nz, nnx, nny, nnz);

    illum = -rdnx * nnx - rdny * nny - rdnz * nnz;
    // std::cout <<"ray dir: " <<  rdx << " " << rdy << " " << rdz << " " << "\n";
    // std::cout <<"normal: " <<  nx << " " << ny << " " << nz << " " << "\n";
    // std::cout << "illum: " << illum << "\n";
    illum = fabs(illum);
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