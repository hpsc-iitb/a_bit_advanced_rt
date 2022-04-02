#include <utils.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>

namespace RayTrace{

    // void writeImage(std::vector<std::vector<FL_TYPE>> &img, std::string &fname, int max_val)
    void writeImage(FL_TYPE *img, std::string fname, int max_val)
    {
        /* Write a greyscale PPM image
        Header format:
        P6
        {width} {height}
        {maxval}
        {r g b r g b ...}
        */
        // size_t r = img.size();
        // if (r < 1)
        // {
        //     throw std::invalid_argument("image row size < 1");
        // }
        // size_t c = img.at(0).size();
        // if (c < 1)
        // {
        //     throw std::invalid_argument("image col size < 1");
        // }

        // open the file
        std::ofstream image_file(fname, std::ios::binary);
        // use 8bits
        image_file << "P6\n" << h << " " << w << "\n" << max_val << "\n";
        for (size_t _i = 0; _i < h; _i++)
        {
            for (size_t _j = 0; _j < w; _j++)
            {
                uint8_t r, g, b;
                // bound pixel intensity in [0, 1]
                FL_TYPE pp = img[_i * w + _j];
                pp = (pp > 1.0)?1.0:pp;
                pp = (pp < 0.0)?0.0:pp;

                r = g = b = (uint8_t) 255 * pp;
                image_file << r << g << b;
            }
        }
        image_file.close();
    }

    void updateRays(std::vector<FL_TYPE> &camera, FL_TYPE *ray)
    {
        FL_TYPE cx = camera[0];
        FL_TYPE cy = camera[1];
        FL_TYPE cz = camera[2];
        
        FL_TYPE dx = w_r / w;
        FL_TYPE dy = h_r / h;

        FL_TYPE dxb2 = dx * 0.5;
        FL_TYPE dyb2 = dy * 0.5;

        FL_TYPE xmin = - 0.5 * w_r + dxb2;
        FL_TYPE ymin = - 0.5 * h_r + dyb2;

        for(size_t _i = 0; _i < h; _i++)
        {
            size_t hidx = _i * w * 6;
            for(size_t _j = 0; _j < w; _j++)
            {
                size_t idx = hidx + _j * 6;
                ray[idx] = cx;
                ray[idx+1] = cy;
                ray[idx+2] = cz;
                // std::cout << idx << "\n";

                FL_TYPE x_i_0 = _j * dx + xmin; 
                FL_TYPE y_i_0 = _i * dy + ymin;
                FL_TYPE z_i_0 = image_plane_camera_distance;
                // std::cout << x_i_0 << " " << y_i_0 << " " << z_i_0 << " " << "\n";

                // FL_TYPE z_i_0 = image_plane_camera_distance * cos(pitch)\
                //     * cos(yaw) + x_i_0 * sin(yaw) + y_i_0 * sin(pitch);

                // x_i_0 *= cos(yaw);
                // y_i_0 *= cos(pitch);

                // x_i_0 += image_plane_camera_distance * sin(yaw) * cos(pitch);
                // y_i_0 += image_plane_camera_distance * sin(pitch);

                ray[idx+3] = x_i_0;
                ray[idx+4] = y_i_0;
                ray[idx+5] = z_i_0;
            }
        }
        
    }

}

void normalize(
    FL_TYPE &x, FL_TYPE &y, FL_TYPE &z,
    FL_TYPE &nx, FL_TYPE &ny, FL_TYPE &nz
)
{
    FL_TYPE l = sqrt(x * x + y * y + z * z);
    nx = x / l;
    ny = y / l;
    nz = z / l;
}