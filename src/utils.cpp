#include "utils.hpp"

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>

namespace RayTrace{

    void writeImage(std::vector<std::vector<FL_TYPE>> &img, std::string &fname, int max_val)
    {
        /* Write a greyscale PPM image
        Header format:
        P6
        {width} {height}
        {maxval}
        {r g b r g b ...}
        */
        size_t r = img.size();
        if (r < 1)
        {
            throw std::invalid_argument("image row size < 1");
        }
        size_t c = img.at(0).size();
        if (c < 1)
        {
            throw std::invalid_argument("image col size < 1");
        }

        // open the file
        std::ofstream image_file(fname);
        // use 8bits
        image_file << "P6\n" << r << " " << c << "\n" << max_val << "\n";
        for (size_t _i = 0; _i < c; _i++)
        {
            for (size_t _j = 0; _j < r; _j++)
            {
                uint8_t r, g, b;
                // bound pixel intensity in [0, 1]
                FL_TYPE pixel_val = (img[r][c] > 1.0)?1.0:img[r][c];
                pixel_val = (pixel_val < 0.0)?0.0:pixel_val;

                r = g = b = (uint8_t) 255 * pixel_val;
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

                FL_TYPE x_i_0 = _j * dx + xmin; 
                FL_TYPE y_i_0 = _i * dy + ymin;
                
                FL_TYPE z_i_0 = image_plane_camera_distance * cos(pitch)\
                    * cos(yaw);

                x_i_0 *= cos(yaw);
                y_i_0 *= cos(pitch);

                x_i_0 += image_plane_camera_distance * sin(yaw) * cos(pitch);
                y_i_0 += image_plane_camera_distance * sin(pitch);

                ray[idx+3] = x_i_0;
                ray[idx+4] = y_i_0;
                ray[idx+5] = z_i_0;
            }
        }
        
    }

}