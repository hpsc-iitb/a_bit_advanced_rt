#include "utils.hpp"

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>

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
        size_t c = img.at[0].size();
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
}