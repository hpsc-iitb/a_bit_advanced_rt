#pragma once

#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstdint>
#include <cmath>

#include "flags.hpp"


namespace RayTrace{

    void writeImage(std::vector<std::vector<FL_TYPE>> &img, std::string &fname, int max_val = 255);
    // + 2 overloads

    // void 

    void updateRays(std::vector<FL_TYPE> &camera, FL_TYPE *ray);

}
#endif