#ifndef COORD_HPP
#define COORD_HPP

#include <vector>
#include <cmath>
#include <stdexcept>

#include "flags.hpp"

namespace RayTrace{
    class Coordinate
    {
        public:
        Coordinate(FL_TYPE x = 0, FL_TYPE y = 0, FL_TYPE z = 0);
        Coordinate(std::vector<FL_TYPE> &cd);
        void normalize();
        void update();
        void cross(Coordinate &rt);
        void dot(Coordinate &rt);
        static Coordinate cross(Coordinate &lt, Coordinate &rt);
        static Coordinate dot(Coordinate &lt, Coordinate &rt);
        static Coordinate normalize(Coordinate &cd);
        FL_TYPE x, y, z;
        FL_TYPE length;
        FL_TYPE nx, ny, nz; // normalized unit vector
    };
}
#endif