#include "coord.hpp"

#include <vector>
#include <cmath>
#include <stdexcept>

namespace RayTrace{
    Coordinate::Coordinate(FL_TYPE x, FL_TYPE y, FL_TYPE z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->update();
    }

    Coordinate::Coordinate(std::vector<FL_TYPE> &cd)
    {
        if(cd.size() != 3)
        {
            throw std::invalid_argument("vector must have exactly 3 values");
        }
        this->x = cd[0];
        this->y = cd[1];
        this->z = cd[2];
        this->update();
    }

    void Coordinate::update()
    {
        this->length = sqrt(
            this->x * this->x + this->y * this->y + this->z * this->z
        );
        this->normalize();
    }

    void Coordinate::normalize()
    {
        Coordinate::normalize(*this);
    }

    void Coordinate::cross(Coordinate &rt)
    {
        Coordinate::cross(*this, rt);
    }

    void Coordinate::dot(Coordinate &rt)
    {
        Coordinate::dot(*this, rt);
    }

    Coordinate Coordinate::normalize(Coordinate &cd)
    {
        cd.nx = cd.x / cd.length;
        cd.ny = cd.y / cd.length;
        cd.nz = cd.z / cd.length;

        return cd;
    }

    Coordinate Coordinate::cross(Coordinate &lt, Coordinate &rt)
    {
        // TODO: implement
        return lt;
    }

    Coordinate Coordinate::dot(Coordinate &lt, Coordinate &rt)
    {
        // TODO: implement
        return lt;
    }
}