#ifndef RAYTRACE_HPP
#define RAYTRACE_HPP
#include <flags.hpp>

void render(
    FL_TYPE *rays, FL_TYPE *nodes, size_t num_elements,
    FL_TYPE *lights, size_t num_lights, FL_TYPE *image_plane
);

inline bool checkIntersection(
     FL_TYPE rox, FL_TYPE roy, FL_TYPE roz,
     FL_TYPE rdx, FL_TYPE rdy, FL_TYPE rdz,
     FL_TYPE ax, FL_TYPE ay, FL_TYPE az,
     FL_TYPE bx, FL_TYPE by, FL_TYPE bz,
     FL_TYPE cx, FL_TYPE cy, FL_TYPE cz,
     FL_TYPE &t, FL_TYPE &illum, FL_TYPE &u,
     FL_TYPE &v
);

inline void cross(
    FL_TYPE &ax, FL_TYPE &ay, FL_TYPE &az,
    FL_TYPE &bx, FL_TYPE &by, FL_TYPE &bz,
    FL_TYPE &cx, FL_TYPE &cy, FL_TYPE &cz
);

inline FL_TYPE distance(
    FL_TYPE ax, FL_TYPE ay, FL_TYPE az,
    FL_TYPE bx, FL_TYPE by, FL_TYPE bz
);

inline void normalize(
    FL_TYPE &x, FL_TYPE &y, FL_TYPE &z,
    FL_TYPE &nx, FL_TYPE &ny, FL_TYPE &nz
);
#endif