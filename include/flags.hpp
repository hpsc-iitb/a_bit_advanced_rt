#pragma once

#ifndef FLAGS_HPP
#define FLAGS_HPP

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159
#endif

typedef double FL_TYPE; // typedef to change types quickly

extern size_t w; // width of image, pixels
extern size_t h; // height of image, pixels

extern FL_TYPE w_r; // width of image plane, real
extern FL_TYPE h_r; // width of image plane, real

extern FL_TYPE image_plane_camera_distance; // distance of image plane from eye

extern FL_TYPE pitch; // inclination from the horizontal, about x axis
extern FL_TYPE yaw; // rotation about y axis

// oner light source
extern FL_TYPE lx;
extern FL_TYPE ly;
extern FL_TYPE lz;

extern size_t element_size; // size of an element in element_vector

#define USE_PRECOMPUTED_NORMALS

#endif