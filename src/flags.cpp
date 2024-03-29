#include "flags.hpp"

size_t w = 400;
size_t h = 400;

FL_TYPE w_r = 0.3;
FL_TYPE h_r = 0.3;

FL_TYPE image_plane_camera_distance = -1;

FL_TYPE pitch = M_PI / 2;
FL_TYPE yaw = - M_PI / 2;

FL_TYPE lx = 45;
FL_TYPE ly = 45;
FL_TYPE lz = 50;

// 3 coords for 3 nodes, +3 surface normal components, +6 edge components
size_t element_size = 19;

// tree tolerances
FL_TYPE tree_plus_tol = 1.2;
FL_TYPE tree_minus_tol = 0.8;

unsigned int max_depth = 5;
