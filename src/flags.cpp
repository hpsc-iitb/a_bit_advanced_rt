#include "flags.hpp"

size_t w = 100;
size_t h = 100;

FL_TYPE w_r = 0.3;
FL_TYPE h_r = 0.3;

FL_TYPE image_plane_camera_distance = -1;

FL_TYPE pitch = M_PI / 2;
FL_TYPE yaw = - M_PI / 2;

FL_TYPE lx = 45;
FL_TYPE ly = 45;
FL_TYPE lz = 60;

// 3 coords for 3 nodes, +3 surface normal components, +6 edge components
size_t element_size = 19;
