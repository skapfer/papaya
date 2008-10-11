#ifndef PAPAYA_GEOMETRY_HEADER 
#define PAPAYA_GEOMETRY_HEADER 

#include "util.h"

namespace PapayaGeometry {

    vec3_t plane_normal (const vec3_t &a, const vec3_t &b, const vec3_t &c);
    double triangle_area (const vec3_t &a, const vec3_t &b, const vec3_t &c);
    vec3_t average (const vec3_t &a, const vec3_t &b, const vec3_t &c);

}

#endif // PAPAYA_GEOMETRY_HEADER
