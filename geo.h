#ifndef PAPAYA_GEOMETRY_HEADER 
#define PAPAYA_GEOMETRY_HEADER 

#include "util.h"

namespace PapayaGeometry {

    // find the normal of a plane defined by three points.
    vec3_t plane_normal (const vec3_t &a, const vec3_t &b, const vec3_t &c);

    // find the area of the triangle defined by three vertices.
    double triangle_area (const vec3_t &a, const vec3_t &b, const vec3_t &c);

    // average of three points
    vec3_t average (const vec3_t &a, const vec3_t &b, const vec3_t &c);

}

#endif // PAPAYA_GEOMETRY_HEADER
