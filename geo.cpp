
#define PAPAYA2
#ifdef DIMENSION
#error "geo.cpp does not use DIMENSION."
#endif
#include "geo.h"


namespace PapayaGeometry {

    // FIXME robustness

    vec3_t plane_normal (const vec3_t &a, const vec3_t &b, const vec3_t &c) {
        vec3_t X = a-b;
        vec3_t Y = a-c;
        vec3_t Z = cross (X, Y);
        Z.normalize ();
        return Z;
    }

    double triangle_area (const vec3_t &a, const vec3_t &b, const vec3_t &c) {
        vec3_t X = a-b;
        vec3_t Y = a-c;
        vec3_t Z = cross (X, Y);
        return Z.norm () / 2;
    }

    vec3_t average (const vec3_t &a, const vec3_t &b, const vec3_t &c) {
        vec3_t X = a;
        X += b;
        X += c;
        return X/3;
    }
}

