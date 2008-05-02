#ifndef INTERSECT_H_INCLUDED
#define INTERSECT_H_INCLUDED

#include "util.h"
#include <vector>

struct intersect_info_t {
    typedef Boundary::edge_iterator edge_iterator;

    vec_t ivtx;
    double inc;
    edge_iterator iedge;
};

typedef std::vector <intersect_info_t> intersect_buffer_t;

bool intersect_ray_boundary (intersect_info_t *,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *);

int  intersect_ray_boundary (intersect_buffer_t *,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *);

int intersect_line_boundary (intersect_buffer_t *,
                             const vec_t &r0, const vec_t &dir,
                             Boundary *);

#endif // INTERSECT_H_INCLUDED
