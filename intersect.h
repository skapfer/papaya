#ifndef INTERSECT_H_INCLUDED
#define INTERSECT_H_INCLUDED

#include "util.h"
#include <vector>

struct intersect_info_t {
    typedef Boundary::edge_iterator edge_iterator;
    typedef intersect_info_t this_t;

    vec_t ivtx;
    double inc;
    edge_iterator iedge;

    static bool by_contour_id (const this_t &, const this_t &);
    static bool by_normal_coordinate (const this_t &, const this_t &);
};

typedef std::vector <intersect_info_t> intersect_buffer_t;

// save the intersect points to a GNUPLOT-readable format
void dump_intersect_buffer (std::ostream &, const intersect_buffer_t &);

bool intersect_ray_boundary (intersect_info_t *,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *);

int  intersect_ray_boundary (intersect_buffer_t *,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *);

int intersect_line_boundary (intersect_buffer_t *,
                             const vec_t &r0, const vec_t &dir,
                             Boundary *);

bool intersect_vertex_rect (const vec_t &v, const rect_t &r);

//
// inline implementation
//
inline bool intersect_info_t::by_contour_id (
    const intersect_info_t &lhs,
    const intersect_info_t &rhs) {
    return lhs.iedge->contour < rhs.iedge->contour;
}

inline bool intersect_info_t::by_normal_coordinate (
    const intersect_info_t &lhs,
    const intersect_info_t &rhs) {
    return lhs.inc < rhs.inc;
}

#endif // INTERSECT_H_INCLUDED
