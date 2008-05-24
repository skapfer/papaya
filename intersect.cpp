// vim:et:sw=4:ts=4
// simple N^2 intersection algorithm

#include "intersect.h"


namespace {
    enum {
        MODE_RAY,
        MODE_LINE,
        MODE_RAY_NEAREST
    };

#define b() b
#define edge_iter() eit
#define cdir() cdir_
#define dir() dir_
    bool does_edge_intersect (intersect_info_t *info,
                              Boundary *b,
                              Boundary::edge_iterator eit,
                              const vec_t &r0, const vec_t &dir_) {
        vec_t cdir_ = rot90_ccw (dir_);
        // determine lambda so that v0 + lambda(v1-v0) in ray
        double cc0 = dot (cdir (), b ()->edge_vertex0 (edge_iter ()));
        double cc1 = dot (cdir (), b ()->edge_vertex1 (edge_iter ()));
        double ccR = dot (cdir (), r0);
        assert (!isnan (cc0));
        assert (!isnan (cc1));
        assert (!isnan (ccR));
        double lambda;
        if (cc0 <= ccR and ccR < cc1) {
            double l1 = ccR - cc0;
            double l2 = cc1 - cc0;
            if (l2 <= l1)        // keep underflows in check
                lambda = 1.;     // (exact value is irrelevant)
            else
                lambda = l1/l2;
        } else if (cc1 < ccR and ccR <= cc0) {
            double l1 = cc0 - ccR;
            double l2 = cc0 - cc1;
            if (l2 <= l1)
                lambda = 1.;
            else
                lambda = l1/l2;
        } else {
            return false;
        }
        assert (lambda <= 1.);
        assert (lambda >= 0.);
        // range of normal coordinates of current edge
        double nc0 = dot (dir (), b ()->edge_vertex0 (edge_iter ()));
        double ncT = dot (dir (), b ()->edge_vertex1 (edge_iter ()));
        double mu = nc0;
        ncT -= nc0;
        mu += lambda * ncT;
        double ncR = dot (dir (), r0);
        mu -= ncR;
        // now: mu: dir r0 + mu = dir (v0 + lambda(v1-v0))
        // update vertex-of-intersection since we have calculated
        // all the necessary things anyway
        vec_t &ivtx = info->ivtx;
        ivtx  = (1-lambda) * b ()->edge_vertex0 (edge_iter ());
        ivtx += lambda     * b ()->edge_vertex1 (edge_iter ());
        info->inc = mu;
        info->iedge = eit;
        return true;
    }
}
#undef b
#undef edge_iter
#undef dir
#undef cdir

void intersect_ray_boundary_impl (intersect_buffer_t *dst,
                                  const vec_t &ray_0, const vec_t &ray_dir,
                                  Boundary *b,
                                  int mode) {
    dst->clear ();
    dst->reserve (200);
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit)
    for (eit = b->edges_begin (cit); eit != b->edges_end (cit); ++eit) {
        intersect_info_t info;
        if (does_edge_intersect (&info, b, eit, ray_0, ray_dir)) {
            switch (mode) {
            /* FIXME broken
            case MODE_RAY_NEAREST:
                if (info.inc >= 0.) {
                    if (dst->size () && (*dst)[0].inc > info.inc)
                        dst->pop_back ();
                    dst->push_back (info);
                }
                break;
            */
            case MODE_RAY:
                if (info.inc >= 0.) {
            case MODE_LINE:
                    dst->push_back (info);
                }
                break;
            default:
                assert (0);
            }
        }
    }
}
                            
bool intersect_ray_boundary (intersect_info_t *dst,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *b) {
    assert (dst);
    intersect_buffer_t tmp;
    intersect_ray_boundary_impl (&tmp, ray_0, ray_dir, b, MODE_RAY_NEAREST);
    assert (tmp.size () == 0 || tmp.size () == 1);
    if (tmp.size ()) {
        *dst = tmp[0];
        return true;
    } else {
        return false;
    }
}

int intersect_ray_boundary  (intersect_buffer_t *dst,
                             const vec_t &ray_0, const vec_t &ray_dir,
                             Boundary *b) {
    assert (dst);
    intersect_ray_boundary_impl (dst, ray_0, ray_dir, b, MODE_RAY);
    return (int)dst->size ();
}

int intersect_line_boundary (intersect_buffer_t *dst,
                             const vec_t &line_0, const vec_t &line_dir,
                             Boundary *b)  {
    assert (dst);
    intersect_ray_boundary_impl (dst, line_0, line_dir, b, MODE_LINE);
    return (int)dst->size ();
}

bool intersect_vertex_rect (const vec_t &v, const rect_t &r) {
    if (r.left < v.x () && r.right > v.x ()) {
        if (r.top > v.y () && r.bottom < v.y ()) {
            return true;
        }
    }
    return false;
}
