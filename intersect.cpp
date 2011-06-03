// vim:et:sw=4:ts=4
// simple N^2 intersection algorithm
// it is tricky in a few points to make the algorithm
// return consistent results, i.e. that every closed contour
// has an even number of intersects with a test ray.
// this is very important for the introduce_divider_w0 function.
// to that end, some intersects are not returned.
//
// consider the signs of the signed distance between vertices
// and the test ray
//
// + + + + + +         0 intersects
// + + - + + +         2 intersects
// + 0 - - + +         2 intersects
// + 0 + + + +         0 intersects (!)
// + 0 - + + +         2 intersects
// + 0 0 0 0 +         0 intersects

#include "intersect.h"

namespace
{
    static
    double sgn (double x)
    {
        if (x > 0.)
            return 1.;
        else if (x < 0.)
            return -1.;
        else
            return x;
    }

    // this function does its best to return a sensible
    // intersection point between a line and a boundary
    // segment.  it does not decide whether there is an intersect,
    // that is decided from an exact predicate before this function is invoked.
    // (otherwise, inconsistent results may be generated)
    static
    void compute_intersection_point (
        intersect_info_t *dst, const Boundary &b,
        Boundary::edge_iterator eit, const vec_t &r0, const vec_t &dir_)
    {
        // solve r0 + lambda dir = mu v0 + (1-mu) v1
        // first, project on ortho(dir)
        const vec_t cdir_ = rot90_ccw (dir_);
        const double cn = dot (r0 - b.edge_vertex1 (eit), cdir_);
        const vec_t edge = b.edge_vertex0 (eit) - b.edge_vertex1 (eit);
        const double cn2 = dot (edge, cdir_);
        double mu = cn/cn2;
        mu = (mu < 1.) ? mu : 1.;  // clip to allowed range
        mu = (mu > 0.) ? mu : 0.;  // note that mu may be NaN before this

        // update inc parameter
        const double dn = mu * dot (edge, dir_)
            - dot (r0 - b.edge_vertex1 (eit), dir_);
        const double dn2 = dot (dir_, dir_);
        const double lambda = dn/dn2;
        const vec_t ivtx = r0 + lambda * dir_;
        const vec_t alternative_ivtx = mu * edge + b.edge_vertex1 (eit);
        dst->inc = lambda;
        dst->ivtx = ivtx;
        dst->iedge = eit;
        assert_not_nan (lambda);
        assert_not_nan (ivtx);

        const vec_t diff = ivtx - alternative_ivtx;
        if (dot (diff, diff) > 0.01*0.01 * dot (edge, edge))
        {
            std::cerr << "Inconsistency found while intersecting line/point\n";
            std::cerr << dot (diff, diff) << std::endl;
            std::abort ();
        }
    }

    // find the intersecting edges in the specified contour,
    // and store them in a STORE object.
    template <typename STORE>
    static
    unsigned find_intersecting_edges (
        STORE *st,
        const Boundary &b, Boundary::contour_iterator cit,
        const vec_t &r0, const vec_t &dir_)
    {
        unsigned ret = 0u;
        const vec_t cdir_ = rot90_ccw (dir_);
        const double offset = dot (cdir_, r0);

        Boundary::edge_iterator it = b.edges_begin (cit),
            end = b.edges_end (cit);
        int sign = sgn (dot (b.edge_vertex0 (it), cdir_) - offset);

        // make sure we don't start _on_ the line
        // zeroes are a problem since we don't know if we actually cross
        // the line here.
        if (sign == 0) {
            Boundary::edge_iterator a = it, aend = end;
            while (sign == 0) {
                ++a, ++aend;
                if (a == end)
                    // contour lies entirely on the ray, your call
                    return 0u;
                sign = sgn (dot (b.edge_vertex0 (a), cdir_) - offset);
            }
            it = a, end = aend;
        }

        // actually find the intersects
        for (; it != end; ++it) {
            int nsign = sgn (dot (b.edge_vertex1 (it), cdir_) - offset);
            // skip over on-the-line points until something interesting happens
            while (nsign == 0) {
                ++it;
                assert (it != end);
                nsign = sgn (dot (b.edge_vertex1 (it), cdir_) - offset);
            }
            if (nsign != sign) {
                intersect_info_t info;
                compute_intersection_point (&info, b, it, r0, dir_);
                info.sign = nsign;
                (*st)(info);
                sign = nsign;
            }
        }

        assert (ret % 2 == 0);
        return ret;
    }

    // collect all the intersects in the boundary
    template <typename STORE>
    static
    void find_intersections (
        STORE *st,
        const Boundary &b,
        const vec_t &r0, const vec_t &dir_)
    {
        assert_complete_boundary (b);
        Boundary::contour_iterator cit;
        for (cit = b.contours_begin (); cit != b.contours_end (); ++cit)
            find_intersecting_edges (st, b, cit, r0, dir_);
    }

    class IntersectCollector : public intersect_buffer_t
    {
    public:
        IntersectCollector (intersect_buffer_t *dst)
        {
            dst_ = dst;
            dst_->clear ();
            dst_->reserve (200u);
        }

        void operator() (const intersect_info_t &info)
        {
            dst_->push_back (info);
        }

    private:
        intersect_buffer_t *dst_;
    };

    class RayIntersectCollector : public IntersectCollector
    {
    public:
        RayIntersectCollector (intersect_buffer_t *dst)
            : IntersectCollector (dst)
        {
        }

        void operator() (const intersect_info_t &info)
        {
            if (info.inc >= 0.)
                static_cast <IntersectCollector &> (*this) (info);
        }
    };

    static
    void sort_intersections (intersect_buffer_t *dst)
    {
        std::sort (dst->begin (), dst->end (),
            intersect_info_t::by_normal_coordinate);
    }
}

unsigned
intersect_line_boundary (intersect_buffer_t *dst,
                         const vec_t &line_0, const vec_t &line_dir,
                         Boundary *b)  {
    IntersectCollector st (dst);
    find_intersections (&st, *b, line_0, line_dir);
    sort_intersections (dst);
    return (int)dst->size ();
}

unsigned
intersect_ray_boundary (intersect_buffer_t *dst,
                        const vec_t &line_0, const vec_t &line_dir,
                        Boundary *b)  {
    RayIntersectCollector st (dst);
    find_intersections (&st, *b, line_0, line_dir);
    sort_intersections (dst);
    return (int)dst->size ();
}



// implement assert_complete_boundary

void dump_intersect_buffer (std::ostream &os, const intersect_buffer_t &buf) {
    for (int i = 0; i != (int)buf.size (); ++i)
        os << buf[i].ivtx[0] << " " << buf[i].ivtx[1] << "\n";
}

bool intersect_vertex_rect (const vec_t &v, const rect_t &r) {
    if (r.left < v.x () && r.right > v.x ()) {
        if (r.top > v.y () && r.bottom < v.y ()) {
            return true;
        }
    }
    return false;
}
