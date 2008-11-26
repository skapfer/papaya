#ifndef PAPAYA_LINT_H_INCLUDED 
#define PAPAYA_LINT_H_INCLUDED 

// lint.h contains visitors which prepare the body data for the 
// actual calculation, doing sanity checks and precomputing stuff.

#include "FacetBoundedBody.h"
#ifndef DIMENSION
#error "Need DIMENSION to be defined."
#endif
#include <map>
#include "geo.h"

namespace Papaya2 {

// This visitor is used to recover neighbourship information from data
// lacking it.  However, it is not always possible to do so for certain
// singular bodies (see WARNING below).  Avoid if possible.
// This can be run on a subset of the data safely.
struct AssignNeighboursVisitor {
    template <typename BODY>
    void operator() (BODY &body, int f) {
#if (DIMENSION==3)
        for (int i = 0; i != 3; ++i) {
            int v0 = body.facet_vertex_id (f, i);
            int v1 = body.facet_vertex_id (f, (i+1) % 3);
            int nf = request_unpaired_edge (v1, v0);
            if (nf != INVALID_FACET) {
                pair_edges (body, f, i, v0, v1, nf);
            } else {
                announce_unpaired_edge (f, v0, v1);
            }
            body.primar_[f].neigh[i] = nf;
        }
#else
        #error "NYI!"
#endif
    }

private:
    struct ann_edge_t {
        int v1;
        int f;
        bool used;
    };

    template <typename BODY>
    void pair_edges (BODY &body, int f1, int i, int v0, int v1, int f2) {
        assert (body.facet_vertex_id (f1, i) == v0);
        int j = body_facet_which_edge (body, f2, v1, v0);
        assert (j >= 0);
        assert (body.primar_[f1].neigh[i] == INVALID_FACET);
        assert (body.primar_[f2].neigh[j] == INVALID_FACET);
        body.primar_[f1].neigh[i] = f2;
        body.primar_[f2].neigh[j] = f1;
    }

    // this is a waste of memory, but the root of all evil theorem
    // applies...
    std::map <int, std::vector <ann_edge_t> > ann_e_;

    void announce_unpaired_edge (int f, int v0, int v1) {
        ann_edge_t e;
        e.v1 = v1;
        e.f = f;
        e.used = false;
        ann_e_[v0].push_back (e);
    }

    int request_unpaired_edge (int v0, int v1) {
        std::vector <ann_edge_t> &es = ann_e_[v0];
        for (int i = 0; i != (int)es.size (); ++i) {
            if (es[i].v1 == v1) {
                if (es[i].used++)
                    std::cerr << "WARNING:  edge " << v0 << "--"
                              << v1 << "found more than once.\n";
                else
                    return es[i].f;
            }
        }
        return INVALID_FACET;
    }

public:
    ~AssignNeighboursVisitor () {}
};

// PrecomputeVisitor computes some of the precomputed quantities in the
// FacetDerivedData struct.
struct PrecomputeVisitor {
    template <typename BODY>
    void operator() (BODY &body, int facet) const {
        // check that nothing is NaN
        for (int i = 0; i != DIMENSION; ++i) {
            typename BODY::vec_t v = body.facet_vertex (facet, i);
            if (!is_finite (v)) {
                std::cerr << "ERROR:  Broken input data at vertex "
                          << body.facet_vertex_id (facet, i) << ": "
                          << v << "\n";
                abort ();
            }
        }

#if (DIMENSION==3)
        // check that the vertex id's are different
        if (body.facet_vertex_id (facet, 0) == body.facet_vertex_id (facet, 1) ||
            body.facet_vertex_id (facet, 1) == body.facet_vertex_id (facet, 2) ||
            body.facet_vertex_id (facet, 0) == body.facet_vertex_id (facet, 2)) {
            std::cerr << "ERROR:  triangle carries at least two identical vertices:\n"
                      << body.facet_vertex_id (facet, 0) << "\n"
                      << body.facet_vertex_id (facet, 1) << "\n"
                      << body.facet_vertex_id (facet, 2) << "\n";
            abort ();
        }

        // check the proper assignment of neighbour links
        for (int i = 0; i != 3; ++i) {
            int nf = body.primar_[facet].neigh[i];
            if (nf == INVALID_FACET)
                continue;
            int v0 = body.facet_vertex_id (facet, i);
            int v1 = body.facet_vertex_id (facet, (i+1) % 3);
            int oe = body_facet_which_edge (body, nf, v1, v0);
            if (oe < 0 || body.primar_[nf].neigh[oe] != facet) {
                std::cerr << "ERROR:  neighbourship information is broken.  Aborting.\n";
                abort ();
            }
        }

        const vec3_t &v0 = body.facet_vertex (facet, -1);
        const vec3_t &v1 = body.facet_vertex (facet, 1);
        const vec3_t &v2 = body.facet_vertex (facet, 2);
        vec3_t n = PapayaGeometry::plane_normal (v0, v1, v2);
        double v = PapayaGeometry::triangle_area (v0, v1, v2);
        body.deriv_[facet].normal = n;
        assert_not_nan (n);
        body.deriv_[facet].volume = v;
        assert_not_nan (v);
        vec3_t c = PapayaGeometry::average (v0, v1, v2);
        body.deriv_[facet].center = c;
        assert_not_nan (c);
        body.deriv_[facet].flag = 42; // dummy to find uninitialised stuff
#else
        #error "NYI!"
#endif
    }
};

// Precompute2Visitor computes the rest of the precomputed quantities in the
// FacetDerivedData struct.
// It assumes that PrecomputeVisitor has run for all the facets
// (because it needs normals) and that neighbourship information is available.
struct Precompute2Visitor {
    static bool triple_ineq (double x, double y, double z) {
        return x < y and y < z;
    }

    // this function fails if the two facets aren't neighbours
    // it returns some nonsense then.
    // otherwise, it returns one edge in facet0 shared with facet1.
    // (there may be multiple shared edges in case of a seriously
    // bodged body)
    // FIXME we might be able to optimize this if the neighbours
    // are stored in the proper order.
    template <typename BODY>
    static
    vec3_t find_shared_edge (const BODY &body,
                             int facet0, int facet1) {
        int v0 = body.facet_vertex_id (facet1, 0);
        int v1 = body.facet_vertex_id (facet1, 1);
        int v2 = body.facet_vertex_id (facet1, 2);
        vec3_t ret;
        int vX = body.facet_vertex_id (facet0, 0);
        if (vX != v0 && vX != v1 && vX != v2) {
            ret  = body.facet_vertex (facet0, 2);
            ret -= body.facet_vertex (facet0, 1);
        } else {
            vX = body.facet_vertex_id (facet0, 1);
            if (vX != v0 && vX != v1 && vX != v2) {
                ret  = body.facet_vertex (facet0, 0);
                ret -= body.facet_vertex (facet0, 2);
            } else {
                ret  = body.facet_vertex (facet0, 1);
                ret -= body.facet_vertex (facet0, 0);
            }
        }
        return ret;
    }

    template <typename BODY>
    void operator() (BODY &body, int f) const {
#if (DIMENSION==3)
        // pre-calculate mean curvature contribution for each edge of this
        // facet.
        for (int i = 0; i != 3; ++i) {
            int n = body.facet_neighbour_facet (f, i);
            double mc = NAN;
            if (n != INVALID_FACET) {
                vec3_t n1xn2 = cross (body.facet_normal (f), body.facet_normal (n));
                vec3_t edge  = find_shared_edge (body, f, n);
                double edge_len = edge.norm ();
                double dp = dot (n1xn2, edge);
                dp /= edge_len;
                // these should be parallel or antiparallel; if not,
                // we have messed up.
                assert (dp > .95 || dp < -.95);
                mc = asin (dp);
                assert (dp * mc >= 0);
                if (! triple_ineq (-.8, mc/M_PI, .8)) {
                    std::cerr << "WARNING: very large inflection between facets "
                              << f << " and " << n << " (" << mc/M_PI << "*pi)\n";
                }
                mc *= edge_len;
                assert (!isnan (mc));
            }
            body.deriv_[f].meanc[i] = mc;
        }
        // FIXME does not compute gaussian curvature.
        body.deriv_[f].gaussian[0] = NAN;
        body.deriv_[f].gaussian[1] = NAN;
        body.deriv_[f].gaussian[2] = NAN;
#else
        #error "NYI!"
#endif
    }
};

void FacetBoundedBody::postcompute_neighbourships () {
    AssignNeighboursVisitor vis;
    this->visit_each_facet (vis);
}

void FacetBoundedBody::body_complete () {
    PrecomputeVisitor vis;
    Precompute2Visitor vis2;
    this->visit_each_facet (vis);
    this->visit_each_facet (vis2);
}

}

#endif // PAPAYA_LINT_H_INCLUDED
