#ifndef PAPAYA_LINT_H_INCLUDED 
#define PAPAYA_LINT_H_INCLUDED 

// lint.h contains visitors which prepare the body data for the 
// actual calculation, doing sanity checks and precomputing stuff.

#include "FacetBoundedBody.h"
#ifndef DIMENSION
#error "Need DIMENSION to be defined."
#endif
#include "geo.h"

namespace Papaya2 {

struct InspectionVisitor {
    template <typename BODY>
    void operator() (const BODY &body, int facet) const {
        // FIXME
    }
};

// PrecomputeVisitor computes some of the precomputed quantities in the
// FacetDerivedData struct.
struct PrecomputeVisitor {
    template <typename BODY>
    void operator() (BODY &body, int facet) const {
#if (DIMENSION==3)
        const vec3_t &v0 = body.facet_vertex (facet, 0);
        const vec3_t &v1 = body.facet_vertex (facet, 1);
        const vec3_t &v2 = body.facet_vertex (facet, 2);
        vec3_t n = PapayaGeometry::plane_normal (v0, v1, v2);
        double v = PapayaGeometry::triangle_area (v0, v1, v2);
        body.deriv_[facet].normal = n;
        body.deriv_[facet].volume = v;
        vec3_t c = PapayaGeometry::average (v0, v1, v2);
        body.deriv_[facet].center = c;
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
    bool triple_ineq (double x, double y, double z) {
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
    vec3_t find_shared_edge (const BODY &body,
                             int facet0, int facet1) {
        int v0 = body.facet_vertex_id (facet1, 0);
        int v1 = body.facet_vertex_id (facet1, 1);
        int v2 = body.facet_vertex_id (facet1, 2);
        int ret1, ret0;
        int vX = body.facet_vertex_id (facet0, 0);
        if (vX != v0 && vX != v1 && vX != v2) {
            ret0 = body.facet_vertex_id (facet0, 1);
            ret1 = body.facet_vertex_id (facet0, 2);
        } else {
            vX = body.facet_vertex_id (facet0, 1);
            if (vX != v0 && vX != v1 && vX != v2) {
                ret0 = body.facet_vertex_id (facet0, 2);
                ret1 = body.facet_vertex_id (facet0, 0);
            } else {
                ret0 = body.facet_vertex_id (facet0, 0);
                ret1 = body.facet_vertex_id (facet0, 1);
            }
        }
        return body.vertex_coords (ret1) - body.vertex_coords (ret0);
    }

    template <typename BODY>
    void operator() (BODY &body, int f) const {
#if (DIMENSION==3)
        for (int i = 0; i != 3; ++i) {
            int n = body.facet_neighbour_facet (f, i);
            double inf = NAN;
            if (n != INVALID_FACET) {
                vec3_t n1xn2 = cross (body.facet_normal (f), body.facet_normal (n));
                vec3_t edge  = find_shared_edge (body, f, n);
                edge.normalize ();
                double dp = dot (n1xn2, edge);
                // these should be parallel or antiparallel; if not,
                // we have messed up.
                assert (dp > .95 || dp < -.95);
                inf = asin (dp);
                assert (dp * inf >= 0);
                if (! triple_ineq (-.8, inf/M_PI, .8)) {
                    std::cerr << "WARNING: very large inflection between facets "
                              << f << " and " << n << " (" << inf/M_PI << "*pi)\n";
                }
            }
            body.deriv_[f].meanc[i] = inf;
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

void FacetBoundedBody::body_complete () {
    this->visit_each_facet (PrecomputeVisitor ());
}

}

#endif // PAPAYA_LINT_H_INCLUDED
