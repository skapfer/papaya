
#include "minkval.h"
#include <math.h>


static const double TWOPI = 2*M_PI;


// W100 = surface integral
//
// tensor is translation invariant.
//
// exact formulas for primitive bodies:
// * circle: 2 pi R
// * square: 4 a 
class W100 : public ScalarMinkowskiFunctional, public SurfaceIntegral {
public:
    W100 ()
        : ScalarMinkowskiFunctional ("W100") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            double &v = acc (b.edge_label (pos));
            v += b.edge_length (pos);
        }
    }
};

// W200 = surface integral weighted with curvature
// Inflection is calculated at end-of-edge, i.e. at vertex1 of the edge
// and counted towards the current label.
// This doesn't matter because we introduce inflection-zero vertices
// at label boundaries.
// @GERD Muessen wir das symmetrisieren?
class W200 : public ScalarMinkowskiFunctional, public SurfaceIntegral {
public:
    W200 ()
        : ScalarMinkowskiFunctional ("W200 / 2pi") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            double &v = acc (b.edge_label (pos));
            v += b.inflection_after_edge (pos) / TWOPI;
        }
    }
};

// W110 = surface integral weighted with location
//
// exact formulas for primitive bodies:
// * circle centered at origin: 0
// * square centered at origin: 0
class W110 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W110 ()
        : VectorMinkowskiFunctional ("W110") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            // center of gravity for edge
            vec_t avgvert = b.edge_vertex1 (pos);
            avgvert += b.edge_vertex0 (pos);
            avgvert /= 2;
            avgvert -= ref_vertex ();
            // weighted by edge length
            acc_ += avgvert * b.edge_length (pos);
        }
    }
};

// W101 = surface integral weighted with normal
// zero by definition for closed surfaces
class W101 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W101 ()
        : VectorMinkowskiFunctional ("W101") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            // normal of edge
            vec_t nrml = b.edge_normal (pos);
            // weighted by edge length
            acc_ += nrml * b.edge_length (pos);
        }
    }
};

// W210 = surface integral weighted with curvature and location
// Inflection is calculated at end-of-edge, i.e. at vertex1 of the edge
// and counted towards the current label.
// This doesn't matter because we introduce inflection-zero vertices
// at label boundaries.
//     label0    label1
// --X--------X---------X----
//            ]
// @GERD Muessen wir das symmetrisieren?
class W210 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W210 ()
        : VectorMinkowskiFunctional ("W210") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            vec_t vert = b.edge_vertex1 (pos);
            vert -= ref_vertex ();
            // weighted by edge length
            acc_ += vert * b.inflection_after_edge (pos);
        }
    }
};

// W201 = surface integral weighted with curvature and normal
// the only contributions come from vertices where the "label" changes,
// i.e. non-closed geometries
// these are neglected
class W201 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W201 ()
        : VectorMinkowskiFunctional ("W201") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            acc (b.edge_label (pos));
        }
    }
};

// W211 = surface integral weighted with curvature and location (tensor) normal
// calcuation is symmetrized to account for cases where the label changes.
//
// contribution per edge is:
// (vert1 - vert0) (tensor) (vert1 - vert0)  /  |vert1 - vert0|
//
// tensor is translation invariant.
//           symmetric.
//
// exact formulas for primitive bodies:
// * circle: R \pi IE
// * square: 2 a IE
//   IE being the 2x2 unit matrix.
class W211 : public MatrixMinkowskiFunctional, public SurfaceIntegral {
public:
    W211 ()
        : MatrixMinkowskiFunctional ("W211") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            mat_t incr;
            vec_t edgevec = b.edge_vertex1 (pos);
            edgevec -= b.edge_vertex0 (pos);
            dyadic_prod_self (&incr, edgevec);
            incr /= edgevec.norm ();
            acc (b.edge_label (pos)) += incr;
        }
    }
};

// W220 = surface integral weighted with curvature and location (tensor) location
// calcuation is symmetrized to account for cases where the label changes.
//
// contribution per edge is:
// 1/2 sum_i inflection_i * vertex_i (tensor) vertex_i
//
// tensor is symmetric.
//
// exact formulas for primitive bodies:
// * circle centered at origin:    \pi R^2 IE       << DISAGREES +20% >>
// * square in positive quadrant:  \pi a^2 IE  
//   square centered at origin:    \pi/2 a^2 IE                   -1%
//   IE being the 2x2 unit matrix.
class W220 : public MatrixMinkowskiFunctional, public SurfaceIntegral {
public:
    W220 ()
        : MatrixMinkowskiFunctional ("W220") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        while (pos != end) {
            mat_t incr;
            vec_t loc = b.edge_vertex1 (pos);
            loc -= ref_vertex ();
            dyadic_prod_self (&incr, loc);
            incr *= .5 * b.inflection_after_edge (pos);
            acc (b.edge_label (pos)) += incr;
            ++pos;
            acc (b.edge_label (pos)) += incr;
        }
    }
};

// centre of mass of vertices, i.e. average over all vertices
// (for testing)
static vec_t centre_of_mass (const Boundary &b) {
    vec_t acc = vec_t (0., 0.);
    int ctr = 0;
    Boundary::edge_iterator eit;
    Boundary::contour_iterator cit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        for (eit = b.edges_begin (cit); eit != b.edges_end (cit); ++eit) {
            acc += b.edge_vertex0 (eit);
            ++ctr;
        }
    }
    return acc / ctr;
}

// driver function for testing the functionals
void calculate_all_surface_integrals (const Boundary &b) {
    std::vector <SurfaceIntegral *> sfints;
    std::vector <SurfaceIntegral *>::iterator iit;
    Boundary::contour_iterator cit;

    W100 w100;
    sfints.push_back (&w100);
    W200 w200;
    sfints.push_back (&w200);
    W110 w110;
    sfints.push_back (&w110);
    W101 w101;
    sfints.push_back (&w101);
    W201 w201;
    sfints.push_back (&w201);
    W211 w211;
    sfints.push_back (&w211);
    W210 w210;
    sfints.push_back (&w210);
    W220 w220;
    sfints.push_back (&w220);

    {
        // set ref. vertex to center-of-mass
        vec_t refvert = centre_of_mass (b);
        for (iit = sfints.begin (); iit != sfints.end (); ++iit)
            (*iit)->ref_vertex (refvert);
    }

    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        for (iit = sfints.begin (); iit != sfints.end (); ++iit)
            (*iit)->add_contour (b, b.edges_begin (cit), b.edges_end (cit));
    }

    w100.dump (std::cout);
    w200.dump (std::cout);
    w110.dump (std::cout);
    w101.dump (std::cout);
    w210.dump (std::cout);
    w201.dump (std::cout);
    w211.dump (std::cout);
    w220.dump (std::cout);
}

