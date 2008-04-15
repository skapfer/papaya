
#include "minkval.h"
#include <math.h>


static const double TWOPI = 2*M_PI;


// W100 = surface integral
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




void calculate_all_surface_integrals (const Boundary &b) {
    std::vector <SurfaceIntegral *> sfints;

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
    W210 w210;
    sfints.push_back (&w210);

    Boundary::contour_iterator cit;
    std::vector <SurfaceIntegral *>::iterator iit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        for (iit = sfints.begin (); iit != sfints.end (); ++iit)
            (*iit)->add_contour (b, b.edges_begin (cit), b.edges_end (cit));
    }

    w100.dump (std::cout);
    w200.dump (std::cout);
    w110.dump (std::cout);
    w101.dump (std::cout);
    w201.dump (std::cout);
    w210.dump (std::cout);
}

