
#include "minkval.h"
#include <math.h>


static const double TWOPI = 2*M_PI;
static const double W0_NORMALIZATION = 1.;
static const double W1_NORMALIZATION = 1.;
static const double W2_NORMALIZATION = 1.;


// W000 = volume integral
// converted to a surface integral by GaussGreen theorem.
//
// tensor is motion invariant.
//           rank-0.
//           hom. degree 2.
//
// exact formulas for primitive bodies:
// * ellipsis:  pi a b
// * rectangle: a b
class W000 : public ScalarMinkowskiFunctional, public SurfaceIntegral {
public:
    W000 ()
        : ScalarMinkowskiFunctional ("W000") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t edge_grav = b.edge_vertex0 (pos);
            edge_grav += b.edge_vertex1 (pos);
            edge_grav *= .25 * b.edge_length (pos);
            double sc = dot (b.edge_normal (pos), edge_grav);
            acc (b.edge_label (pos)) += sc;
        }
    }
};

// W100 = surface integral
//
// tensor is motion invariant.
//           rank-0.
//           hom. degree 1.
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
//
// tensor is motion invariant.
//           rank-0.
//           hom. degree 0.
//
// exact formulas for primitive bodies:
// * circle: 2 pi
// * square: 2 pi
class W200 : public ScalarMinkowskiFunctional, public SurfaceIntegral {
public:
    W200 ()
        : ScalarMinkowskiFunctional ("W200") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            double &v = acc (b.edge_label (pos));
            v += b.inflection_after_edge (pos);
        }
    }
};

// W010 = volume integral weighted with location
// converted to a surface integral by GaussGreen theorem.
//
// tensor is motion covariant.
//           rank-1.
//           hom. degree 3.
//
// volume center of gravity times total volume.
//
// exact formulas for primitive bodies:
// * circle centered at origin: 0
// * square centered at origin: 0
class W010 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W010 ()
        : VectorMinkowskiFunctional ("W010") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            const vec_t &v1 = b.edge_vertex1 (pos);
            const vec_t &v0 = b.edge_vertex0 (pos);
            vec_t first_factor = v1;
            first_factor -= v0;
            first_factor /= 6.;
            vec_t second_factor;
            second_factor[0] = v1[0]*v1[0] + v0[0]*v0[0] + v0[0]*v1[0];
            second_factor[1] = v1[1]*v1[1] + v0[1]*v0[1] + v0[1]*v1[1];
            vec_t &acc_ = acc (b.edge_label (pos));
            acc_[0] +=  first_factor[1] * second_factor[0];
            acc_[1] += -first_factor[0] * second_factor[1];
        }
    }
};

// W110 = surface integral weighted with location
//
// tensor is motion covariant.
//           rank-1.
//           hom. degree 2.
//
// surface center of gravity times total surface.
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

// W210 = surface integral weighted with curvature and location
// Inflection is calculated at end-of-edge, i.e. at vertex1 of the edge
// and counted towards the current label.
// This doesn't matter because we introduce inflection-zero vertices
// at label boundaries.
//     label0    label1
// --X--------X---------X----
//            ]
// @GERD Muessen wir das symmetrisieren?
//
// tensor is motion covariant.
//           rank-1.
//           hom. degree 1.
//
// surface center of gravity times total curvature
class W210 : public VectorMinkowskiFunctional, public SurfaceIntegral {
public:
    W210 ()
        : VectorMinkowskiFunctional ("W210") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            vec_t vert = b.edge_vertex1 (pos);
            vert -= ref_vertex ();
            acc_ += vert * b.inflection_after_edge (pos);
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
//           rank-2.
//           hom. degree 1.
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
//           rank-2.
//           hom. degree 2.
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


// W120 = surface integral weighted with location (tensor) location
//
// contribution per edge is:
// 1/3 sum_i vertex_i (tensor) vertex_i + vertex_i+1 (tensor) vertex_(i+1)
//           + 1/2 vertex_i (tensor) vertex_i+1 + 1/2 vertex_i+1 (tensor) vertex_i
//
// tensor is symmetric.
//           rank-2.
//           hom. degree 3.
//
// exact formulas for primitive bodies:
// * circle centered at origin:    \pi R^3 IE     
// * square in positive quadrant: 
//   square centered at origin:  
//   IE being the 2x2 unit matrix.
class W120 : public MatrixMinkowskiFunctional, public SurfaceIntegral {
public:
    W120 ()
        : MatrixMinkowskiFunctional ("W120") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = W1_NORMALIZATION / 3.;
        for (; pos != end; ++pos) {
            mat_t &acc_ = acc (b.edge_label (pos));
            mat_t incr;
            double l_prefactor = prefactor * b.edge_length (pos);
            // vertex 1
            vec_t loc1 = b.edge_vertex1 (pos);
            loc1 -= ref_vertex ();
            dyadic_prod_self (&incr, loc1);
            incr *= l_prefactor;
            acc_ += incr;
            // vertex 0
            vec_t loc0 = b.edge_vertex0 (pos);
            loc0 -= ref_vertex ();
            dyadic_prod_self (&incr, loc0);
            incr *= l_prefactor;
            acc_ += incr;
            // mixed term
            dyadic_prod_symmetrized (&incr, loc0, loc1);
            incr *= l_prefactor;
            acc_ += incr;
        }
    }
};

// W020 = volume integral weighted with location (tensor) location
// converted to a surface integral by GaussGreen theorem.
//
// contribution per edge is:
//
// tensor is symmetric.
//           rank-2.
//           hom. degree 4.
//
// exact formulas for primitive bodies:
// * circle centered at origin:       .25 * R^4 * pi * IE
// * square in positive quadrant: 
//   square centered at origin:  
//   IE being the 2x2 unit matrix.
class W020 : public MatrixMinkowskiFunctional, public SurfaceIntegral {
public:
    W020 ()
        : MatrixMinkowskiFunctional ("W020") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = W0_NORMALIZATION;
        for (; pos != end; ++pos) {
            const vec_t &v1 = b.edge_vertex1 (pos);
            const vec_t &v0 = b.edge_vertex0 (pos);
            mat_t &acc_ = acc (b.edge_label (pos));
            // xx element
            double prefactor = W0_NORMALIZATION / 12. * (v1[1] - v0[1]);
            acc_(0,0) += prefactor * (v0[0] + v1[0]) * (v0[0]*v0[0] + v1[0]*v1[0]);
            // xy element
            double t = v0[0]*v0[0] * (3.*v0[1] + v1[1]);
            t += v1[0]*v1[0] * (3.*v1[1] + v0[1]);
            t *= .5;
            t += (v0[0] + v1[0]) * (v0[1] + v1[1]);
            acc_(0,1) += prefactor * t;
            acc_(1,0) = acc_(0,1);
            // yy element
            prefactor = W0_NORMALIZATION / 12. * (v0[0] - v1[0]);
            acc_(1,1) += prefactor * (v0[1] + v1[1]) * (v0[1]*v0[1] + v1[1]*v1[1]);
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

    W000 w000;
    sfints.push_back (&w000);
    W100 w100;
    sfints.push_back (&w100);
    W200 w200;
    sfints.push_back (&w200);
    W010 w010;
    sfints.push_back (&w010);
    W110 w110;
    sfints.push_back (&w110);
    W211 w211;
    sfints.push_back (&w211);
    W210 w210;
    sfints.push_back (&w210);
    W220 w220;
    sfints.push_back (&w220);
    W120 w120;
    sfints.push_back (&w120);
    W020 w020;
    sfints.push_back (&w020);

    {
        // set ref. vertex to center-of-mass
        vec_t refvert = centre_of_mass (b);
        std::cout << "ref. vertex: " << refvert << "\n";
        for (iit = sfints.begin (); iit != sfints.end (); ++iit)
            (*iit)->ref_vertex (refvert);
    }

    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        for (iit = sfints.begin (); iit != sfints.end (); ++iit)
            (*iit)->add_contour (b, b.edges_begin (cit), b.edges_end (cit));
    }

    std::cout << "]]]]] SCALARS ]]]]]]]]]]]]]]]]]]]]]]\n";
    w000.dump (std::cout);
    w100.dump (std::cout);
    w200.dump (std::cout);
    std::cout << "]]]]] VECTORS ]]]]]]]]]]]]]]]]]]]]]]\n";
    w010.dump (std::cout);
    w110.dump (std::cout);
    w210.dump (std::cout);
    std::cout << "]]]]] TENSORS ]]]]]]]]]]]]]]]]]]]]]]\n";
    std::cout << "Hom. degree 1\n";
    w211.dump (std::cout);
    std::cout << "Hom. degree 2\n";
    w220.dump (std::cout);
    std::cout << "Hom. degree 3\n";
    w120.dump (std::cout);
    std::cout << "Hom. degree 4\n";
    w020.dump (std::cout);
}

