
#include "minkval.h"
#include <math.h>


// change the 1 to a 0 in order to use the normalization
// described in the paper
#if 1
// Breidenbach normalization
static const double W0_NORMALIZATION = 1.;
static const double W1_NORMALIZATION = 1.;
static const double W2_NORMALIZATION = 1.;
#else
static const double W0_NORMALIZATION = 1.;
static const double W1_NORMALIZATION = .5;
static const double W2_NORMALIZATION = .5;
#endif


void print_version_header (std::ostream &of) {
    of << "# papaya version " << VERSION <<
        "\n# (c) 2008-2010 Sebastian Kapfer "
        "<sebastian.kapfer@physik.uni-erlangen.de>\n";
}

// W000 = volume integral
// converted to a surface integral by GaussGreen theorem.
//
// tensor is motion invariant.
//           rank-0.
//           hom. degree 2.
//           linear dependancy with W111.
//
// exact formulas for primitive bodies:
// * ellipsis:  pi a b
// * rectangle: a b
class W000 : public ScalarMinkowskiFunctional {
public:
    W000 ()
        : ScalarMinkowskiFunctional ("W000") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            vec_t edge_grav = b.edge_vertex0 (pos);
            edge_grav += b.edge_vertex1 (pos);
            edge_grav *= .25 * b.edge_length (pos);
            double sc = dot (b.edge_normal (pos), edge_grav);
            sc *= W0_NORMALIZATION;
            acc (b.edge_label (pos)) += sc;
        }
    }
};

// W100 = surface integral
//
// tensor is motion invariant.
//           rank-0.
//           hom. degree 1.
//           together with W211, determines W102.
//
// exact formulas for primitive bodies:
// * circle: 2 pi R
// * square: 4 a 
//   rectangle: Lx Ly
class W100 : public ScalarMinkowskiFunctional {
public:
    W100 ()
        : ScalarMinkowskiFunctional ("W100") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            double &v = acc (b.edge_label (pos));
            v += b.edge_length (pos) * W1_NORMALIZATION;
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
//           linear dependancy with W202.
//
// exact formulas for primitive bodies:
// * circle: 2 pi
// * square: 2 pi
class W200 : public ScalarMinkowskiFunctional {
public:
    W200 ()
        : ScalarMinkowskiFunctional ("W200") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            double &v = acc (b.edge_label (pos));
            v += b.inflection_after_edge (pos) * W2_NORMALIZATION;
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
class W010 : public VectorMinkowskiFunctional {
public:
    W010 ()
        : VectorMinkowskiFunctional ("W010") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            int l = b.edge_label (pos);
            const vec_t &v1 = b.edge_vertex1 (pos) - ref_vertex (l);
            const vec_t &v0 = b.edge_vertex0 (pos) - ref_vertex (l);
            vec_t first_factor = v1;
            first_factor -= v0;
            first_factor *= W0_NORMALIZATION / 6.;
            vec_t second_factor;
            second_factor[0] = v1[0]*v1[0] + v0[0]*v0[0] + v0[0]*v1[0];
            second_factor[1] = v1[1]*v1[1] + v0[1]*v0[1] + v0[1]*v1[1];
            vec_t &acc_ = acc (l);
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
class W110 : public VectorMinkowskiFunctional {
public:
    W110 ()
        : VectorMinkowskiFunctional ("W110") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = W1_NORMALIZATION;
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            // center of gravity for edge
            vec_t avgvert = b.edge_vertex1 (pos);
            avgvert += b.edge_vertex0 (pos);
            avgvert /= 2;
            avgvert -= ref_vertex (b.edge_label (pos));
            // weighted by edge length
            acc_ += (prefactor * b.edge_length (pos)) * avgvert;
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
class W210 : public VectorMinkowskiFunctional {
public:
    W210 ()
        : VectorMinkowskiFunctional ("W210") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = W2_NORMALIZATION;
        for (; pos != end; ++pos) {
            vec_t &acc_ = acc (b.edge_label (pos));
            vec_t vert = b.edge_vertex1 (pos);
            vert -= ref_vertex (b.edge_label (pos));
            acc_ += (prefactor * b.inflection_after_edge (pos)) * vert;
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
//           together with W100, determines W102.
//
// exact formulas for primitive bodies:
// * circle: R \pi IE
// * square: 2 a IE
//   IE being the 2x2 unit matrix.
class W211 : public MatrixMinkowskiFunctional {
public:
    W211 ()
        : MatrixMinkowskiFunctional ("W211") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            mat_t incr;
            vec_t edgevec = b.edge_vertex1 (pos);
            edgevec -= b.edge_vertex0 (pos);
            dyadic_prod_self (&incr, edgevec);
            incr *= W2_NORMALIZATION / edgevec.norm ();
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
//   rectangle centered at origin  \pi/2 diag (Lx^2, Ly^2)
//   IE being the 2x2 unit matrix.
class W220 : public MatrixMinkowskiFunctional {
public:
    W220 ()
        : MatrixMinkowskiFunctional ("W220") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = .5 * W2_NORMALIZATION;
        while (pos != end) {
            mat_t incr;
            vec_t loc = b.edge_vertex1 (pos);
            loc -= ref_vertex (b.edge_label (pos));
            dyadic_prod_self (&incr, loc);
            incr *= prefactor * b.inflection_after_edge (pos);
            acc (b.edge_label (pos)) += incr;
            loc = b.edge_vertex0 (pos);
            loc -= ref_vertex (b.edge_label (pos));
            dyadic_prod_self (&incr, loc);
            incr *= prefactor * b.inflection_before_edge (pos);
            acc (b.edge_label (pos)) += incr;
            ++pos;
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
class W120 : public MatrixMinkowskiFunctional {
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
            loc1 -= ref_vertex (b.edge_label (pos));
            dyadic_prod_self (&incr, loc1);
            incr *= l_prefactor;
            acc_ += incr;
            // vertex 0
            vec_t loc0 = b.edge_vertex0 (pos);
            loc0 -= ref_vertex (b.edge_label (pos));
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

// W102 = surface integral weighted with normal (tensor) normal
//
// contribution per edge is:
// |vert1 - vert0| * normal (tensor) normal
//
// tensor is symmetric.
//           translation invariant.
//           rank-2.
//           hom. degree 1.
//
// exact formulas for primitive bodies:
class W102 : public MatrixMinkowskiFunctional {
public:
    W102 ()
        : MatrixMinkowskiFunctional ("W102") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        const double prefactor = W1_NORMALIZATION;
        for (; pos != end; ++pos) {
            mat_t &acc_ = acc (b.edge_label (pos));
            assert_not_nan (acc_);
            mat_t incr;
            dyadic_prod_self (&incr, b.edge_normal (pos));
            incr *= b.edge_length (pos) * prefactor;
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
class W020 : public MatrixMinkowskiFunctional {
public:
    W020 ()
        : MatrixMinkowskiFunctional ("W020") { }

    virtual void add_contour (const Boundary &b, edge_iterator pos, edge_iterator end) {
        for (; pos != end; ++pos) {
            int l = b.edge_label (pos);
            const vec_t &v1 = b.edge_vertex1 (pos) - ref_vertex (l);
            const vec_t &v0 = b.edge_vertex0 (pos) - ref_vertex (l);
            mat_t &acc_ = acc (l);
            assert_not_nan (acc_);
            // xx element
            double prefactor = W0_NORMALIZATION / 12. * (v1[1] - v0[1]);
            acc_(0,0) += prefactor * (v0[0] + v1[0]) * (v0[0]*v0[0] + v1[0]*v1[0]);
            // xy element
            double t = v0[0]*v0[0] * (3.*v0[1] + v1[1]);
            t += v1[0]*v1[0] * (3.*v1[1] + v0[1]);
            t *= .5;
            t +=  v0[0] * v1[0]  * (v0[1] + v1[1]);
            acc_(0,1) += prefactor * t;
            acc_(1,0) = acc_(0,1);
            // yy element
            prefactor = W0_NORMALIZATION / 12. * (v0[0] - v1[0]);
            acc_(1,1) += prefactor * (v0[1] + v1[1]) * (v0[1]*v0[1] + v1[1]*v1[1]);
        }
    }
};


template <> mat_t MatrixMinkowskiFunctional::dummy_acc  = mat_t ();
template <> vec_t  VectorMinkowskiFunctional::dummy_acc = vec_t ();
template <> double ScalarMinkowskiFunctional::dummy_acc = double ();

ScalarMinkowskiFunctional *create_w000 () { return new W000; }
ScalarMinkowskiFunctional *create_w100 () { return new W100; }
ScalarMinkowskiFunctional *create_w200 () { return new W200; }
VectorMinkowskiFunctional *create_w010 () { return new W010; }
VectorMinkowskiFunctional *create_w110 () { return new W110; }
VectorMinkowskiFunctional *create_w210 () { return new W210; }
MatrixMinkowskiFunctional *create_w020 () { return new W020; }
MatrixMinkowskiFunctional *create_w120 () { return new W120; }
MatrixMinkowskiFunctional *create_w102 () { return new W102; }
MatrixMinkowskiFunctional *create_w220 () { return new W220; }
MatrixMinkowskiFunctional *create_w211 () { return new W211; }
