// vim: et:sw=4:ts=4
// general utilities.
#ifndef UTIL_H_INCLUDED 
#define UTIL_H_INCLUDED 

#include <string>
#include <assert.h>
#include <limits.h>
#include "vector.h"
#include "matrix.h"
#include <algorithm>
#include <vector>
#include <ostream>

void die (const char *fmt, ...);

typedef Eigen::Vector2d vec_t;
typedef Eigen::Matrix2d mat_t;

// return eigenvalues of a
//  (a1 a2)
//  (b1 b2)    quadratic real matrix.
void eigensystem (struct EigenSystem *, double a1, double a2, double b1, double b2);
void eigensystem (struct EigenSystem *, const mat_t &);

struct EigenSystem {
    vec_t evec[2];
    double eval[2];

#ifndef NDEBUG
    void dump (std::ostream &os);
#endif
};

// return dyadic product  lhs (tensor) rhs in out,
// as a 2x2 matrix.
void dyadic_prod      (mat_t *out, const vec_t &lhs, const vec_t &rhs);
// dyadic product of vector with itself.
void dyadic_prod_self (mat_t *out, const vec_t &lhs);


//
// primitive Pixmap encap
//
class Pixmap {
public:
    Pixmap ();

    void resize (int dimx, int dimy);
    int size1 () const;
    int size2 () const;
    void init_zero ();

    typedef unsigned char val_t;
    static val_t min_val () { return 0; }
    static val_t max_val () { return UCHAR_MAX; }

    val_t &operator() (int x, int y);
    const val_t &operator() (int x, int y) const;

private:
    std::vector <val_t> my_data;
    int my_xdim, my_ydim;
};

void load_pgm (Pixmap *, const std::string &pgmfilename);
void invert (Pixmap *);


//
// data structure & toolbox for boundaries
//
class Boundary {
public:
    Boundary ();
    void clear ();

    enum { INVALID_EDGE = -1, INVALID_VERTEX = -2, INVALID_CONTOUR = -3 };

    typedef ::vec_t vec_t;
    int insert_vertex (const vec_t &);
    int insert_vertex (double x, double y);
    const vec_t &vertex (int) const;
    int num_vertices () const;

    struct edge_t {
        // neighbour edges (given anticlockwise for outer boundary)
        // these are indices into my_edge
        int prev, next;
        // vertices, indices into my_vert
        // vert0 == prev->vert1, vert1 == next->vert0
        int vert0, vert1;
        // contour identifier (arbitrary edge identifier on the contour)
        int contour;
        int label;
    };

    int insert_edge (int prev, int vert0, int vert1, int next = INVALID_EDGE);
    int split_edge (); // FIXME
    const edge_t &edge (int) const;
    int num_edges () const;

    // iterate over contours
    typedef std::vector <int>::const_iterator contour_iterator;
    contour_iterator contours_begin () const;
    contour_iterator contours_end () const;
    // iterator for traversing a contour
    // note that no "end" provided since contours
    // are assumed to be cycles
    // note that in case of an open (incomplete) contour, edges_begin ()
    // won't necessarily point at the first edge.
    class edge_iterator;
    edge_iterator    edges_begin (contour_iterator) const;
    edge_iterator    edges_end   (contour_iterator) const;

    class edge_iterator {
    public:
        edge_iterator ();
        edge_iterator &operator++ ();
        edge_iterator &operator-- ();
        const edge_t &operator* () const;
        const edge_t *operator-> () const;

        friend bool operator == (const edge_iterator &, const edge_iterator &);
        friend bool operator != (const edge_iterator &, const edge_iterator &);
    private:
        friend class Boundary;
        edge_iterator (const Boundary *, int, int);
        const Boundary *my_boundary;
        int my_position;
        int my_period;
    };

    double edge_length (edge_iterator it) const;
    // returns inflection angle between it and it+1,
    // in radians, in range [-pi/2, pi/2]
    // positive means convex edge
    double inflection_after_edge (edge_iterator it) const;
    double inflection_before_edge (edge_iterator it) const;
    // return outward-pointing normal vector
    vec_t edge_normal (edge_iterator) const;
    vec_t edge_tangent (edge_iterator) const;
    // return vertices
    vec_t edge_vertex0 (edge_iterator) const;
    vec_t edge_vertex1 (edge_iterator) const;
    // edge label
    int edge_label (edge_iterator) const;

private:
    vec_t &vertex (int i);
    edge_t &edge (int);
    void merge_contour_asc (int ed, int newc);
    void erase_contour (int c);
    edge_iterator edges_begin (int) const;
    edge_iterator edges_end (int) const;
    std::vector <vec_t> my_vert;
    std::vector <edge_t> my_edge;
    std::vector <int> my_contours;
};

void marching_squares (Boundary *, const Pixmap &);
void dump_contours (std::ostream &, const Boundary &);
void load_test_pixmap (Pixmap *);


//
// inline implementation
//
inline int Pixmap::size1 () const {
    return my_xdim;
}

inline int Pixmap::size2 () const {
    return my_ydim;
}

inline Pixmap::val_t &Pixmap::operator() (int x, int y) {
    const Pixmap *this_ = this;
    return const_cast <val_t &> ((*this_)(x, y));
}

inline const Pixmap::val_t &Pixmap::operator() (int x, int y) const {
#ifndef NDEBUG
    if (x < 0 || x >= my_xdim) {
        die ("invalid x coordinate %i [0:%i)", x, my_xdim);
    } else if (y < 0 || y >= my_ydim) {
        die ("invalid y coordinate %i [0:%i)", y, my_ydim);
    }
#endif
    return my_data[y * my_xdim + x];
}

inline const vec_t &Boundary::vertex (int i) const {
#ifndef NDEBUG
    return my_vert.at (i);
#else
    return my_vert[i];
#endif
}

inline int Boundary::insert_vertex (double x, double y) {
    return insert_vertex (vec_t (x, y));
}

inline vec_t &Boundary::vertex (int i) {
    const Boundary *this_ = this;
    return const_cast <vec_t &> (this_->vertex (i));
}

inline int Boundary::num_vertices () const {
    return (int)my_vert.size ();
}

inline const Boundary::edge_t &Boundary::edge (int i) const {
#ifndef NDEBUG
    return my_edge.at (i);
#else
    return my_edge[i];
#endif
}

inline Boundary::edge_t &Boundary::edge (int i) {
    const Boundary *this_ = this;
    return const_cast <Boundary::edge_t &> (this_->edge (i));
}

inline int Boundary::num_edges () const {
    return my_edge.size ();
}

inline Boundary::contour_iterator Boundary::contours_begin () const {
    return my_contours.begin ();
}

inline Boundary::contour_iterator Boundary::contours_end () const {
    return my_contours.end ();
}

inline Boundary::edge_iterator Boundary::edges_begin (Boundary::contour_iterator it) const {
    return edges_begin (*it);
}

inline Boundary::edge_iterator Boundary::edges_end (Boundary::contour_iterator it) const {
    return edges_end (*it);
}

inline Boundary::edge_iterator Boundary::edges_begin (int edge) const {
    return edge_iterator (this, edge, 0);
}

inline Boundary::edge_iterator Boundary::edges_end (int edge) const {
    return edge_iterator (this, edge, 1);
}

inline Boundary::edge_iterator::edge_iterator ()
    : my_position (Boundary::INVALID_EDGE), my_boundary (0), my_period (-1) { }

inline Boundary::edge_iterator::edge_iterator (const Boundary *b, int e, int p)
    : my_position (e), my_boundary (b), my_period (p) { }

inline const Boundary::edge_t &Boundary::edge_iterator::operator* () const {
    return my_boundary->edge(my_position);
}

inline const Boundary::edge_t *Boundary::edge_iterator::operator-> () const {
    return &(operator* ());
}

inline double Boundary::inflection_before_edge (Boundary::edge_iterator it) const {
    --it;
    return inflection_after_edge (it);
}

inline int Boundary::edge_label (Boundary::edge_iterator it) const {
    return it->label;
}

inline void eigensystem (EigenSystem *sys, const mat_t &mat) {
    eigensystem (sys, mat(0,0), mat(0,1), mat(1,0), mat(1,1));
}

inline void dyadic_prod (mat_t *mat, const vec_t &lhs, const vec_t &rhs) {
    (*mat)(0,0) = lhs[0] * rhs[0];
    (*mat)(0,1) = lhs[0] * rhs[1];
    (*mat)(1,0) = lhs[1] * rhs[0];
    (*mat)(1,1) = lhs[1] * rhs[1];
}

inline void dyadic_prod_self (mat_t *mat, const vec_t &lhs) {
    (*mat)(0,0) = lhs[0] * lhs[0];
    (*mat)(0,1) = 
    (*mat)(1,0) = lhs[0] * lhs[1];
    (*mat)(1,1) = lhs[1] * lhs[1];
}
#endif /* UTIL_H_INCLUDED */
