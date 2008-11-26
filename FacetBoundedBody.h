#ifndef FACETBOUNDEDBODY 
#define FACETBOUNDEDBODY 

#ifndef PAPAYA2
#error "This header is PAPAYA2 code and cannot be mixed with PAPAYA1."
#endif // PAPAYA2
#include "util.h"


namespace Papaya2 {


enum { INVALID_FACET = -1, INVALID_DOMAIN = -5 };

// data as read in
struct FacetPrimaryData {
    int vertex[DIMENSION];
    int neigh[DIMENSION];
    int label;
    int domain;
};

/* temporary data which is built by PAPAYA in order to compute Minkowski
 * valuations.  This is used by the current implementation of FacetBoundedBody.
 * This class does not initialize any of its data members.
 */
struct FacetDerivedData {
    vec_t normal;
    double volume;
#if (DIMENSION == 3)
    double gaussian[3];
    double meanc[3];
#elif (DIMENSION == 2)
    double inflect[2];
#endif
    vec_t center;
    int flag; // internal flag for graph code
};

/* This class is _one_ implementation of the data storage back-end
 * for PAPAYA.  Other implementations can be written (e.g. for
 * parallelized use) and plugged into the templated visitor objects.
 * most likely, these will reuse FacetDerivedData but not FacetPrimaryData.
 */
class FacetBoundedBody {
public:
    typedef ::vec_t vec_t;

    // vertices are mostly internal objects which exist
    // to conserve a little storage.
    // they should not be used widely by the visitors.
    // however, equal vertex id's are guaranteed at the intersection
    // of two neighbouring facets.
    // note that two topologically inconnected facets may share a
    // vertex id at singular points of the body boundary.
    // "edges" as such don't exist (except in the 2D case, where
    // the facets are in fact edges.)
    int num_vertices () const;
    int num_facets () const;
    const vec_t &facet_normal (int) const;
          double facet_volume (int) const;
    const vec_t &facet_vertex (int, int) const;
          int    facet_vertex_id (int, int) const;
    const vec_t  facet_center (int) const;
          int    facet_domain (int) const;
          void   facet_domain (int, int);
          int    facet_label  (int) const;
          void   facet_label  (int, int);
          double facet_inflection (int, int);
    #if (DIMENSION==3)
          double facet_gausscurv (int, int) const;
    #endif
          int    facet_neighbour_facet (int, int) const;

          int    insert_vertex (const vec_t &);
          int    insert_facet (const int []);

    // FIXME Domains interface

public:
    // see note in implementation
    void postcompute_neighbourships ();
    // signals that all facets have been inserted.
    // this also triggers the "inspection run" to find possible
    // problems and fills in the FacetDerivedData structures.
    // before this is called, several of the accessors above
    // return nonsense.
    void body_complete ();

    template <typename VISITOR>
    void visit_each_facet (VISITOR &) const;
    template <typename VISITOR>
    void visit_each_facet (VISITOR &);
    template <typename VISITOR>
    void visit_each_facet (const VISITOR &) const;
    template <typename VISITOR>
    void visit_each_facet (const VISITOR &);



private:
    std::vector <vec_t> coords_;
    std::vector <FacetPrimaryData> primar_;
    std::vector <FacetDerivedData> deriv_;
    friend struct AssignNeighboursVisitor;
    friend struct PrecomputeVisitor;
    friend struct Precompute2Visitor;
};


} // namespace Papaya2


//
// INLINES AND TEMPLATES
//

namespace Papaya2 {

    inline int FacetBoundedBody::num_facets () const {
        return primar_.size ();
    }

    inline int FacetBoundedBody::facet_label (int f) const {
        assert (f < num_facets ());
        assert (f >= 0);
        return primar_[f].label;
    }

    inline const vec_t &FacetBoundedBody::facet_normal (int f) const {
        assert (f < num_facets ());
        assert (f >= 0);
        return deriv_[f].normal;
    }

    inline int FacetBoundedBody::facet_neighbour_facet (int f, int n) const {
        assert (f < num_facets ());
        assert (f >= 0);
        assert (n >= 0);
        assert (n < DIMENSION);
        int ret = primar_[f].neigh[n];
        assert (ret == INVALID_FACET || (ret >= 0 && ret < num_facets ()));
        return ret;
    }

    inline int FacetBoundedBody::facet_vertex_id (int f, int i) const {
        assert (i >= 0);
        assert (i < DIMENSION);
        return primar_[f].vertex[i];
    }

    inline const vec_t &FacetBoundedBody::facet_vertex (int f, int i) const {
        return coords_[facet_vertex_id (f, i)];
    }

    inline int FacetBoundedBody::insert_vertex (const vec_t &coords) {
        int ret = coords_.size ();
        coords_.push_back (coords);
        return ret;
    }

    inline int FacetBoundedBody::insert_facet (const int verts_[]) {
        int ret = primar_.size ();
        FacetPrimaryData nf;
        nf.label = 0;
        nf.domain = INVALID_DOMAIN;
        for (int i = 0; i != DIMENSION; ++i) {
            assert (verts_[i] >= 0);
            assert (verts_[i] < (int)coords_.size ());
            nf.vertex[i] = verts_[i];
            nf.neigh[i] = INVALID_FACET;
        }
        primar_.push_back (nf);
        deriv_.push_back (FacetDerivedData ());
        return ret;
    }

    template <typename VISITOR>
    void FacetBoundedBody::visit_each_facet (VISITOR &vis) {
        for (int f = 0; f != num_facets (); ++f) {
            vis (*this, f);
        }
    }

    template <typename VISITOR>
    void FacetBoundedBody::visit_each_facet (const VISITOR &vis) {
        for (int f = 0; f != num_facets (); ++f) {
            vis (*this, f);
        }
    }

    template <typename VISITOR>
    void FacetBoundedBody::visit_each_facet (VISITOR &vis) const {
        for (int f = 0; f != num_facets (); ++f) {
            vis (*this, f);
        }
    }

    template <typename VISITOR>
    void FacetBoundedBody::visit_each_facet (const VISITOR &vis) const {
        for (int f = 0; f != num_facets (); ++f) {
            vis (*this, f);
        }
    }

    // return which of the edges of facet f is v0--v1 (in that order)
    // if none is, returns a negative integer.
    template <typename BODY>
    int body_facet_which_edge (const BODY &body, int f, int v0, int v1) {
        if (body.facet_vertex_id (f, 2) == v0) {
            if (body.facet_vertex_id (f, 0) == v1)
                return 2;
        } else if (body.facet_vertex_id (f, 1) == v0) {
            if (body.facet_vertex_id (f, 2) == v1)
                return 1;
        } else if (body.facet_vertex_id (f, 0) == v0) {
            if (body.facet_vertex_id (f, 1) == v1)
                return 0;
        }
        return -42;
    }
}

#endif // FACETBOUNDEDBODY
