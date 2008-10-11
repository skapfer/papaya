#ifndef FACETBOUNDEDBODY 
#define FACETBOUNDEDBODY 

#ifndef PAPAYA2
#error "This header is PAPAYA2 code and cannot be mixed with PAPAYA1."
#endif // PAPAYA2
#include "util.h"


namespace Papaya2 {


enum { INVALID_FACET = -1 };

// data as read in
struct FacetPrimaryData {
    int vertex[DIMENSION];
    int neigh[DIMENSION];
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

    // FIXME Domains interface

public:
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
    friend struct PrecomputeVisitor;
    friend struct Precompute2Visitor;
};


} // namespace Papaya2


//
// INLINES AND TEMPLATES
//

namespace Papaya2 {

    int FacetBoundedBody::num_facets () const {
        return primar_.size ();
    }

    int FacetBoundedBody::facet_neighbour_facet (int f, int n) const {
        assert (f < num_facets ());
        assert (f >= 0);
        assert (n >= 0);
        assert (n < DIMENSION);
        int ret = primar_[f].neigh[n];
        assert (ret == INVALID_FACET || (ret >= 0 && ret < num_facets ()));
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
}

#endif // FACETBOUNDEDBODY
