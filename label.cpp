// vim:et:sw=4:ts=4
// attach labels

#include "util.h"

// update the label of each edge on the contour to be "label"
static void relabel_contour (Boundary *b, Boundary::contour_iterator cit, int label) {
    Boundary::edge_iterator eit;
    for (eit = b->edges_begin (cit); eit != b->edges_end (cit); ++eit) {
        b->edge_label (eit, label);
    }
}

// labels by contour index, i.e. running number [0...num_contours-1]
// this is _not_ the contour id in the edge_t::contour field
// FIXME this is ugly, merge these two concepts
void label_by_contour_index (Boundary *b) {
    Boundary::contour_iterator cit;
    int l = 0;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit, ++l) {
        relabel_contour (b, cit, l);
    }
}

namespace {
    double total_inflection_for_contour (Boundary *b, Boundary::contour_iterator cit) {
        Boundary::edge_iterator eit = b->edges_begin (cit),
            eit_end = b->edges_end (cit);
        double acc = 0.;
        for (; eit != eit_end; ++eit) {
            acc += b->inflection_after_edge (eit);
        }
        // total inflection should be +/-2pi.
        assert (fabs (fabs (acc) - 2*M_PI) < 1e-10);
        return acc;
    }

    class LineSearchObject {
    public:
        // construct a LineSearchObject.
        // this is an expensive operation which should speed up
        // searching later on
        LineSearchObject (const Boundary *bnd, const vec_t &ray_direction) {
            b = bnd;
            assert (b);
            dir = ray_direction;
            dir.normalize ();
            assert ((dir.norm () - 1.) < 1e-10);
            cdir = rot90_ccw (dir);
            data.reserve (bnd->num_edges ());
            Boundary::contour_iterator cit;
            Boundary::edge_iterator eit;
            edge_data_t d;
            for (cit = bnd->contours_begin (); cit != bnd->contours_end (); ++cit)
            for (eit = bnd->edges_begin (cit); eit != bnd->edges_end (cit); ++eit) {
                double nc0 = dot (dir, bnd->edge_vertex0 (eit));
                double nc1 = dot (dir, bnd->edge_vertex0 (eit));
                d.it = eit;
                d.max_normalcoord = std::max (nc0, nc1);
                this->data.push_back (d);
            }
            std::sort (data.begin (), data.end (), by_max_normalcoord);
        }

        typedef Boundary::edge_t edge_t;

        struct edge_data_t {
            Boundary::edge_iterator it;
            double max_normalcoord;
        };

        class iterator;

        iterator edges_intersecting_ray_begin (const vec_t &ray_0);
        iterator end () { return iterator (0, data.end (), vec_t (0., 0.)); }

    private:
        const Boundary *b;
        vec_t dir, cdir;
        std::vector <edge_data_t> data;
        static inline bool by_max_normalcoord (const edge_data_t &lhs, const edge_data_t &rhs) {
            return lhs.max_normalcoord < rhs.max_normalcoord;
        }
    public:
        // this should implement the basic requirements for a forward iterator.
        // as internal state, it keeps an iterator into the "data" field of
        // LineSearchObject.
        // it steps through data, returning only edges that intersect
        // a ray starting in point ray_0 and having the direction
        // specified by parent->dir.
        // the normal coordinate is the coordinate in ray direction,
        // the conjugate coordinate is orthogonal to it.
        class iterator : private std::vector <edge_data_t>::const_iterator {
            typedef std::vector <edge_data_t>::const_iterator base_t;
        public:
            iterator ()
                : base_t (), parent (0) { }

            iterator (LineSearchObject *p, const base_t &itcore, const vec_t &ray_0)
                : base_t (itcore), r0 (ray_0), parent (p) { }

            iterator &operator++ () {
                check_valid ();
                do {
                    base_this ()->operator++ ();
                } while (*this != parent->data.end () && !does_edge_intersect ());
                return *this;
            }

            void operator++ (int) {
                this->operator++ ();
            }

            const edge_t *operator-> () const {
                check_valid ();
                return &operator* ();
            }

            const edge_t &operator* () const {
                check_valid ();
                return *(edge_data ()->it);
            }

            double max_normalcoord () const {
                check_valid ();
                return edge_data ()->max_normalcoord;
            }

            const vec_t &intersection_vertex () const {
                check_valid ();
                return ivtx;
            }

            friend bool operator!= (const iterator &lhs, const iterator &rhs) {
                // comparison ignores all iterator fields for reasons of
                // practicality. does not matter because only iterators within a
                // search should be compared
                return (*lhs.base_this ()) != rhs;
            }

        private:
            LineSearchObject *parent;
            const vec_t r0;
            vec_t ivtx;
            void check_valid () const {
                assert (parent);
                assert (*this != parent->data.end ());
            }
            const Boundary *b () const {
                return parent->b;
            }
            const vec_t &dir () const {
                return parent->dir;
            }
            const vec_t &cdir () const {
                return parent->cdir;
            }
            Boundary::edge_iterator edge_iter () const {
                return edge_data ()->it;
            }
            const edge_data_t *edge_data () const {
                return base_t::operator->();
            }
            const base_t *base_this () const {
                return this;
            }
            base_t *base_this () {
                return this;
            }
            bool does_edge_intersect ();
        };
    };

    LineSearchObject::iterator LineSearchObject::edges_intersecting_ray_begin (
        const vec_t &ray_0) {
        std::vector <edge_data_t>::const_iterator X;
        edge_data_t crit;
        crit.max_normalcoord = dot (ray_0, dir);
        X = std::lower_bound (data.begin (), data.end (),
                              crit, by_max_normalcoord);
        return iterator (this, X, ray_0);
    }

    bool LineSearchObject::iterator::does_edge_intersect () {
        // determine lambda so that v0 + lambda(v1-v0) in ray
        double cc0 = dot (cdir (), b ()->edge_vertex0 (edge_iter ()));
        double cc1 = dot (cdir (), b ()->edge_vertex1 (edge_iter ()));
        double ccR = dot (cdir (), r0);
        assert (!isnan (cc0));
        assert (!isnan (cc1));
        assert (!isnan (ccR));
        double lambda;
        if (cc0 <= ccR and ccR <= cc1) {
            double l1 = ccR - cc0;
            double l2 = cc1 - cc0;
            if (l2 <= l1)        // keep underflows in check
                lambda = 1.;     // (exact value is irrelevant)
            else
                lambda = l1/l2;
        } else if (cc1 <= ccR and ccR <= cc0) {
            double l1 = cc0 - ccR;
            double l2 = cc0 - cc1;
            if (l2 <= l1)
                lambda = 1.;
            else
                lambda = l1/l2;
        } else {
            return false;
        }
        assert (lambda <= 1.);
        assert (lambda >= 0.);
        // range of normal coordinates of current edge
        double nc0 = dot (dir (), b ()->edge_vertex0 (edge_iter ()));
        double ncT = dot (dir (), b ()->edge_vertex1 (edge_iter ()));
        double mu = nc0;
        ncT -= nc0;
        mu += lambda * ncT;
        double ncR = dot (dir (), r0);
        mu -= ncR;
        // now: mu: dir r0 + mu = dir (v0 + lambda(v1-v0))
        if (mu < 0.)
            return false;
        // update vertex-of-intersection since we have calculated
        // all the necessary things anyway
        ivtx = b ()->edge_vertex0 (edge_iter ());
        ivtx = lambda * (b ()->edge_vertex1 (edge_iter ()) 
                         - b ()->edge_vertex0 (edge_iter ()));
        return true;
    }
}

void label_by_component (Boundary *b) {
    // first, label everything by contour index.
    label_by_contour_index (b);
    // then, find clockwise contours which correspond to
    // interior boundary segments. these are assigned to the
    // first counterclockwise contour that is found via an
    // outward search.
    std::vector <int> ccw (b->max_contour_id (), -1);
    {
        // init the ccw array to indicate whether a contour is
        // counterclockwise
        Boundary::contour_iterator cit;
        int label = 0;
        for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
            bool is_ccw = (total_inflection_for_contour (b, cit) > 0.);
            if (is_ccw)
                // ccw contours have a label assigned now.
                // cw contours are merged into their containing ccw contour in
                // step 2.
                relabel_contour (b, cit, label++);
            ccw.at (*cit) = is_ccw;
        }
    }
    {
        // step 2
        vec_t ray_direction = vec_t (1., 0.);
        LineSearchObject ls (b, ray_direction);
        Boundary::contour_iterator cit;
        for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
            int ccwflag = ccw.at (*cit);
            assert (ccwflag == 0 || ccwflag == 1);
            if (!ccwflag) {
                // is clockwise
                vec_t ray_begin = b->edge_vertex0 (b->edges_begin (cit));
                LineSearchObject::iterator it =
                    ls.edges_intersecting_ray_begin (ray_begin);
                LineSearchObject::iterator end = ls.end ();
                for (; it != end; ++it) {
                    int found_ccwflag = ccw.at (it->contour);
                    assert (found_ccwflag == 0 || found_ccwflag == 1);
                    if (found_ccwflag) {
                        // found a counterclockwise countour which contains us.
                        // done.
                        relabel_contour (b, cit, it->label);
                        goto next_contour;
                    }
                }
                // we hit the outer boundary without finding a ccw contour.
                // most probably, this is a user error.
                die ("label_by_component: hit dataset boundary while searching for exterior boundary. "
                     "most probably, your polygons' vertices are in inverse order.");
            }
        next_contour:;
        }
    }
}

