// vim: et:sw=4:ts=4

#include "util.h"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <iomanip>

static const double VERTEX_MERGE_TOLERANCE = 1e-6;


void die (const char *fmt, ...) {
    va_list al;
    va_start (al, fmt);
    vfprintf (stderr, fmt, al);
    fputs ("\n", stderr);
    abort ();
}

Pixmap::Pixmap ()
    : my_xdim (0), my_ydim (0) { }

void Pixmap::resize (int dimx, int dimy) {
    my_data.resize ((my_xdim = dimx) * (my_ydim = dimy));
}

void Pixmap::init_zero () {
    memset ((void *)(&my_data[0]), 0, sizeof (val_t) * my_data.size ());
}

Boundary::Boundary () {
}

int Boundary::insert_vertex (const vec_t &v) {
    int ret = my_vert.size ();
    my_vert.push_back (v);
    return ret;
}

int Boundary::insert_edge (int a, int b, int c, int d) {
    int ret = my_edge.size ();
    // look up missing arguments
    if (c == INVALID_VERTEX) {
        assert (d != INVALID_EDGE); // need second vertex or next edge
        c = edge(d).vert0;
    }
    if (b == INVALID_VERTEX) {
        assert (a != INVALID_EDGE);  // need first vertex or prev edge
        b = edge(a).vert1;
    }
    // connect neighbours, join contours
    int contour = INVALID_CONTOUR;
    if (a != INVALID_EDGE) {
        assert (edge(a).next == INVALID_EDGE);
        edge(a).next = ret;
        contour = edge(a).contour;
    }
    if (d != INVALID_EDGE) {
        assert (edge(d).prev == INVALID_EDGE);
        edge(d).prev = ret;
        if (contour != INVALID_CONTOUR) {
            // merge contours edge(ret).contour and edge(d).contour
            if (contour != edge(d).contour) {
                erase_contour_by_id (edge(d).contour);
                merge_contour_asc (d, contour);
            }
        } else {
            contour = edge(d).contour;
        }
    }
    if (d == INVALID_EDGE && a == INVALID_EDGE) {
        my_contours.push_back (ret);
        contour = ret;
    }
    edge_t E;
    E.prev = a;
    E.next = d;
    E.vert0 = b;
    E.vert1 = c;
    E.contour = contour;
    E.label = 0;
    my_edge.push_back (E);
    return ret;
}

void Boundary::merge_contour_asc (int ed, int newc) {
    assert (edge(ed).prev == INVALID_EDGE);
    // go through everything from ed up, and change contour
    // do not call on closed contours
    for (;;) {
        edge(ed).contour = newc;
        ed = edge(ed).next;
        if (ed == INVALID_EDGE)
            return;
    }
}

#ifndef NDEBUG
// for fuzzy vector comparison
static double normed_diff (const vec_t &a, const vec_t &b) {
    vec_t d = a;
    d -= b;
    return d.norm ();
}
#endif // NDEBUG

// perform sanity checking on a Boundary.
// this is not public code for now.
static void assert_boundary (const Boundary &b) {
#ifndef NDEBUG
    Boundary::contour_iterator cit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        Boundary::edge_iterator eit     = b.edges_begin (cit),
                                eit_end = b.edges_end (cit);
        while (eit != eit_end) {
            Boundary::edge_iterator next_eit = eit;
            ++next_eit;
            if (! (normed_diff (b.edge_vertex1 (eit),
                                b.edge_vertex0 (next_eit))
                   < VERTEX_MERGE_TOLERANCE)) {
                die ("assert_boundary: contour is not continuous.");
            }

            eit = next_eit;
        }
    }
#endif // NDEBUG
}

Boundary::edge_iterator Boundary::remove_vertex1 (edge_iterator eit) {
    assert (this == eit.my_boundary);
    assert (eit->next != INVALID_EDGE);
    // normally, delete the next edge.
    // if that edge is the magic edge which gives the name to this contour,
    // delete this edge instead.
    // if the contour consists of just one vertex (degenerate
    // contour, fail loudly.
    if (is_self_referential (eit))
        die ("Deleting vertex of degenerate contour");
    if (eit->next == eit->contour) {
        // delete this edge object and its metadata.
        // this code path is probably totally untested.
        int saved_vert0 = eit->vert0;
        ++eit;
        take_edge_out (eit->prev);
        edge_t *peit = &edge(eit);
        peit->vert0 = saved_vert0;
        assert_valid_link_structure ();
        return eit;
    } else {
        // delete next edge object, but this object's metadata.
        edge_t saved = edge(eit->next);
        take_edge_out (eit->next);
        // copy data back.
        edge_t *peit = &edge(eit);
        peit->label = saved.label;
        peit->vert1 = saved.vert1;
        assert_valid_link_structure ();
        return eit;
    }
}

// take an edge object out of the loop
// the object is assumed to have successors on both sides.
// this function just does the re-wiring.
inline void Boundary::take_edge_out (int ed) {
    edge_t &E = edge(ed);
    edge_t &nextE = edge(E.next);
    edge_t &prevE = edge(E.prev);
    nextE.prev = E.prev;
    prevE.next = E.next;

#ifndef NDEBUG
    // write nonsense to a edge that should no longer be in use
    E.prev = E.next = INVALID_EDGE;
    E.vert0 = E.vert1 = INVALID_VERTEX;
    E.contour = INVALID_EDGE;
    E.label = -1;
#endif
}

void Boundary::split_edge (edge_iterator eit, const vec_t &newv) {
    int nedge_id = my_edge.size ();
    int nvert_id = insert_vertex (newv);
    edge_t E = *eit;
    E.vert0 = nvert_id;
    E.prev  = eit.my_position;
    if (E.next != INVALID_EDGE) {
        edge (E.next).prev = nedge_id;
    }
    my_edge.push_back (E);
    edge (eit).vert1 = nvert_id;
    edge (eit).next = nedge_id;
    // check that everything is alright
    assert_valid_link_structure ();
}

void Boundary::split_contour_inserting_edge (edge_iterator it0,
                                             edge_iterator it1) {
    it0.my_period = 0;
    it1.my_period = 0;
    // remove from lookup table
    int del_contour = it0->contour;
    erase_contour_by_id (del_contour);
    assert (it1->contour == del_contour);
    // open the contour
    edge_iterator it0_succ = it0;
    ++it0_succ;
    open_link_ (it0);
    edge_iterator it1_succ = it1;
    ++it1_succ;
    open_link_ (it1);
    // splice two new edges into the contours
    int upper_cid = insert_edge (it0.my_position, INVALID_VERTEX, INVALID_VERTEX, it1_succ.my_position);
    int lower_cid = insert_edge (it1.my_position, INVALID_VERTEX, INVALID_VERTEX, it0_succ.my_position);
    ++it0;
    ++it0;
    ++it1;
    ++it1;
    assert (it1_succ == it0);
    assert (it0_succ == it1);
    ++it0.my_period;
    ++it1.my_period;
    // set the proper contour id.
    set_contour_id_ (it1_succ, it0, upper_cid);
    set_contour_id_ (it0_succ, it1, lower_cid);
    assert (it1_succ.my_position != it0_succ.my_position);
    my_contours.push_back (lower_cid);
    my_contours.push_back (upper_cid);
    // check that everything is alright
    assert_valid_link_structure (upper_cid);
    assert_valid_link_structure (lower_cid);
}

void Boundary::merge_contours_inserting_edge (edge_iterator it0,
                                              edge_iterator it1) {
    it0.my_period = 0;
    it1.my_period = 0;
    // remove from lookup table
    int del_contour = it0->contour;
    erase_contour_by_id (del_contour);
    assert (it1->contour != del_contour);
    erase_contour_by_id (it1->contour);
    // open the contours
    edge_iterator it0_succ = it0;
    ++it0_succ;
    open_link_ (it0);
    edge_iterator it1_succ = it1;
    ++it1_succ;
    open_link_ (it1);
    // splice two new edges into the contours
    int upper_cid = insert_edge (it0.my_position, INVALID_VERTEX, INVALID_VERTEX, it1_succ.my_position);
                    insert_edge (it1.my_position, INVALID_VERTEX, INVALID_VERTEX, it0_succ.my_position);
    it1 = it0;
    ++it1.my_period;
    // set the proper contour id.
    set_contour_id_ (it0, it1, upper_cid);
    my_contours.push_back (upper_cid);
    // check that everything is alright
    assert_valid_link_structure (upper_cid);
}


void Boundary::set_contour_id_ (edge_iterator a, edge_iterator b, int cc) {
    for (; a != b; ++a)
        edge (a).contour = cc;
}

// remove the connection between *it and it successor
// (internal helper function)
void Boundary::open_link_ (edge_iterator it) {
    edge_t &E = edge (it);
    assert (E.next != INVALID_EDGE);
    edge_t &nE = edge (E.next);
    assert (nE.prev == it.my_position);
    E.next = INVALID_EDGE;
    nE.prev = INVALID_EDGE;
}

bool Boundary::is_self_referential (edge_iterator eit) const {
    edge_iterator eit2 = eit;
    ++eit;
    return (&*eit) == (&*eit2);
}

// check whether prev/next fields in the structure are set up
// correctly and whether each edge has two correct vertices.
void Boundary::assert_valid_link_structure () const {
#ifndef NDEBUG
    contour_iterator cit = contours_begin ();
    for (; cit != contours_end (); ++cit) {
        assert_valid_link_structure (*cit);
    }
#endif
}

void Boundary::assert_valid_link_structure (int edge_id) const {
#ifndef NDEBUG
    assert (edge_id != INVALID_EDGE);
    int terminal_edge = edge_id;
    bool found_naming_edge = false;
    int contour_id = edge (edge_id).contour;
    for (;;) {
        // check prev/next pointers
        assert (edge (edge_id).contour == contour_id);
        if (edge_id == contour_id)
            found_naming_edge = true;
        int prev_edge = edge_id;
        edge_id = edge (edge_id).next;
        if (edge_id == INVALID_EDGE)
            break;
        assert (edge (edge_id).prev == prev_edge);
        if (edge_id == terminal_edge)
            break;
    }
    assert (found_naming_edge);
    // the other direction (contour could be still open)
    edge_id = terminal_edge;
    for (;;) {
        int prev_edge = edge_id;
        edge_id = edge (edge_id).prev;
        if (edge_id == INVALID_EDGE)
            break;
        assert (edge (edge_id).next == prev_edge);
        if (edge_id == terminal_edge)
            break;
    }
#endif
}

void Boundary::erase_contour_by_id (int c) {
    std::vector <int>::iterator i = std::find (my_contours.begin (), my_contours.end (), c);
    assert (i != my_contours.end ());
    std::swap (*i, my_contours.back ());
    my_contours.erase (my_contours.end () - 1);
}

void Boundary::erase_contour_by_index (int cindex) {
    assert (cindex >= 0);
    assert (cindex < (int)my_contours.size ());
    std::vector <int>::iterator i = my_contours.begin () + cindex;
    std::swap (*i, my_contours.back ());
    my_contours.erase (my_contours.end () - 1);
}

Boundary::edge_iterator &Boundary::edge_iterator::operator++ () {
    assert (my_position != INVALID_EDGE);
    my_position = my_boundary->edge(my_position).next;
    if (my_position == my_boundary->edge(my_position).contour)
        ++my_period;
    return *this;
}

Boundary::edge_iterator &Boundary::edge_iterator::operator-- () {
    assert (my_position != INVALID_EDGE);
    if (my_position == my_boundary->edge(my_position).contour)
        --my_period;
    my_position = my_boundary->edge(my_position).prev;
    return *this;
}

bool operator!= (const Boundary::edge_iterator &lhs, const Boundary::edge_iterator &rhs) {
    return lhs.my_position != rhs.my_position || lhs.my_boundary != rhs.my_boundary
        || lhs.my_period != rhs.my_period;
}

bool operator== (const Boundary::edge_iterator &lhs, const Boundary::edge_iterator &rhs) {
    return !(lhs != rhs);
}

std::string Boundary::edge_iterator::to_string () const {
    char buf[200] = { 0 };
    int cont = INVALID_CONTOUR;
    if (my_position != INVALID_EDGE) {
        cont = my_boundary->edge(my_position).contour;
    }
    snprintf (buf, 199, "Boundary::edge_iterator (b = %p, contour = %i, edge = %i, period = %i)",
        (void *)my_boundary,
        cont, my_position, my_period);
    return buf;
}

double Boundary::edge_length (Boundary::edge_iterator it) const {
    vec_t dist = vertex(it->vert0);
    dist -= vertex(it->vert1);
    return dist.norm ();
}

double Boundary::inflection_after_edge (Boundary::edge_iterator it) const {
    // normalized tangent vectors
    vec_t tang0 = vertex(it->vert1);
    tang0 -= vertex(it->vert0);
    tang0 /= tang0.norm ();
    ++it;
    vec_t tang1 = vertex(it->vert1);
    tang1 -= vertex(it->vert0);
    tang1 /= tang1.norm ();
    // z component of cross prod.
    double sinphi = tang0[0] * tang1[1] - tang1[0] * tang0[1];
#ifndef NDEBUG
    if (isnan (sinphi)) {
        std::string msg;
        if (! (tang0.norm () > 1e-6))
            msg += "norm0 is close to zero or NaN\n";
        if (! (tang1.norm () < 1e-6))
            msg += "norm1 is close to zero or NaN\n";
        msg += "nan in Boundary::inflection_after_edge (%s), vector product\n";
        std::string itstr = it.to_string ();
        die (msg.c_str (), itstr.c_str ());
    }
#endif
    double cosphi = dot (tang0, tang1);
    double ret = atan2 (sinphi, cosphi);
#ifndef NDEBUG
    if (isnan (ret)) {
        die ("Boundary::inflection_after_edge: NaN in atan2");
    }
#endif // NDEBUG
    assert (ret <= M_PI);
    assert (ret >= -M_PI);
    return ret;
}

Boundary::vec_t Boundary::edge_tangent (Boundary::edge_iterator it) const {
    vec_t ret = vertex(it->vert1);
    ret -= vertex(it->vert0);
    ret /= ret.norm ();
    return ret;
}

Boundary::vec_t Boundary::edge_normal (Boundary::edge_iterator it) const {
    vec_t tmp = vertex(it->vert1);
    tmp -= vertex(it->vert0);
    tmp /= tmp.norm ();
    return vec_t (tmp[1], -tmp[0]);
}

Boundary::vec_t Boundary::edge_vertex0 (Boundary::edge_iterator it) const {
    return vertex (it->vert0);
}

Boundary::vec_t Boundary::edge_vertex1 (Boundary::edge_iterator it) const {
    return vertex (it->vert1);
}

void Boundary::fix_contours (bool silent) {
    assert_boundary (*this);
    if (!silent)
        std::cerr << "fix_contours";
    int deg_edges = 0;
    int deg_spikes = 0;
    std::vector <int> deg_con_indices;
    Boundary::contour_iterator cit = contours_begin ();
    Boundary::contour_iterator cit_end = contours_end ();
    Boundary::edge_iterator eit, eit_end;
    while (cit != cit_end) {
        bool fixed_something = false;
        eit = edges_begin (cit);
        eit_end = edges_end (cit);
        // remove spikes with exterior angle = pi
        while (eit != eit_end) {
            edge_iterator eit_next = eit;
            ++eit_next;

            vec_t d = edge_vertex1 (eit_next);
            d -= edge_vertex0 (eit);

            if (d.norm () < VERTEX_MERGE_TOLERANCE) {
                if (is_self_referential (eit)) {
                    // further processing is dangerous
                    goto deg_contour;
                }
                // remove spike.  the norm-zero edge created by this
                // operation is removed later
                eit = remove_vertex1 (eit);
                ++deg_spikes;
                --eit;
                // remove_vertex1 could have removed our end
                eit_end = edges_end (cit);
                fixed_something = true;
            } else {
                eit = eit_next;
            }
        }
        eit = edges_begin (cit);
        // remove norm-zero edges
        while (eit != eit_end) {
            double len = edge_length (eit);
            if (len < VERTEX_MERGE_TOLERANCE) {
                if (is_self_referential (eit)) {
                    // further processing is dangerous
                    goto deg_contour;
                }
                ++deg_edges;
                eit = remove_vertex1 (eit);
                // remove_vertex1 could have removed our end
                eit_end = edges_end (cit);
                fixed_something = true;
            } else {
                ++eit;
            }
        }

        if (!fixed_something)
            ++cit;
        continue;

    deg_contour:
        // current contour has been reduced to or was a single-vertex
        // contour.  we mark it for later removal.
        deg_con_indices.push_back (cit - contours_begin ());
        ++cit;
    }
    int deg_contours = (int)deg_con_indices.size ();
    // remove degenerate contours
    {
        std::vector <int>::reverse_iterator it;
        it = deg_con_indices.rbegin ();
        for (; it != deg_con_indices.rend (); ++it) {
            erase_contour_by_index (*it);
        }
    }

    // final "status report"
    if (!silent) {
        std::cerr << "\n"
                  << std::setw (4) << deg_contours << " deg. contours (removed)\n"
                  << std::setw (4) << deg_spikes << " deg. spikes (removed)\n"
                  << std::setw (4) << deg_edges << " deg. edges (removed)\n";
    }
    assert_boundary (*this);
}

void force_counterclockwise_contours (Boundary *b) {
    Boundary::contour_iterator cit;
    double total_inflection_for_contour (const Boundary &b, Boundary::contour_iterator);
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
        if (total_inflection_for_contour (*b, cit) < 0)
            b->reverse_contour (cit);
    }
}

void Boundary::reverse_contour (Boundary::contour_iterator cit) {
    edge_iterator eit     = edges_begin (cit);
    edge_iterator eit_end = edges_end (cit);
    for (; eit != eit_end; ++eit) {
        edge_t &E = edge (eit);
        std::swap (E.prev, E.next);
        std::swap (E.vert0, E.vert1);
    }
}

void dump_vertex (std::ostream &os, int vertex, const Boundary &b) {
    const Boundary::vec_t &v = b.vertex (vertex);
    os << std::setprecision (18) << v.x () << " "
       << std::setprecision (18) << v.y () << "\n";
}

static double total_inflection_for_contour_unchecked (
        const Boundary *b, Boundary::contour_iterator cit) {
    Boundary::edge_iterator eit = b->edges_begin (cit),
        eit_end = b->edges_end (cit);
    double acc = 0.;
    for (; eit != eit_end; ++eit) {
        acc += b->inflection_after_edge (eit);
    }
    //fprintf (stderr, "%i: %f\n", *cit, acc/2/M_PI);
    return acc;
}

double total_inflection_for_contour (const Boundary *b, Boundary::contour_iterator cit) {
    double acc = total_inflection_for_contour_unchecked (b, cit);
    // total inflection should be +/-2pi.
    if (! (fabs (fabs (acc) - 2*M_PI) < 1e-3)) {
        Boundary::edge_iterator eit = b->edges_begin (cit),
            eit_end = b->edges_end (cit);
        int ctr = 0;
        for (eit = b->edges_begin (cit); eit != eit_end; ++eit) {
            vec_t v0 = b->edge_vertex0 (eit);
            vec_t v1 = b->edge_vertex1 (eit);
            fprintf (stderr, "infl at edge %.4i = %f\n"
                "(%f,%f) (%f,%f)\n", ctr++,
                b->inflection_after_edge (eit),
                v0[0], v0[1], v1[0], v1[1]);
        }
        die ("total_inflection_for_contour (contour_id = %i):\n%.20e\n%.20e",
            *cit, acc, 2*M_PI);
    }
    return acc;
}

double total_inflection_for_contour (const Boundary &b, Boundary::contour_iterator cit) {
    return total_inflection_for_contour (&b, cit);
}

#ifndef NDEBUG

template <typename BOUNDARY, typename VISITOR>
static void visit_possibly_unclosed_contour (BOUNDARY &b, 
                                             Boundary::contour_iterator cit,
                                             VISITOR &vis) {
    Boundary::edge_iterator eit_begin = b.edges_begin (cit);
    Boundary::edge_iterator eit = b.edges_end (cit);
    // if the contour is not closed, we have to find its start edge
    for (;;) {
        if (eit == eit_begin) {
            // contour is closed, we looped once around
            eit = b.edges_begin (cit);
            for (; eit != b.edges_end (cit); ++eit) {
                vis (b, eit);
            }
            return;
        }
        if (b.edge_has_predecessor (eit)) {
            --eit;
        } else {
            // contour is not closed, we got to the beginning
            for (;;) {
                vis (b, eit);
                if (b.edge_has_successor (eit))
                    ++eit;
                else
                    return;
            }
        }
    }
}

template <typename BOUNDARY, typename VISITOR>
static void visit_possibly_incomplete_boundary (BOUNDARY &b,
                                                VISITOR &vis) {
    typename BOUNDARY::contour_iterator cit = b.contours_begin ();
    for (; cit != b.contours_end (); ++cit) {
        visit_possibly_unclosed_contour (b, cit, vis);
    }
}

static void assert_sensible_visitor (const Boundary &b, 
                                     Boundary::edge_iterator eit) {
    assert (b.edge_length (eit) > VERTEX_MERGE_TOLERANCE);
    // FIXME other checks.
}

void assert_sensible_boundary (const Boundary &b) {
    visit_possibly_incomplete_boundary (b, assert_sensible_visitor);
}

#endif // NDEBUG

void dump_contours (const std::string &filename, const Boundary &a, int flags) {
    std::ofstream os ((filename + ".out").c_str ());
    dump_contours (os, a);
    std::ofstream scr ((filename + ".gp").c_str ());
    scr << "plot \"" << filename << ".out\" w lp\n"; 
}

enum { BY_DIRECTION = 1 };

void dump_contours (std::ostream &os, const Boundary &a, int flags) {
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit, eit_end;
    for (cit = a.contours_begin (); cit != a.contours_end (); ++cit) {
        if (flags & BY_DIRECTION) {
            if (total_inflection_for_contour (&a, cit) < 0)
                continue;
        }
        eit = a.edges_begin (cit);
        eit_end = a.edges_end (cit);
        ++eit_end;
        for (; eit != eit_end; ++eit) {
            dump_vertex (os, eit->vert0, a);
        }
        os << "\n"; // contour sep. (for gnuplot)
    }
    os << "\n"; // index sep. (for gnuplot)
    if (! (flags & BY_DIRECTION))
        return;
    for (cit = a.contours_begin (); cit != a.contours_end (); ++cit) {
        if (flags & BY_DIRECTION) {
            if (total_inflection_for_contour (&a, cit) > 0)
                continue;
        }
        eit = a.edges_begin (cit);
        eit_end = a.edges_end (cit);
        ++eit_end;
        for (; eit != eit_end; ++eit) {
            dump_vertex (os, eit->vert0, a);
        }
        os << "\n"; // contour sep. (for gnuplot)
    }
    os << "\n"; // index sep. (for gnuplot)
}

void invert (Pixmap *p) {
    Pixmap::val_t off = p->max_val () - p->min_val ();
    for (int y = 0; y != p->size2 (); ++y)
    for (int x = 0; x != p->size1 (); ++x)
        (*p)(x,y) = off - (*p)(x,y);
}

void eigensystem_symm (EigenSystem *sys, double a1, double off, double b2) {
    assert_not_nan (a1);
    assert_not_nan (b2);
    assert_not_nan (off);
    double evp = a1*b2 - off*off;
    double b = a1 + b2;
    double q = b*b - 4.*evp;
    // should be non-negative for symmetric matrices.
    assert_not_nan (q);
    //assert (q >= -1e-12);
    if (! (q >= -1e-12))
        // this qty sometimes has a relatively large error, especially for 
        // near-diagonal matrices. we probaly need sth. better here.
        fprintf (stderr, "warning: a1 = %f; b2 = %f; off = %f; => q = %f\n", a1, b2, off, q);
    if (q < 0.)
        q = 1e-42;
    else 
        q = sqrt (q) + 1e-42;
    double evr;
    if (b < 0.)
        evr = b - q;
    else
        evr = b + q;
    evr *= .5;
    sys->eval[0] = evr;
    sys->eval[1] = evp/evr;
    if (b > 0.)
        std::swap (sys->eval[0], sys->eval[1]);
    vec_t v;
    v[0] = a1 - b2;
    if (v[0] < 0.)
        v[0] -= q;
    else {
        v[0] += q;
        std::swap (sys->eval[0], sys->eval[1]);
    }
    v[1] = 2.*off;
    v.normalize ();
    sys->evec[0] = v;
    sys->evec[1] = rot90_ccw (v);
}

#ifndef NDEBUG
void EigenSystem::dump (std::ostream &os) {
    for (int i = 0; i != 2; ++i) {
        os << std::setw (7) << eval[i] << "("
           << std::setw (7) << evec[i][0] << ", "
           << std::setw (7) << evec[i][1] << ")\n";
    }
}
#endif


