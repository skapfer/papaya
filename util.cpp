// vim: et:sw=4:ts=4

#include "util.h"
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iomanip>

void die (const char *fmt, ...) {
    va_list al;
    va_start (al, fmt);
    vfprintf (stderr, fmt, al);
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
                erase_contour (edge(d).contour);
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

inline Boundary::edge_t &Boundary::edge (edge_iterator eit) {
    assert (this == eit.my_boundary);
    return edge(eit.my_position);
}

Boundary::edge_iterator Boundary::remove_vertex1 (edge_iterator eit) {
    assert (this == eit.my_boundary);
    // normally, delete the next edge.
    // if that edge is the magic edge which gives the name to this contour,
    // delete this edge instead.
    // if the contour consists of just one vertex (degenerate
    // contour, fail loudly.
    if (eit->next == eit->contour) {
        // delete this edge object and its metadata.
        int saved_vert0 = eit->vert0;
        ++eit;
        take_edge_out (eit->prev);
        edge_t *peit = &edge(eit);
        peit->vert0 = saved_vert0;
    } else {
        // delete next edge object, but this object's metadata.
        edge_t saved = edge(eit->next);
        take_edge_out (eit->next);
        // copy data back.
        edge_t *peit = &edge(eit);
        peit->label = saved.label;
        peit->vert1 = saved.vert1;
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

bool Boundary::is_self_referential (edge_iterator eit) const {
    edge_iterator eit2 = eit;
    ++eit;
    return (&*eit) == (&*eit2);
}

void Boundary::erase_contour (int c) {
    std::vector <int>::iterator i = std::find (my_contours.begin (), my_contours.end (), c);
    assert (i != my_contours.end ());
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
    double cz = tang0[0] * tang1[1] - tang1[0] * tang0[1];
#ifndef NDEBUG
    if (isnan (cz)) {
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
    double ret = asin (cz);
#ifndef NDEBUG
    if (isnan (ret)) {
        die ("nan in Boundary::inflection_after_edge, arcsine\n");
    }
#endif
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
    if (!silent)
        std::cerr << "fix_contours";
    int deg_contours = 0;
    int deg_edges = 0;
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit, eit_end;
    for (cit = contours_begin (); cit != contours_end (); ++cit) {
        eit = edges_begin (cit);
        eit_end = edges_end (cit);
        // remove norm-zero edges
        while (eit != eit_end) {
            double len = edge_length (eit);
            if (len < 1e-6) {
                if (is_self_referential (eit)) {
                    ++deg_contours;
                    break;
                }
                ++deg_edges;
                eit = remove_vertex1 (eit);
                // remove_vertex1 could have removed our end
                eit_end = edges_end (cit);
            } else {
                ++eit;
            }
        }
    }
    if (!silent) {
        std::cerr << "\n "
                  << deg_contours << " deg. contours (NOT REMOVED)\n "
                  << deg_edges << " deg. edges (removed)\n";
    }
}

void dump_vertex (std::ostream &os, int vertex, const Boundary &b) {
    const Boundary::vec_t &v = b.vertex (vertex);
    os << std::setprecision (18) << v.x () << " "
       << std::setprecision (18) << v.y () << "\n";
}

double total_inflection_for_contour (const Boundary *b, Boundary::contour_iterator cit) {
    Boundary::edge_iterator eit = b->edges_begin (cit),
        eit_end = b->edges_end (cit);
    double acc = 0.;
    for (; eit != eit_end; ++eit) {
        acc += b->inflection_after_edge (eit);
    }
    // total inflection should be +/-2pi.
    if (! (fabs (fabs (acc) - 2*M_PI) < 1e-3)) {
        die ("total_inflection_for_contour:\n%.20e\n%.20e", acc, 2*M_PI);
    }
    return acc;
}

void dump_contours (std::ostream &os, const Boundary &a) {
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit, eit_end;
    for (cit = a.contours_begin (); cit != a.contours_end (); ++cit) {
        if (total_inflection_for_contour (&a, cit) < 0)
            continue;
        eit = a.edges_begin (cit);
        eit_end = a.edges_end (cit);
        ++eit_end;
        for (; eit != eit_end; ++eit) {
            dump_vertex (os, eit->vert0, a);
        }
        os << "\n"; // contour sep. (for gnuplot)
    }
    os << "\n"; // index sep. (for gnuplot)
    for (cit = a.contours_begin (); cit != a.contours_end (); ++cit) {
        if (total_inflection_for_contour (&a, cit) > 0)
            continue;
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

void load_test_pixmap (Pixmap *p) {
    const char *bla =
       /*01234567890123456789012345678901234567890123456789012345678901234567890123456789*/
        "                                                                                "
        "             XXXXXXXXXXXX                                                       "
        "               XXXXXXXXXXXXXXX                                                  "
        "               XXXX     XXXXXXXX                    XXXXXXXXXXXXXXXX            "
        "              XXXXXXXXXXXXXX                        XXX          XXX            "
        "                  XXXXXXX       XXXXX               XXX          XXX            "
        "                         XXXXXXXXXXXXXXXX           XXXXXXXXXXXXXXXX            "
        "                           XXXX XXXX XXX                                        "
        "                          XXXX XXXXXX X                                         "
        "                              XXXXXXXX                                          "
        "                                                                                "
        "                                   XXX                                          "
        "                               XXXXXXXXXX                                       "
        "                                XXXXXXX                                         "
        "                                                                                "
       /*01234567890123456789012345678901234567890123456789012345678901234567890123456789*/
        "                                                                                "
        "                                                                                ";
    int xdim = 80;
    int ydim = strlen (bla);
    assert (! (ydim % xdim));
    ydim /= xdim;
    p->resize (xdim, ydim);
    for (int y = 0; y != ydim; ++y)
    for (int x = 0; x != xdim; ++x)
        (*p)(x,y) = (bla[xdim*y + x] != ' ') ? Pixmap::max_val () : Pixmap::min_val ();
}

void invert (Pixmap *p) {
    long off = p->max_val () - p->min_val ();
    for (int y = 0; y != p->size2 (); ++y)
    for (int x = 0; x != p->size1 (); ++x)
        (*p)(x,y) = off - (*p)(x,y);
}

// get eigenvalues of a 2x2 matrix [if they are real].
static vec_t eigenvalues (double a1, double a2, double b1, double b2) {
    double c = a1*b2 - a2*b1;
    double b = a1 + b2;
    b *= -1.;
    double q = b*b - 4.*c;
    if (q < 0.) {
        // ignore small complex contributions.
        if (q > -1e-12)
            q = 0.;
        else
            die ("complex eigenvalues in eigenvalues (q = %g)", q);
    }
    q = sqrt (q);
    if (b < 0.) q *= -1.;
    q += b;
    q *= -.5;
    return vec_t (q, c/q);
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

void eigensystem (EigenSystem *sys, double a1, double a2, double b1, double b2) {
    vec_t ev = eigenvalues (a1, a2, b1, b2);
    sys->eval[0] = ev[0];
    sys->eval[1] = ev[1];
    for (int i = 0; i != 2; ++i) {
        double X = sys->eval[i];
        double d1 = a2 / (X-a1);
        double d2 = b1 / (X-b2);
        if (fabs (d1) < fabs (d2)) {
            sys->evec[i][0] = d1;
            sys->evec[i][1] = 1.;
        } else {
            sys->evec[i][0] = 1.;
            sys->evec[i][1] = d2;
        }
    }
}


