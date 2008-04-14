
#include "minkval.h"

SurfaceMinkowskiValuator::SurfaceMinkowskiValuator () {
    my_ref = vec_t (0, 0);
    my_b = 0;
}

SurfaceMinkowskiValuator::~SurfaceMinkowskiValuator () {
}

void SurfaceMinkowskiValuator::set_reference (vec_t ref) {
    my_ref = ref;
}

void SurfaceMinkowskiValuator::add_contourseg (Boundary *b,
                                               edge_iterator begin,
                                               edge_iterator end) {
    my_b = b;
    // do sth.
    b = 0;
}

double SurfaceMinkowskiValuator::scalar (const string &name) {
    std::map <string, double>::const_iterator it;
    it = sacc.find (name);
    if (it != sacc.end ()) {
        return it->second;
    } else {
        die ("Invalid scalar key");
    }
}

vec_t SurfaceMinkowskiValuator::vector (const string &name) {
    std::map <string, vec_t>::const_iterator it;
    it = vacc.find (name);
    if (it != vacc.end ()) {
        return it->second;
    } else {
        die ("Invalid vector key");
    }
}

mat_t SurfaceMinkowskiValuator::matrix (const string &name) {
    std::map <string, mat_t>::const_iterator it;
    it = macc.find (name);
    if (it != macc.end ()) {
        return it->second;
    } else {
        die ("Invalid matrix key");
    }
}

vec_t SurfaceMinkowskiValuator::edge_normal (edge_iterator it) {
    assert (my_b);
    return my_b->edge_normal (it);
}

vec_t SurfaceMinkowskiValuator::edge_vertex0 (edge_iterator it) {
    assert (my_b);
    return my_b->edge_vertex0 (it);
}

vec_t SurfaceMinkowskiValuator::edge_vertex1 (edge_iterator it) {
    assert (my_b);
    return my_b->edge_vertex1 (it);
}




#if 0
// W100 = boundary length
double W100 (const Boundary &b) {
    double ret = 0.;
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        double acc = 0.;
        for (eit = b.edges_begin (cit); eit != b.edges_end (cit); ++eit) {
            acc += b.edge_length (eit);
        }
        ret += acc;
        std::cerr << " interm. result: " << acc << "\n";
    }
    return ret;
}

// W200 = surface integral curvature
double W200 (const Boundary &b) {
    double ret = 0.;
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        double acc = 0.;
        for (eit = b.edges_begin (cit); eit != b.edges_end (cit); ++eit) {
            acc += b.inflection_after_edge (eit);
        }
        ret += acc;
        std::cerr << " interm. result: " << acc/2/M_PI << " * 2pi\n";
    }
    return ret;
}

// W101 = surface integral normal -- zero by def.
vec_t W101 (const Boundary &b) {
    vec_t ret = vec_t (0., 0.);
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        vec_t acc = vec_t (0., 0.);
        for (eit = b.edges_begin (cit); eit != b.edges_end (cit); ++eit) {
            acc += b.edge_normal (eit) * b.edge_length (eit);
        }
        ret += acc;
        std::cerr << " interm. result: " << acc << "\n";
    }
    return ret;
}

// W110 = surface integral place
vec_t W110 (const Boundary &b) {
    vec_t ret = vec_t (0., 0.);
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        vec_t acc = vec_t (0., 0.);
        for (eit = b.edges_begin (cit); eit != b.edges_end (cit); ++eit) {
            vec_t avgvert = b.vertex (eit->vert1);
            avgvert += b.vertex (eit->vert0);
            avgvert /= 2;
            acc += avgvert * b.edge_length (eit);
        }
        ret += acc;
        std::cerr << " interm. result: " << acc << "\n";
    }
    return ret;
}

#endif


