
#include "minkval.h"

//namespace {
//}

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

// W110 = surface integral place -- zero by def.
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




