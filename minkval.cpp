
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

// W200 = integral curvature
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


