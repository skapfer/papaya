// vim:et:sw=4:ts=4
// attach labels

#include "util.h"


void label_by_contour_index (Boundary *b) {
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    int l = 0;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit, ++l) {
        for (eit = b->edges_begin (cit); eit != b->edges_end (cit); ++eit) {
            b->edge_label (eit, l);
        }
    }
}

