// vim:et:sw=4:ts=4
// attach labels

#include "util.h"
#include "intersect.h"

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
    inline bool even (int i) {
        return ! (i%2);
    }
}

void label_by_component (Boundary *b) {
    double total_inflection_for_contour (const Boundary *b, Boundary::contour_iterator cit);

    // find clockwise contours which correspond to
    // interior boundary segments. these are assigned to the
    // first counterclockwise contour that is found via an
    // outward search.
    std::vector <int> ccw (b->max_contour_id (), -1);
    {
        // step 1
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
        Boundary::contour_iterator cit;
        for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
            int ccwflag = ccw.at (*cit);
            assert (ccwflag == 0 || ccwflag == 1);
            if (!ccwflag) {
                // is clockwise
                fprintf (stderr, "contour %i is clockwise.\n", *cit);
                vec_t ray_begin = b->edge_vertex0 (b->edges_begin (cit));
                vec_t ray_direction = b->edge_normal (b->edges_begin (cit));
                intersect_buffer_t info;
                intersect_ray_boundary (&info,
                    ray_begin, ray_direction,
                    b);
                // the buffer now contains a lot of intersections,
                // most of which are useless because they are caused by
                // * this contour
                // * other clockwise contours
                // * ccw contours which do not enclose this contour
                //   (e.g. enclosed interior components)
                // * ccw contours which do enclose this contour but are not the
                //   closest match

                // sort the hits so hits in the same contour are adjacent
                std::sort (info.begin (), info.end (), intersect_info_t::by_contour_id);
                intersect_buffer_t::iterator it = info.begin ();
                intersect_info_t *found_ = 0;
                while (it != info.end ()) {
                    intersect_buffer_t::iterator it_up =
                        std::upper_bound (it+1, info.end (), *it,
                            intersect_info_t::by_contour_id);
                    // ignore even contours
                    int num_hits = it_up - it;
                    if (even (num_hits))
                        goto next_intersect;
                    // ignore cw contours.
                    // [the contour we're processing is cw too, so self-intersects
                    //  will be detected here]
                    int this_contour = it->iedge->contour;
                    if (!ccw[this_contour])
                        goto next_intersect;
                    // find closest match in this contour
                    {
                        intersect_buffer_t::iterator it_candidate =
                            std::min_element (it, it_up,
                                    intersect_info_t::by_normal_coordinate);
                        if (found_) {
                            fprintf (stderr, "previous intersect exists\n");
                            if (found_->inc > it_candidate->inc)
                                found_ = &*it_candidate;
                        } else {
                            fprintf (stderr, "found intersect\n");
                            found_ = &*it_candidate;
                        }
                    }
                next_intersect:
                    it = it_up;
                }

                if (found_) {
                    relabel_contour (b, cit, found_->iedge->label);
                    goto next_contour;
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

#ifndef NDEBUG

#include <fstream>

static void dump_data (std::ostream &os, const Boundary &a,
                       Boundary::contour_iterator cit) {
    Boundary::edge_iterator eit, eit_end;
    eit = a.edges_begin (cit);
    eit_end = a.edges_end (cit);
    ++eit_end;
    for (; eit != eit_end; ++eit) {
        dump_vertex (os, eit->vert0, a);
    }
}

void dump_components (const std::string &filename, Boundary &b) {
    label_by_component (&b);
    std::ofstream ofdat ((filename + ".out").c_str (), std::ios::out);
    std::ofstream ofscr ((filename + ".gp").c_str (),  std::ios::out);
    ofscr << "unset key\n";
    ofscr << "plot \\\n";
    Boundary::contour_iterator cit = b.contours_begin ();
    int index = 0;
    dump_data (ofdat, b, cit);
    ofdat << "\n\n";  // index sep. (for gnuplot)
    int label = b.edges_begin (cit)->label;
    ofscr << "\t\"" << filename << ".out\" index "
          << index++ << " w lp lt " << label+1 << " pt " << label+1 << "\\\n";
    for (++cit; cit != b.contours_end (); ++cit) {
        dump_data (ofdat, b, cit);
        ofdat << "\n\n";  // index sep. (for gnuplot)
        int label = b.edges_begin (cit)->label;
        ofscr << "\t,\\\n";
        ofscr << "\t\"" << filename << ".out\" index "
              << index++ << " w lp lt " << label+1 << " pt " << label+1 << "\\\n";
    }
}
#endif // NDEBUG
