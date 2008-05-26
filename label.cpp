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

// reset all the labels to zero.
int label_none (Boundary *b) {
    Boundary::contour_iterator cit;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
        relabel_contour (b, cit, 0);
    }
    return 1;
}

// labels by contour index, i.e. running number [0...num_contours-1]
// this is _not_ the contour id in the edge_t::contour field
// FIXME this is ugly, merge these two concepts
int label_by_contour_index (Boundary *b) {
    Boundary::contour_iterator cit;
    int l = 0;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit, ++l) {
        relabel_contour (b, cit, l);
    }
    return b->num_contours ();
}

namespace {
    inline bool even (int i) {
        return ! (i%2);
    }
}

int label_by_component (Boundary *b) {
    double total_inflection_for_contour (const Boundary *b, Boundary::contour_iterator cit);

    // find clockwise contours which correspond to
    // interior boundary segments. these are assigned to the
    // first counterclockwise contour that is found via an
    // outward search.
    std::vector <int> ccw (b->max_contour_id (), -1);
    int ret;
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
        ret = label;
    }
    {
        // step 2
        Boundary::contour_iterator cit;
        for (cit = b->contours_begin (); cit != b->contours_end (); ++cit) {
            int ccwflag = ccw.at (*cit);
            assert (ccwflag == 0 || ccwflag == 1);
            if (!ccwflag) {
                // is clockwise
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
                            if (found_->inc > it_candidate->inc)
                                found_ = &*it_candidate;
                        } else {
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
    return ret;
}

static void introduce_divider (Boundary *b, const vec_t &line_0, const vec_t &line_dir) {
    intersect_buffer_t buff;
    intersect_line_boundary (&buff, line_0, line_dir, b);
    intersect_buffer_t::iterator it;
    for (it = buff.begin (); it != buff.end (); ++it) {
        b->split_edge (it->iedge, it->ivtx);
    }
}

int label_by_domain (Boundary *b, const rect_t &bbox, int divx, int divy) {
    // split lines reaching across domains
    double xstrip = bbox.right - bbox.left;
    assert (xstrip > 0.);
    xstrip /= divx;
    for (int i = 0; i <= divx; ++i)
        introduce_divider (b, vec_t (bbox.left + xstrip*i, 0.), vec_t (0., 1.));
    double ystrip = bbox.top - bbox.bottom;
    assert (ystrip > 0.);
    ystrip /= divy;
    for (int i = 0; i <= divy; ++i)
        introduce_divider (b, vec_t (0., bbox.bottom + ystrip*i), vec_t (1., 0.));
    // remove nonsense we just introduced
    fix_contours (b);
    // label everything
    Boundary::contour_iterator cit;
    Boundary::edge_iterator eit;
    for (cit = b->contours_begin (); cit != b->contours_end (); ++cit)
    for (eit = b->edges_begin (cit); eit != b->edges_end (cit); ++eit) {
        vec_t p = b->edge_vertex0 (eit);
        p += b->edge_vertex1 (eit);
        p /= 2.;
        if (intersect_vertex_rect (p, bbox)) {
            p -= vec_t (bbox.left, bbox.bottom);
            p.x () /= xstrip;
            p.y () /= ystrip;
            b->edge_label (eit, int (p.x ()) + int (p.y ()) * divy);
        } else {
            b->edge_label (eit, Boundary::NO_LABEL);
        }
    }
    return divx*divy;
}

#include <fstream>

void dump_labels (const std::string &filename, const Boundary &b) {
    std::ofstream ofdat ((filename + ".out").c_str (), std::ios::out);
    std::ofstream ofscr ((filename + ".gp").c_str (),  std::ios::out);
    ofscr << "unset key\n";
    ofscr << "plot \\\n";
    Boundary::contour_iterator cit = b.contours_begin ();
    int index = 0;
    for (; cit != b.contours_end (); ++cit) {
        Boundary::edge_iterator eit, eit_end;
        eit = b.edges_begin (cit);
        eit_end = b.edges_end (cit);
        ++eit_end;
        int label = b.edges_begin (cit)->label;

        while (eit != eit_end) {
            if (label != eit->label) {
                dump_vertex (ofdat, eit->vert0, b);
                ofdat << "\n\n";  // index sep. (for gnuplot)
                ofscr << "\t\"" << filename << ".out\" index "
                      << index++ << " w lp lt " << label+1 << " pt " << label+1 << "\\\n";
                ofscr << "\t,\\\n";
                label = eit->label;
            }
            dump_vertex (ofdat, eit->vert0, b);
            ++eit;
        }
        ofdat << "\n\n";  // index sep. (for gnuplot)
        ofscr << "\t\"" << filename << ".out\" index "
              << index++ << " w lp lt " << label+1 << " pt " << label+1 << "\\\n";
        ofscr << "\t,\\\n";
    }
    ofscr << "\t(1./0)\n";
}
