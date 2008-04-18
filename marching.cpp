// implement marching squares algorithm
// for FIXME
// see http://en.wikipedia.org/wiki/Marching_squares
//


#include "util.h"
#include <stdarg.h>
#define MARCSQ_NO_LOG

class MarchingSquares {
public:
    void run (Boundary *, const Pixmap &dataset, Pixmap::val_t threshold, bool connectblack);

private:
    enum {
        UPPERLEFT  = 1,
        UPPERRIGHT = 2,
        LOWERLEFT  = 4,
        LOWERRIGHT = 8,
    };

    typedef Pixmap::val_t val_t;

    Boundary *boundary;
    Pixmap dataset;
    val_t threshold;
    bool connectblack;
    Pixmap visited;
    int dualxmax;
    int dualymax;
    int square_type (int x, int y);
    void trace_contour (int x, int y);
    void log (const char *fmt, ...);
    val_t lowerright (int x, int y);
    val_t upperright (int x, int y);
    val_t lowerleft  (int x, int y);
    val_t upperleft  (int x, int y);
};

void MarchingSquares::run (Boundary *b, const Pixmap &dataset_, Pixmap::val_t threshold_,
                           bool connectblack_) {
    assert (b);
    boundary = b;
    threshold = threshold_;
    connectblack = connectblack_;
    // create a copy of dataset, padded by a one pixel border
    // so that contours definitely end there
    log ("init\n");
    dataset.resize (dataset_.size1 () + 2, dataset_.size2 () + 2);
    for (int j = 0; j != dataset.size2 (); ++j)
    for (int i = 0; i != dataset.size1 (); ++i)
        dataset(i,j) = Pixmap::min_val ();
    for (int j = 0; j != dataset_.size2 (); ++j)
    for (int i = 0; i != dataset_.size1 (); ++i)
        dataset(i+1,j+1) = dataset_(i,j);
    // dual lattice sites are centered between pixel centers
    dualxmax = dataset.size1 () - 1;
    dualymax = dataset.size2 () - 1;
    // 1 = site has been visited, 0 = not
    visited.resize (dualxmax, dualymax);
    for (int j = 0; j != visited.size2 (); ++j)
    for (int i = 0; i != visited.size1 (); ++i)
        visited(i,j) = 0;
    // go through the data and look for a boundary that has not yet been
    // treated
    log ("traces\n");
    for (int ys = 0; ys != dualymax; ++ys)
    for (int xs = 0; xs != dualxmax; ++xs) {
        if (visited(xs,ys))
            continue;
        int type = square_type (xs, ys);
        switch (type) {
        // interior sites
        case 0:
        case UPPERLEFT|UPPERRIGHT|LOWERLEFT|LOWERRIGHT:
            continue;
        // never start on ambiguous sites.
        case UPPERLEFT|LOWERRIGHT:
        case UPPERRIGHT|LOWERLEFT:
            continue;
        }
        trace_contour (xs, ys);
    }
}

void MarchingSquares::trace_contour (int thisx, int thisy) {
    log ("trace_contour (%i, %i)\n", thisx, thisy);
    int prevx = -1;
    int prevy = -1;
    int prevvertex = Boundary::INVALID_VERTEX;
    int prevedge = Boundary::INVALID_EDGE;
    int initialedge = Boundary::INVALID_EDGE;
    while (!visited(thisx,thisy)) {
#define movedown  void ((++nexty, vertx = thisx+.5, verty = thisy+1.))
#define moveup    void ((--nexty, vertx = thisx+.5, verty = thisy))
#define moveright void ((++nextx, verty = thisy+.5, vertx = thisx+1.))
#define moveleft  void ((--nextx, verty = thisy+.5, vertx = thisx))
        int type = square_type (thisx, thisy);
        visited(thisx,thisy) = 1;
        int nextx = thisx;
        int nexty = thisy;
        double vertx, verty;
        // find out where to go, and where to put the next vertex.
        switch (type) {
        // straightforward cases
        case UPPERLEFT:
            moveup;
            break;
        case UPPERRIGHT:
            moveright;
            break;
        case LOWERLEFT:
            moveleft;
            break;
        case LOWERRIGHT:
            movedown;
            break;
        case UPPERLEFT|UPPERRIGHT:
            moveright;
            break;
        case LOWERLEFT|LOWERRIGHT:
            moveleft;
            break;
        case LOWERLEFT|UPPERLEFT:
            moveup;
            break;
        case LOWERRIGHT|UPPERRIGHT:
            movedown;
            break;
        case UPPERLEFT|LOWERLEFT|LOWERRIGHT:
            moveup;
            break;
        case UPPERRIGHT|LOWERLEFT|LOWERRIGHT:
            moveleft;
            break;
        case UPPERLEFT|UPPERRIGHT|LOWERLEFT:
            moveright;
            break;
        case UPPERLEFT|UPPERRIGHT|LOWERRIGHT:
            movedown;
            break;
        // ambiguous ones
        // these are never start sites, so prevx, prevy have sensible values.
        case UPPERLEFT|LOWERRIGHT:
            // ambiguous squares are traversed twice
            visited(thisx,thisy) = 0;
            if (connectblack) {
                if (prevx < thisx) {
                    moveup;
                } else if (prevx > thisx) {
                    movedown;
                } else {
                    die ("Entered UPPERLEFT|LOWERRIGHT from above or below.\n");
                }
            } else {
                if (prevx < thisx) {
                    movedown;
                } else if (prevx > thisx) {
                    moveup;
                } else {
                    die ("Entered UPPERLEFT|LOWERRIGHT from above or below.\n");
                }
            }
            break;
        case UPPERRIGHT|LOWERLEFT:
            // ambiguous squares are traversed twice
            visited(thisx,thisy) = 0;
            if (connectblack) {
                if (prevy > thisy) {
                    moveleft;
                } else if (prevy < thisy) {
                    moveright;
                } else {
                    die ("Entered UPPERRIGHT|LOWERLEFT from left or right.\n");
                }
            } else {
                if (prevy > thisy) {
                    moveright;
                } else if (prevy < thisy) {
                    moveleft;
                } else {
                    die ("Entered UPPERRIGHT|LOWERLEFT from left or right.\n");
                }
            }
            break;
        default:
            die ("Lost contact to contour. This should not have happened (type = %i).", type);
        }
        // add new vertex
        int thisvertex = boundary->insert_vertex (vertx, verty);
        // add new edge
        if (prevvertex != Boundary::INVALID_VERTEX) {
            int thisedge = boundary->insert_edge (
                prevedge, prevvertex, thisvertex, Boundary::INVALID_EDGE);
            if (initialedge == Boundary::INVALID_EDGE)
                initialedge = thisedge;
            prevedge = thisedge;
        }
        prevx = thisx;
        prevy = thisy;
        thisx = nextx;
        thisy = nexty;
        prevvertex = thisvertex;
#undef moveleft
#undef moveright
#undef moveup
#undef movedown
    }
    // encountered start site of this contour.
    // close contour.
    log ("close contour, initialedge = %i\n", initialedge);
    boundary->insert_edge (prevedge, prevvertex, Boundary::INVALID_VERTEX, initialedge);
}

inline int MarchingSquares::square_type (int x, int y) {
    // determine which kind of boundary we found
    // in type, all "white" pixels are set
    int type;
    type  = (upperleft  (x, y) > threshold) * UPPERLEFT;
    type |= (upperright (x, y) > threshold) * UPPERRIGHT;
    type |= (lowerleft  (x, y) > threshold) * LOWERLEFT;
    type |= (lowerright (x, y) > threshold) * LOWERRIGHT;
    return type;
}

inline MarchingSquares::val_t MarchingSquares::upperleft (int x, int y) {
    return dataset (x, y);
}

inline MarchingSquares::val_t MarchingSquares::upperright (int x, int y) {
    return dataset (x+1, y);
}

inline MarchingSquares::val_t MarchingSquares::lowerleft (int x, int y) {
    return dataset (x, y+1);
}

inline MarchingSquares::val_t MarchingSquares::lowerright (int x, int y) {
    return dataset (x+1, y+1);
}

inline void MarchingSquares::log (const char *fmt, ...) {
#ifndef MARCSQ_NO_LOG
    va_list va;
    va_start (va, fmt);
    vfprintf (stderr, fmt, va);
#endif
}


#include <iostream>
#include <fstream>
#include "minkval.h"
#include "tinyconf.h"

int main () {
    Configuration conf ("test.conf");
    Pixmap p;
    //load_test_pixmap (&p);
    load_pgm (&p, conf.string ("input", "filename"));
    if (conf.boolean ("input", "invert"))
        invert (&p);
    Boundary b;
    MarchingSquares m;
    bool connectblack = conf.boolean ("segment", "connectblack");
    //m.run (&b, p, 50, connectblack);
    load_poly (&b, "testdata/PolyExFuerSeb.poly");
    b.fix_contours ();
    std::ofstream contfile ("contours.out");
    dump_contours (contfile, b);
    calculate_all_surface_integrals (b);
    /*
    std::cerr << "W100 " << W100 (b) << std::endl;
    std::cerr << "W200 " << W200 (b) << std::endl;
    std::cerr << "W101 " << W101 (b) << std::endl;
    std::cerr << "W110 " << W110 (b) << std::endl;
    EigenSystem es;
    eigensystem (&es, 9., 4., 5, 1e-9);
    es.dump (std::cout);
    */
    return 0;
}

