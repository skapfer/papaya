// implement marching squares algorithm
// for converting a pixelized digital image into a set of contours
// see http://en.wikipedia.org/wiki/Marching_squares


#include "util.h"
#include <stdarg.h>
#define MARCSQ_NO_LOG

namespace {

class MarchingSquares {
    // just a dummy class to hold the algorithm's state
public:
    MarchingSquares (bool connect_body, bool periodic_data);

    void run (Boundary *, const Pixmap &dataset,
              Pixmap::val_t threshold);

private:
    enum {
        UPPERLEFT  = 1,
        UPPERRIGHT = 2,
        LOWERLEFT  = 4,
        LOWERRIGHT = 8
    };

    typedef Pixmap::val_t val_t;

    Boundary *boundary;
    Pixmap dataset;
    val_t threshold;
    const bool connect_void;
    const bool periodic_data;
    const int padding_shift;
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

MarchingSquares::MarchingSquares (bool connect_void_, bool periodic_data_)
    : connect_void (connect_void_), periodic_data (periodic_data_),
      padding_shift (periodic_data_ ? 2 : 1)  {
}

void MarchingSquares::run (Boundary *b, const Pixmap &dataset_,
                           Pixmap::val_t threshold_) {
    assert (b);
    boundary = b;
    threshold = threshold_;
    assert (threshold <= Pixmap::max_val ());
    assert (threshold >= Pixmap::min_val ());
    // create a copy of dataset, padded by a one pixel border
    // so that contours definitely end there
    log ("init\n");
    dataset.resize (dataset_.size1 () + 2*padding_shift,
                    dataset_.size2 () + 2*padding_shift);
    for (int j = 0; j != dataset.size2 (); ++j)
    for (int i = 0; i != dataset.size1 (); ++i)
        dataset(i,j) = Pixmap::min_val ();
    for (int j = 0; j != dataset_.size2 (); ++j)
    for (int i = 0; i != dataset_.size1 (); ++i)
        dataset(i+padding_shift,j+padding_shift) = dataset_(i,j);
    if (periodic_data) {
        // copy periodic wrappings in (incredibly messy)
        const int frow = dataset_.size2 () - 1;
        assert (frow+3 == dataset.size2 () - 2);
        for (int i = 0; i != dataset_.size1 (); ++i) {
            dataset(i+padding_shift,1)      = dataset_(i,frow);
            dataset(i+padding_shift,frow+3) = dataset_(i,0);
        }
        const int fcol = dataset_.size1 () - 1;
        assert (fcol+3 == dataset.size1 () - 2);
        for (int j = 0; j != dataset_.size2 (); ++j) {
            dataset(1,j+padding_shift)      = dataset_(fcol,j);
            dataset(fcol+3,j+padding_shift) = dataset_(0,j);
        }
        dataset(1,1)           = dataset_(fcol,frow);
        dataset(1,frow+3)      = dataset_(fcol,0);
        dataset(fcol+3,frow+3) = dataset_(0,0);
        dataset(fcol+3,1)      = dataset_(0,frow);
        for (int i = 0; i != dataset.size1 (); ++i) {
            assert (dataset(i,0) == Pixmap::min_val ());
            assert (dataset(i,dataset.size2 () - 1) == Pixmap::min_val ());
        }
        for (int j = 0; j != dataset_.size2 (); ++j) {
            assert (dataset(0,j) == Pixmap::min_val ());
            assert (dataset(dataset.size1 () - 1,j) == Pixmap::min_val ());
        }
    }
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
            if (connect_void) {
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
            if (connect_void) {
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
        // bitmap is top-down, but we need right-handed coordinates later
        int thisvertex = boundary->insert_vertex (
            vertx + .5 - padding_shift,
            dataset.size2() - verty - .5 - padding_shift);
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

}

void marching_squares (Boundary *b, const Pixmap &p,
                       Pixmap::val_t threshold,
                       bool connect_void, bool periodic_data) {
    MarchingSquares m (connect_void, periodic_data);
    m.run (b, p, threshold);
}
