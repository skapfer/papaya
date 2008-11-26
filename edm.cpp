
#include "util.h"

void euclidean_distance_map (Pixmap *edm, const Pixmap &p, double threshold) {
    assert (edm);

    // fill *edm with EDM of p.

    // p is not binarised as such, use p(x,y) > threshold to find
    // out whether the pixel is in the white phase.

#define binary_value(x,y) (p(x,y) > threshold)

    //int value_at_0_15 = binary_value (0, 15);

    abort ();
}

