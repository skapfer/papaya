// vim: et:sw=4:ts=4
// driver code

#include <iostream>
#include <fstream>
#include "util.h"
#include "minkval.h"
#include "tinyconf.h"

int main () {
    Configuration conf ("test.conf");
    Pixmap p;
    load_test_pixmap (&p);
    //load_pgm (&p, conf.string ("input", "filename"));
    if (conf.boolean ("input", "invert"))
        invert (&p);
    Boundary b;
    int  threshold    = conf.integer ("segment", "threshold");
    bool connectblack = conf.boolean ("segment", "connectblack");
    marching_squares (&b, p, threshold, connectblack);
    //load_poly (&b, "testdata/circle-d=1k.poly");
    //b.fix_contours ();
    std::ofstream contfile ("contours.out");
    dump_contours (contfile, b);
    dump_components ("components", b);
    calculate_all_surface_integrals (b);
    return 0;
}
