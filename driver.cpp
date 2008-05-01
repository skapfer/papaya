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
    //load_test_pixmap (&p);
    load_pgm (&p, conf.string ("input", "filename"));
    if (conf.boolean ("input", "invert"))
        invert (&p);
    Boundary b;
    bool connectblack = conf.boolean ("segment", "connectblack");
    marching_squares (&b, p, 50, connectblack);
    load_poly (&b, "testdata/circle-d=1k.poly");
    //b.fix_contours ();
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
