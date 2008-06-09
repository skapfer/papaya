// vim: et:sw=4:ts=4
// study dependancy of tensors on aspect ratio of simple shapes


#include <iostream>
#include <fstream>
#include "tinyconf.h"
#include "minkval.h"

int main (int argc, const char **argv) {
    std::string infile = argv[1];
    std::string outfile = argv[2];
    Boundary b;
    load_poly (&b, infile);
    label_by_contour_index (&b);

    // calculate the functionals
    MatrixMinkowskiFunctional &w211 = *create_w211 ();
    MatrixMinkowskiFunctional &w220 = *create_w220 ();
    MatrixMinkowskiFunctional &w120 = *create_w120 ();
    MatrixMinkowskiFunctional &w020 = *create_w020 ();
    Boundary::contour_iterator cit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit) {
        w211.add_contour (b, b.edges_begin (cit), b.edges_end (cit));
        w220.add_contour (b, b.edges_begin (cit), b.edges_end (cit));
        w120.add_contour (b, b.edges_begin (cit), b.edges_end (cit));
        w020.add_contour (b, b.edges_begin (cit), b.edges_end (cit));
    }
    w211.dump (std::cout);
    w220.dump (std::cout);
    w120.dump (std::cout);
    w020.dump (std::cout);

    EigenSystem es;
    std::ofstream os (outfile.c_str ());
    //os << std::scientific;
    double rat;
    int j = 0;
    double rat_step = 1. / (b.contours_end () - b.contours_begin ());
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit, ++j) {
        rat = (j+1) * rat_step;
        os << std::setw (20) << std::setprecision (12) << rat;
        eigensystem_symm (&es, w211.value (j));
        rat = es.eval[0] / es.eval[1];
        if (rat > 1) rat = 1. / rat;
        os << std::setw (20) << std::setprecision (12) << rat;
        eigensystem_symm (&es, w220.value (j));
        rat = es.eval[0] / es.eval[1];
        if (rat > 1) rat = 1. / rat;
        os << std::setw (20) << std::setprecision (12) << rat;
        eigensystem_symm (&es, w120.value (j));
        rat = es.eval[0] / es.eval[1];
        if (rat > 1) rat = 1. / rat;
        os << std::setw (20) << std::setprecision (12) << rat;
        eigensystem_symm (&es, w020.value (j));
        rat = es.eval[0] / es.eval[1];
        if (rat > 1) rat = 1. / rat;
        os << std::setw (20) << std::setprecision (12) << rat;
        os << "\n";
    }
    delete &w211;
    delete &w220;
    delete &w120;
    delete &w020;
    return 0;
}

