// vim: et:sw=4:ts=4
// driver code

#include <iostream>
#include "util.h"
#include "minkval.h"
#include "tinyconf.h"

static bool ends_with (const std::string &s1, const std::string &s2) {
    if (s1.size () < s2.size ())
        return false;
    return std::string (s1, s1.size () - s2.size (), s2.size ()) == s2;
}

int main () {
    Configuration conf ("test.conf");
    Pixmap p;

    std::string filename = conf.string ("input", "filename");
    Boundary b;

    if (ends_with (filename, ".poly")) {
        load_poly (&b, filename);
        bool runfix   = conf.boolean ("polyinput", "fix_contours");
        bool forceccw = conf.boolean ("polyinput", "force_counterclockwise");
        if (runfix)
            fix_contours (&b, conf.boolean ("polyinput", "silent_fix_contours"));
        if (forceccw)
            force_counterclockwise_contours (&b);
    } else if (ends_with (filename, ".pgm")) {
        load_pgm (&p, filename);
        if (conf.boolean ("segment", "invert"))
            invert (&p);
        int  threshold    = conf.integer ("segment", "threshold");
        bool connectblack = conf.boolean ("segment", "connectblack");
        marching_squares (&b, p, threshold, connectblack);
    }

    std::string labcrit = conf.string ("output", "labels");
    if (labcrit == "none")
        label_none (&b);
    else if (labcrit == "by_contour")
        label_by_contour_index (&b);
    else if (labcrit == "by_component")
        label_by_component (&b);
    else if (labcrit == "by_domain") {
        rect_t r;
        r.top    = conf.floating ("domains", "clip_top");
        r.right  = conf.floating ("domains", "clip_right");
        r.bottom = conf.floating ("domains", "clip_bottom");
        r.left   = conf.floating ("domains", "clip_left");
        int xdomains = conf.integer ("domains", "xdomains");
        int ydomains = conf.integer ("domains", "ydomains");
        label_by_domain (&b, r, xdomains, ydomains);
    }


    std::string contfile ("contours.out");
    dump_contours (contfile, b, 1);
    dump_labels ("labels", b);
    calculate_all_surface_integrals (b);

    return 0;
}
