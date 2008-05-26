// vim: et:sw=4:ts=4
// driver code

#include <iostream>
#include <fstream>
#include <getopt_pp_standalone.h>
#include "util.h"
#include "minkval.h"
#include "tinyconf.h"
using namespace GetOpt;

typedef AbstractMinkowskiFunctional **func_iterator;

static bool ends_with (const std::string &s1, const std::string &s2) {
    if (s1.size () < s2.size ())
        return false;
    return std::string (s1, s1.size () - s2.size (), s2.size ()) == s2;
}

static std::string basename (const std::string &str) {
    int i, j = -1;
    for (i = 0; i != (int)str.size (); ++i) {
        if (str[i] == '/')
            j = i;
    }
    if (j == (int)str.size ())
        return str;
    else
        return std::string (str, j+1, str.npos);
}

static void calculate_functional (AbstractMinkowskiFunctional *p,
                                  const Boundary &b) {
    Boundary::contour_iterator cit;
    for (cit = b.contours_begin (); cit != b.contours_end (); ++cit)
        p->add_contour (b, b.edges_begin (cit), b.edges_end (cit));
}

static void set_refvert_com (func_iterator begin, func_iterator end,
                             const Boundary &b, int num_labels) {
    VectorMinkowskiFunctional *w010 = create_w010 ();
    calculate_functional (w010, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w010->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
}

static void set_refvert_cos (func_iterator begin, func_iterator end,
                             const Boundary &b, int num_labels) {
    VectorMinkowskiFunctional *w110 = create_w110 ();
    calculate_functional (w110, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w110->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
}

static void set_refvert_coc (func_iterator begin, func_iterator end,
                             const Boundary &b, int num_labels) {
    VectorMinkowskiFunctional *w210 = create_w210 ();
    calculate_functional (w210, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w210->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
}

static void set_refvert_origin (func_iterator begin, func_iterator end,
                                int num_labels) {
    vec_t refvert = vec_t (0., 0.);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            (*begin)->ref_vertex (l, refvert);
        }
    }
}

static void set_refvert_domain_center (func_iterator begin, func_iterator end,
                                       const rect_t &r,
                                       int xdomains, int ydomains) {
    for (; begin != end; ++begin) {
        for (int l = 0; l != xdomains*ydomains; ++l) {
            (*begin)->ref_vertex (l, label_domain_center (
                l, r, xdomains, ydomains));
        }
    }
}

// the gigantic main function of the program.
int main (int argc, char **argv) {
    GetOpt_pp ops (argc, argv);

    // preparing the config file, and the command line.
    std::string configfile = basename (argv[0]) + ".conf";
    if (ops >> OptionPresent ('c', "--config")) {
        ops >> Option ('c', "--config", configfile);
    }
    std::cerr << "Using config file " << configfile << "\n";
    Configuration conf (configfile);

    std::string filename = conf.string ("input", "filename");
    if (ops >> OptionPresent ('i', "--input")) {
        // override the input specified in config file
        ops >> Option ('i', "--input", filename);
    }

    std::string output_prefix = "out/";
    if (ops >> OptionPresent ('o', "--output")) {
        // override the output specified in config file
        ops >> Option ('o', "--output", output_prefix);
    }

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
        Pixmap p;
        load_pgm (&p, filename);
        if (conf.boolean ("segment", "invert"))
            invert (&p);
        int  threshold    = conf.integer ("segment", "threshold");
        bool connectblack = conf.boolean ("segment", "connectblack");
        marching_squares (&b, p, threshold, connectblack);
    }

    std::string contfile (output_prefix + "contours.out");
    dump_contours (contfile, b, 1);

    ScalarMinkowskiFunctional *w000 = create_w000 ();
    ScalarMinkowskiFunctional *w100 = create_w100 ();
    ScalarMinkowskiFunctional *w200 = create_w200 ();
    VectorMinkowskiFunctional *w010 = create_w010 ();
    VectorMinkowskiFunctional *w110 = create_w110 ();
    VectorMinkowskiFunctional *w210 = create_w210 ();
    MatrixMinkowskiFunctional *w020 = create_w020 ();
    MatrixMinkowskiFunctional *w120 = create_w120 ();
    MatrixMinkowskiFunctional *w220 = create_w220 ();
    MatrixMinkowskiFunctional *w211 = create_w211 ();

    AbstractMinkowskiFunctional *all_funcs[]
        = { w000, w100, w200, w020, w120, w220, w211, w010, w110, w210 };
    func_iterator all_funcs_begin = all_funcs;
    func_iterator all_funcs_end   = all_funcs
                                    + sizeof (all_funcs)/sizeof (*all_funcs);

    std::string labcrit = conf.string ("output", "labels");
    std::string point_of_ref = conf.string ("output", "point_of_reference");
    int num_labels;
    if (labcrit == "none")
        num_labels = label_none (&b);
    else if (labcrit == "by_contour") {
        num_labels = label_by_contour_index (&b);
        if (point_of_ref == "contour_com")
            set_refvert_com (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "contour_cos")
            set_refvert_cos (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "contour_coc")
            set_refvert_coc (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "origin")
            set_refvert_origin (all_funcs_begin, all_funcs_end, num_labels);
        else
            die ("option \"point_of_reference\" in section [output] has illegal value");
    } else if (labcrit == "by_component") {
        num_labels = label_by_component (&b);
        if (point_of_ref == "component_com")
            set_refvert_com (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "component_cos")
            set_refvert_cos (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "component_coc")
            set_refvert_coc (all_funcs_begin, all_funcs_end, b, num_labels);
        else if (point_of_ref == "origin")
            set_refvert_origin (all_funcs_begin, all_funcs_end, num_labels);
        else
            die ("option \"point_of_reference\" in section [output] has illegal value");
    } else if (labcrit == "by_domain") {
        rect_t r;
        r.top    = conf.floating ("domains", "clip_top");
        r.right  = conf.floating ("domains", "clip_right");
        r.bottom = conf.floating ("domains", "clip_bottom");
        r.left   = conf.floating ("domains", "clip_left");
        int xdomains = conf.integer ("domains", "xdomains");
        int ydomains = conf.integer ("domains", "ydomains");
        num_labels =
            label_by_domain (&b, r, xdomains, ydomains);
        if (point_of_ref == "origin")
            set_refvert_origin (all_funcs_begin, all_funcs_end, num_labels);
        else if (point_of_ref == "domain_center")
            set_refvert_domain_center (all_funcs_begin, all_funcs_end, r, xdomains, ydomains);
        else 
            die ("option \"point_of_reference\" in section [output] has illegal value");
    } else {
        die ("option \"labels\" in section [output] has illegal value");
    }

    dump_labels (output_prefix + "labels", b);

    {
        // calculate all functionals
        func_iterator it;

        for (it = all_funcs_begin; it != all_funcs_end; ++it) {
            calculate_functional (*it, b);
        }
    }

    int precision = conf.integer ("output", "precision");

    {
        // output scalars
        ScalarMinkowskiFunctional *all_sca_begin[] = { w000, w100, w200 };
        ScalarMinkowskiFunctional **all_sca_end = all_sca_begin + 3;
        ScalarMinkowskiFunctional **it;
        std::string filename = output_prefix + "scalar" + ".out";
        std::ofstream of (filename.c_str ());
        of << std::setw (20) << "#   1          label";
        int col = 2;
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w000";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w100";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w200";
        of << "\n";

        for (int l = 0; l != num_labels; ++l) {
            of << " " << std::setw (19) << l;
            for (it = all_sca_begin; it != all_sca_end; ++it) {
                ScalarMinkowskiFunctional *p = *it;
                of << " " << std::setw (19) << std::setprecision (precision) << p->value (l);
            }
            of << "\n";
        }
    }

    {
        // output vectors
        VectorMinkowskiFunctional *all_sca_begin[] = { w010, w110, w210 };
        VectorMinkowskiFunctional **all_sca_end = all_sca_begin + 3;
        VectorMinkowskiFunctional **it;
        std::string filename = output_prefix + "vector" + ".out";
        std::ofstream of (filename.c_str ());
        of << std::setw (20) << "#   1          label";
        int col = 2;
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w010.x";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w010.y";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w110.x";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w110.y";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w210.x";
        of << std::setw ( 4) << col++;
        of << std::setw (16) << "w210.y";
        of << "\n";

        for (int l = 0; l != num_labels; ++l) {
            of << " " << std::setw (19) << l;
            for (it = all_sca_begin; it != all_sca_end; ++it) {
                vec_t val = (*it)->value (l);
                of << " " << std::setw (19) << std::setprecision (precision) << val.x ()
                   << " " << std::setw (19) << std::setprecision (precision) << val.y ();
            }
            of << "\n";
        }
    }

    {
        // output tensors
        MatrixMinkowskiFunctional *all_mat_begin[] = { w020, w120, w220, w211 };
        MatrixMinkowskiFunctional **all_mat_end = all_mat_begin + 4;
        MatrixMinkowskiFunctional **it;
        for (it = all_mat_begin; it != all_mat_end; ++it) {
            MatrixMinkowskiFunctional *p = *it;
            std::string filename = output_prefix + "tensor_" + p->name () + ".out";
            std::ofstream of (filename.c_str ());
            of << std::setw (20) << "#   1          label";
            int col = 2;
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "a11";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "a12";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "a21";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "a22";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "eval1";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "eval2";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "eval1/eval2";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec1";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec2";
            of << "\n";
            for (int l = 0; l != num_labels; ++l) {
                of << " " << std::setw (19) << l;
                mat_t val = p->value (l);
                of << " " << std::setw (19) << std::setprecision (precision) << val(0,0);
                of << " " << std::setw (19) << std::setprecision (precision) << val(0,1);
                of << " " << std::setw (19) << std::setprecision (precision) << val(1,0);
                of << " " << std::setw (19) << std::setprecision (precision) << val(1,1);
                EigenSystem esys;
                eigensystem (&esys, val);
                if (fabs (esys.eval[0]) < fabs (esys.eval[1])) {
                    swap_eigenvalues (&esys);
                }
                double ratio = esys.eval[1]/esys.eval[0];
                of << " " << std::setw (19) << std::setprecision (precision) << esys.eval[0];
                of << " " << std::setw (19) << std::setprecision (precision) << esys.eval[1];
                of << " " << std::setw (19) << std::setprecision (precision) << ratio;
                of << " " << std::setw (19) << std::setprecision (precision) << esys.evec[0][0];
                of << " " << std::setw (19) << std::setprecision (precision) << esys.evec[0][1];
                of << " " << std::setw (19) << std::setprecision (precision) << esys.evec[1][0];
                of << " " << std::setw (19) << std::setprecision (precision) << esys.evec[1][1];
                of << "\n";
            }
        }
    }

    delete w000;
    delete w100;
    delete w200;
    delete w010;
    delete w110;
    delete w210;
    delete w020;
    delete w120;
    delete w220;
    delete w211;


    return 0;
}
