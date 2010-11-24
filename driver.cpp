// vim: et:sw=4:ts=4
// driver code, the Papaya command-line tool.
#include <iostream>
#include <fstream>
#include <getopt_pp_standalone.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <errno.h>
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

static
std::string detect_fileformat (std::string filename)
{
    if (ends_with (filename, ".poly")) {
        return "poly";
    } else if (ends_with (filename, ".pgm") || ends_with (filename, ".pbm")) {
        return "pgm";
    } else {
        die ("Cannot deduce the input file format from filename %s.  "
            "Please specify --format=poly or --format=pgm.",
            filename.c_str ());
    }
}

static
std::string my_basename (const std::string &str) {
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

typedef std::vector <std::string> string_vector;

static
bool vector_contains (const string_vector &v, const std::string &value)
{
    return std::find (v.begin (), v.end (), value) != v.end ();
}

static
string_vector split_at_comma (std::string input)
{
    string_vector ret;
    std::string::size_type x = 0u, y = input.find_first_of (',');
    for(;;)
    {
        if (y == input.npos) {
            ret.push_back (input.substr (x));
            return ret;
        } else {
            ret.push_back (input.substr (x, y-x));
            x = y+1;
            y = input.find_first_of (',', x);
        }
    }
}

static
string_vector parse_what_to_compute (std::string what)
{
    string_vector what_to_c = split_at_comma (what);
    const char *const tensors[] = { "W020", "W120", "W220", "W211", "W102" };
    const int num_tensors = sizeof(tensors)/sizeof(*tensors);
    string_vector legal_options (tensors, tensors+num_tensors);
    legal_options.push_back ("contours");
    legal_options.push_back ("labels");
    legal_options.push_back ("scalars");
    legal_options.push_back ("vectors");
    legal_options.push_back ("tensors");

    // check that only legal options are given
    string_vector::iterator it;
    for (it = what_to_c.begin (); it != what_to_c.end (); ++it)
    {
        std::cerr << *it << " ";
        if (!vector_contains (legal_options, *it))
        {
            std::cerr << "invalid value in \"compute\" option: " << *it << std::endl;
            std::abort ();
        }
    }
    std::cerr << "\n";

    // if "tensors" is given, expand this to /all/ the tensors we can compute
    if (vector_contains (what_to_c, "tensors"))
        what_to_c.insert (what_to_c.end (), tensors, tensors+num_tensors);

    return what_to_c;
}

// recursively create directories for the output files
static
void prefix_mkdir (std::string prefix)
{
    std::string::size_type p = prefix.find_last_of ("/");
    if (p == prefix.npos) return;
    prefix.assign (prefix, 0u, p);
    prefix_mkdir (prefix);
    errno = 0;
    std::cerr << "[papaya] mkdir " << prefix << "\n";
    mkdir (prefix.c_str (), 0777);
    if (errno && errno != EEXIST)
        perror ("mkdir");
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
    ScalarMinkowskiFunctional *w000 = create_w000 ();
    w010->global_ref_vertex (vec_t (0., 0.));
    calculate_functional (w010, b);
    calculate_functional (w000, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w010->value (l) / w000->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
    delete w010;
    delete w000;
}

static void set_refvert_cos (func_iterator begin, func_iterator end,
                             const Boundary &b, int num_labels) {
    VectorMinkowskiFunctional *w110 = create_w110 ();
    ScalarMinkowskiFunctional *w100 = create_w100 ();
    w110->global_ref_vertex (vec_t (0., 0.));
    calculate_functional (w110, b);
    calculate_functional (w100, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w110->value (l) / w100->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
    delete w110;
    delete w100;
}

static void set_refvert_coc (func_iterator begin, func_iterator end,
                             const Boundary &b, int num_labels) {
    VectorMinkowskiFunctional *w210 = create_w210 ();
    ScalarMinkowskiFunctional *w200 = create_w200 ();
    w210->global_ref_vertex (vec_t (0., 0.));
    calculate_functional (w210, b);
    calculate_functional (w200, b);
    for (; begin != end; ++begin) {
        for (int l = 0; l != num_labels; ++l) {
            vec_t refvert = w210->value (l) / w200->value (l);
            (*begin)->ref_vertex (l, refvert);
        }
    }
    delete w210;
    delete w200;
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

static void save_ref_vertex_map (AbstractMinkowskiFunctional  *w020,
                                 std::string output_prefix,
                                 int precision,
                                 int xdomains, int ydomains) {
    std::string filename = output_prefix + "by_domain_ref_vertex.out";
    std::ofstream of (filename.c_str ());
    if (!of)
        std::cerr << "[papaya] WARNING unable to open " << filename << "\n";
    print_version_header (of);

    of << std::setw (20) << "#   1          label";
    int col = 2;
    of << std::setw ( 4) << col++;
    of << std::setw (16) << "domain no x";
    of << std::setw ( 4) << col++;
    of << std::setw (16) << "domain no y";
    of << std::setw ( 4) << col++;
    of << std::setw (16) << "refvert x W020";
    of << std::setw ( 4) << col++;
    of << std::setw (16) << "refvert y W020";
    of << "\n";

    for (int l = 0; l != xdomains*ydomains; ++l) {
        of << " " << std::setw (19) << l;
        of << " " << std::setw (19) << (l%xdomains);
        of << " " << std::setw (19) << (l/xdomains);
        of << " " << std::setw (19) << std::setprecision (precision) << (w020->ref_vertex (l)[0]);
        of << " " << std::setw (19) << std::setprecision (precision) << (w020->ref_vertex (l)[1]);
        of << "\n";
    }
}


// the gigantic main function of the program.
int main (int argc, char **argv) {
    GetOpt_pp ops (argc, argv);

    if (ops >> OptionPresent ('v', "version")) {
        std::cerr << "papaya " << VERSION << "\n";
        return 0;
    }

    // preparing the config file, and the command line.
    std::string configfile = my_basename (argv[0]) + ".conf";
    if (ops >> OptionPresent ('c', "config")) {
        ops >> Option ('c', "config", configfile);
    }
    std::cerr << "[papaya] Using config file " << configfile << "\n";
    Configuration conf (configfile);

    std::string filename = conf.string ("input", "filename");
    std::string in_fileformat = conf.string ("input", "format", "deduce_from_filename");
    if (ops >> OptionPresent ('i', "input"))
        // override the input specified in config file
        ops >> Option ('i', "input", filename);
    if (ops >> OptionPresent ('F', "format"))
        ops >> Option ('F', "format", in_fileformat);
    if (in_fileformat == "deduce_from_filename")
        in_fileformat = detect_fileformat (filename);

    std::cerr << "[papaya] Using input file " << filename << "\n";

    std::string output_prefix = conf.string ("output", "prefix");
    if (ops >> OptionPresent ('o', "output")) {
        // override the output specified in config file
        ops >> Option ('o', "output", output_prefix);
    }
    std::cerr << "[papaya] Using output prefix " << output_prefix << "\n";
    if (conf.boolean ("output", "mkdir", false))
        prefix_mkdir (output_prefix);

    double thresh_override = -INFINITY;
    if (ops >> OptionPresent ('\0', "threshold")) {
        ops >> Option ('\0', "threshold", thresh_override);
        std::cerr << "[papaya] Using threshold " << thresh_override << "\n";
    }

    int precision = conf.integer ("output", "precision");

    std::string normalization = conf.string ("output", "normalization", "code_default");
    if (ops >> OptionPresent ('N', "normalization"))
        ops >> Option ('N', "normalization", normalization);
    if (normalization == "breidenbach") {
        W0_NORMALIZATION = W1_NORMALIZATION = W2_NORMALIZATION = 1.;
    } else if (normalization == "new") {
        W0_NORMALIZATION = 1.;
        W1_NORMALIZATION = W2_NORMALIZATION = .5;
    } else if (normalization == "code_default") {
    } else {
        std::cerr << "Invalid normalization setting: " << normalization << std::endl;
        return 1;
    }

    Boundary b, b_for_w0_storage_;
    Boundary *b_for_w0 = &b;

    if (in_fileformat == "poly") {
        load_poly (&b, filename);
        if (thresh_override != -INFINITY) {
            std::cerr << "--threshold is not useful in .poly mode.\n";
            abort ();
        }
        bool runfix   = conf.boolean ("polyinput", "fix_contours");
        bool forceccw = conf.boolean ("polyinput", "force_counterclockwise");
        if (runfix)
            fix_contours (&b, conf.boolean ("polyinput", "silent_fix_contours"));
        if (forceccw)
            force_counterclockwise_contours (&b);
    } else if (in_fileformat == "pgm" || in_fileformat == "pbm") {
        Pixmap p;
        load_pgm (&p, filename);
        if (conf.boolean ("segment", "invert"))
            invert (&p);
        double threshold  = conf.floating ("segment", "threshold");
        if (thresh_override != -INFINITY)
            threshold = thresh_override;
        bool connectblack = conf.boolean ("segment", "connectblack");
        bool periodic_data = conf.boolean ("segment", "data_is_periodic");
        marching_squares (&b, p, threshold, connectblack, periodic_data);
    } else {
        std::cerr << "only .pgm, .pbm and .poly files are valid input"
                  << "(\"" <<  filename << "\")" << std::endl;
        abort ();
    }

    assert_sensible_boundary (b);

    // find out what we're supposed to compute
    std::string default_what = "contours,labels,scalars,vectors,tensors";
    string_vector what_to_compute =
        parse_what_to_compute (conf.string ("output", "compute", default_what));
    if (ops >> OptionPresent (' ', "compute")) {
        std::string what;
        ops >> Option (' ', "compute", what);
        what_to_compute = parse_what_to_compute (what);
    }

    // write contours prior to labelling (in case that crashes...)
    if (vector_contains (what_to_compute, "contours")) {
        std::string contfile (output_prefix + "contours");
        dump_contours (contfile, b, 1);
    }

    ScalarMinkowskiFunctional *w000 = create_w000 ();
    ScalarMinkowskiFunctional *w100 = create_w100 ();
    ScalarMinkowskiFunctional *w200 = create_w200 ();
    VectorMinkowskiFunctional *w010 = create_w010 ();
    VectorMinkowskiFunctional *w110 = create_w110 ();
    VectorMinkowskiFunctional *w210 = create_w210 ();
    MatrixMinkowskiFunctional *w020 = create_w020 ();
    MatrixMinkowskiFunctional *w120 = create_w120 ();
    MatrixMinkowskiFunctional *w102 = create_w102 ();
    MatrixMinkowskiFunctional *w220 = create_w220 ();
    MatrixMinkowskiFunctional *w211 = create_w211 ();

    AbstractMinkowskiFunctional *all_funcs[]
        = { w000, w100, w200, w020, w120, w102, w220, w211, w010, w110, w210 };
    func_iterator all_funcs_begin = all_funcs;
    func_iterator all_funcs_end   = all_funcs
                                    + sizeof (all_funcs)/sizeof (*all_funcs);

    std::string labcrit = conf.string ("output", "labels");
    std::string point_of_ref = conf.string ("output", "point_of_reference");
    int num_labels = -1;
    if (labcrit == "none") {
        num_labels = label_none (&b);
        set_refvert_origin (all_funcs_begin, all_funcs_end, num_labels);
    } else if (labcrit == "by_contour") {
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
        b_for_w0_storage_ = b;
        b_for_w0 = &b_for_w0_storage_;
        num_labels = label_by_domain (&b, r, xdomains, ydomains, false);
                     label_by_domain (b_for_w0, r, xdomains, ydomains, true);
        if (point_of_ref == "origin")
            set_refvert_origin (all_funcs_begin, all_funcs_end, num_labels);
        else if (point_of_ref == "domain_center")
            set_refvert_domain_center (all_funcs_begin, all_funcs_end, r, xdomains, ydomains);
        else 
            die ("option \"point_of_reference\" in section [output] has illegal value");

        save_ref_vertex_map (w020, output_prefix, precision, xdomains, ydomains);
    } else {
        die ("option \"labels\" in section [output] has illegal value");
    }

    assert (num_labels != -1);
    if (vector_contains (what_to_compute, "labels")) {
        dump_labels (output_prefix + "labels", b);
        dump_labels (output_prefix + "nu0labels", *b_for_w0);
    }

    {
        // calculate all functionals
        func_iterator it;

        for (it = all_funcs_begin; it != all_funcs_end; ++it) {
            if (*it == w000 || *it == w010 || *it == w020)
                calculate_functional (*it, *b_for_w0);
            else
                calculate_functional (*it, b);
        }
    }

    if (vector_contains (what_to_compute, "scalars"))
    {
        // output scalars
        ScalarMinkowskiFunctional *all_sca_begin[] = { w000, w100, w200 };
        ScalarMinkowskiFunctional **all_sca_end = all_sca_begin + 3;
        ScalarMinkowskiFunctional **it;
        std::string filename = output_prefix + "scalar.out";
        std::ofstream of (filename.c_str ());
        if (!of)
            std::cerr << "[papaya] WARNING unable to open " << filename << "\n";
        print_version_header (of);
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

    if (vector_contains (what_to_compute, "vectors"))
    {
        // output vectors
        VectorMinkowskiFunctional *all_sca_begin[] = { w010, w110, w210 };
        VectorMinkowskiFunctional **all_sca_end = all_sca_begin + 3;
        VectorMinkowskiFunctional **it;
        std::string filename = output_prefix + "vector.out";
        std::ofstream of (filename.c_str ());
        if (!of)
            std::cerr << "[papaya] WARNING unable to open " << filename << "\n";
        print_version_header (of);
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
        MatrixMinkowskiFunctional *all_mat_begin[] = { w020, w120, w102, w220, w211 };
        MatrixMinkowskiFunctional **all_mat_end = all_mat_begin + 5;
        MatrixMinkowskiFunctional **it;
        for (it = all_mat_begin; it != all_mat_end; ++it) {
            MatrixMinkowskiFunctional *p = *it;
            if (!vector_contains (what_to_compute, p->name ()))
                continue;
            std::string filename = output_prefix + "tensor_" + p->name () + ".out";
            std::ofstream of (filename.c_str ());
            if (!of)
                std::cerr << "[papaya] WARNING unable to open " << filename << "\n";
            print_version_header (of);
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
            of << std::setw (16) << "eval2/eval1";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec1x";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec1y";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec2x";
            of << std::setw ( 4) << col++;
            of << std::setw (16) << "evec2y";
            of << "\n";
            for (int l = 0; l != num_labels; ++l) {
                of << " " << std::setw (19) << l;
                mat_t val = p->value (l);
                of << " " << std::setw (19) << std::setprecision (precision) << val(0,0);
                of << " " << std::setw (19) << std::setprecision (precision) << val(0,1);
                of << " " << std::setw (19) << std::setprecision (precision) << val(1,0);
                of << " " << std::setw (19) << std::setprecision (precision) << val(1,1);
                EigenSystem esys;
                eigensystem_symm (&esys, val);
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
    delete w102;

    return 0;
}
