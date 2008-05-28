// read a pgm graphics file into a Pixmap object.
// reads 16-bit and 8-bit PGM.

#include "util.h"
#include <fstream>
#include <stdexcept>
using namespace std;

static void format_error (const string &msg) {
    throw std::runtime_error (msg);
}

static long read_header (Pixmap *p, istream &is) {
    string magic;
    is >> magic;
    if (magic != "P2")
        format_error ("magic incorrect");
    is >> ws;
    if (is.peek () == '#')
        getline (is, magic); // ignore filename
    int w, h;
    long max_value;
    is >> w >> h >> max_value >> ws;
    // allocate pixmap
    if (w <= 0 || h <= 0)
        format_error ("header damaged");
    p->resize (w, h);
    return max_value;
}

void load_pgm (Pixmap *p, const string &filename) {
    assert (p);
    ifstream is (filename.c_str ());
    if (!is)
        throw std::runtime_error ("Cannot open \"" + filename + "\"");
    is.exceptions (ios::failbit | ios::badbit);
    double nrml = 1. / read_header (p, is);
    for (int j = 0; j != p->size2 (); ++j)
    for (int i = 0; i != p->size1 (); ++i) {
        long tmp;
        is >> tmp;
        (*p)(i,j) = Pixmap::val_t (tmp * nrml);
    }

    // file should be empty now
    is >> ws;
    is.exceptions (ios::badbit);
    is.get (); // try to read something
    if (is)
        format_error ("trailing data in .pgm");
    else
        return; // yup, it's empty
}
