// read a pgm graphics file into a Pixmap object.
// reads 16-bit and 8-bit PGM.

#include "util.h"
#include <fstream>
#include <stdexcept>
using namespace std;

static void format_error (const string &msg) {
    throw std::runtime_error (msg);
}

static long read_header (Pixmap *p, string *magic, std::string *comment, istream &is) {
    is >> *magic;
    if (*magic != "P1" && *magic != "P2" && *magic != "P4" &&
            *magic != "P5")
        format_error ("magic incorrect");
    is >> ws;
    if (is.peek () == '#')
        getline (is, *comment); // ignore comment
    int w, h;
    long max_value = 1;
    is >> w >> h;
    if (*magic == "P2" || *magic == "P5")
        is >> max_value;
    is >> ws;
    // allocate pixmap
    if (w <= 0 || h <= 0)
        format_error ("header damaged");
    p->resize (w, h);
    return max_value;
}

void load_pgm (Pixmap *p, const string &filename) {
    assert (p);
    ifstream is (filename.c_str ());
    string magic, comment;
    if (!is)
        throw std::runtime_error ("Cannot open \"" + filename + "\"");
    is.exceptions (ios::failbit | ios::badbit);
    long max_value = read_header (p, &magic, &comment, is);
    double nrml = 1. / max_value;

    if (magic == "P2") {   // PGM ASCII
        for (int j = 0; j != p->size2 (); ++j)
        for (int i = 0; i != p->size1 (); ++i) {
            long tmp;
            is >> tmp;
            (*p)(i,j) = Pixmap::val_t (tmp * nrml);
        }
    }
    else if (magic == "P1") // PBM ASCII
    {
        for (int j = 0; j != p->size2 (); ++j)
        for (int i = 0; i != p->size1 (); ++i) {
            int tmp;
            is >> tmp;
            (*p)(i,j) = Pixmap::val_t (1-tmp); // PBM: black=1, white=0
        }
    }
    else if (magic == "P4") // PBM binary
    {
        for (int j = 0; j != p->size2 (); ++j)
        for (int i = 0; i != p->size1 (); ) {
            char tmp;
            is.get(tmp); 
            // 1 byte=8 pixel, split accordingly
            unsigned char mask = 0x80;
            while (mask) {
                // extract bit
                int value = !! (static_cast <unsigned char> (tmp) & mask);
                (*p)(i,j) = Pixmap::val_t (1-value); // black=1, white=0
                if (++i == p->size1())
                    break;
                mask >>= 1;
            }
        }
    }
    else if (magic == "P5") // PGM binary
    {
        // <256 == 1byte/pix
        if (max_value<256)
            for (int j = 0; j != p->size2 (); ++j)
            for (int i = 0; i != p->size1 (); ++i) {
                char tmp;
                is.get(tmp);
                const long value=(unsigned char)(tmp);
                (*p)(i,j) = Pixmap::val_t (value * nrml);
            }
        // >255 == 2byte/pix
        else
            for (int j = 0; j != p->size2 (); ++j)
            for (int i = 0; i != p->size1 (); ++i) {
                char tmp1, tmp2;
                is.get(tmp1);
                is.get(tmp2);
                long tmp = long (static_cast <unsigned char> (tmp1)) << 8 |
                           long (static_cast <unsigned char> (tmp2));
                (*p)(i,j) = Pixmap::val_t (tmp * nrml);
            }
    }
    else
        format_error ("magic incorrect");

    // file should be empty now
    is.exceptions (ios::badbit); // Clear eofbit before reading beyond limit
    is >> ws;
    is.get (); // try to read something
    if (is)
        format_error ("trailing data in \'" + filename ="\"");
    else
        return; // yup, it's empty
}
