
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <limits>
#include <math.h>
#include <assert.h>
#include <string>
#include "util.h"

class TsvFile {
public:
    void read_file (const std::string &);
    int num_rows () const;
    int num_cols (int) const;
    double operator() (int i, int j) const;

    struct unreadable_file : public std::runtime_error {
        unreadable_file (const std::string &filename)
            : std::runtime_error ("unable to read \"" + filename + "\"") {}
    };

private:
    std::vector <std::vector <double> > my_data;
    static void eat_line (std::istream &is) {
        for (;;) {
            int c = is.peek ();
            if (c == EOF || c == '\n')
                return;
            else
                is.get ();
        }
    }
    static void eat_spaces (std::istream &is) {
        for (;;) {
            int c = is.peek ();
            if (c == ' ' || c == '\t')
                is.get ();
            else
                return;
        }
    }
};

void TsvFile::read_file (const std::string &filename) {
    std::ifstream is (filename.c_str ());
    if (!is)
        throw unreadable_file (filename);
    // iostreams exceptions just suck.
    //is.exceptions (std::ios::badbit);
    std::vector <double> line;
    for (;;) {
        eat_spaces (is);
        int pk = is.peek ();
        switch (pk) {
        case '\n':
            is.get ();
            if (line.size ()) {
                my_data.push_back (line);
                line.clear ();
            }
            continue;
        case EOF:
            if (line.size ()) {
                my_data.push_back (line);
                line.clear ();
            }
            return;
        case '#':
            eat_line (is);
            continue;
        default:
            double tmp;
            is >> tmp;
            line.push_back (tmp);
        }
    }
}

int TsvFile::num_rows () const {
    return my_data.size ();
}

int TsvFile::num_cols (int i) const {
    return my_data.at (i).size ();
}

double TsvFile::operator() (int i, int j) const {
    return my_data.at (i).at (j);
}

int main (int argc, const char **argv) { 
    const std::string tolerance_o = "--tolerance";
    double tolerance = 1e-3;
    std::string filename1, filename2;
    for (++argv; *argv; ++argv) {
        if (tolerance_o == *argv) {
            ++argv;
            if (!*argv) {
                die ("Missing tolerance");
            } else {
                tolerance = strtod (*argv, NULL);
            }
        } else {
            if (!filename1.size ()) {
                filename1 = *argv;
            } else if (!filename2.size ()) {
                filename2 = *argv;
            } else {
                die ("Too many filenames: %s", *argv);
            }
        }
    }

    TsvFile f1, f2;
    bool emitted_filenames = false;
    f1.read_file (filename1);
    f2.read_file (filename2);
    if (f1.num_rows () != f2.num_rows ()) {
        fprintf (stderr, "[tsvdiff] Number of rows don't match"
                 " comparing %s and %s\n",
                 filename1.c_str (), filename2.c_str ());
        return 1;
    }
    for (int i = 0; i != f1.num_rows (); ++i) {
        if (f1.num_cols (i) != f2.num_cols (i)) {
            fprintf (stderr, "[tsvdiff] Number of columns in row %i"
                     " don't match comparing %s and %s\n", i,
                     filename1.c_str (), filename2.c_str ());
            return 1;
        }
        for (int j = 0; j != f1.num_cols (i); ++j) {
            double delta = fabs (f1(i,j) - f2(i,j));
            if (delta > tolerance) {
                if (!emitted_filenames) {
                    std::cout << filename1 << " " << filename2 << "\n";
                    emitted_filenames = true;
                }
                std::cout << " " << std::setw (4) << i
                          << " " << std::setw (4) << j
                          << " " << std::setw (20) << std::setprecision (15) << f1(i,j)
                          << " " << std::setw (20) << std::setprecision (15) << f2(i,j)
                          << "\n";
            }
        }
    }
    return 0;
}
