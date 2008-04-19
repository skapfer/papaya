// vim: et:sw=4:ts=4
// read some data out of POLY files

#include "util.h"
#include <fstream>
#include <stdexcept>
#include <limits>
#include <stdarg.h>
#include <stdio.h>
typedef std::string string;

static void format_error (const string &message, ...) {
    char buffer[200];
    va_list va;
    va_start (va, message);
    vsnprintf (buffer, 199, message.c_str (), va);
    throw std::runtime_error (buffer);
}

// object keeping the state during the process of reading a POLY file
class PolyFileReader {
public:
    std::ifstream is;
    Boundary *b;
    // mapping POLY vertex indices -> Boundary vertex indices
    std::vector <int> vertex_map;

    // isdigit independant of locale
    static bool is_digit (int c) {
        switch (c) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
            return true;
        default:
            return false;
        }
    }

    PolyFileReader (Boundary *b_, const string &polyfilename)
        : is (polyfilename.c_str ()),
          b (b_)
    {
        if (!is) {
            throw std::runtime_error ("Cannot open \"" +
                polyfilename + "\"");
        }
        is.exceptions (std::ios::failbit | std::ios::badbit);
        vertex_map.reserve (1000);
    }

    void extract_header (const string &tag) {
        std::string tmp;
        std::getline (is, tmp);
        if (tmp != tag) {
            format_error ("format error in POLY file: expected %s, got %s",
                tag.c_str (), tmp.c_str ());
        }
    }

    void ignore_rest_of_line () {
        std::streamsize mx = std::numeric_limits <std::streamsize>::max ();
        is.ignore (mx, '\n');
    }

    void read_vertex () {
        int number;
        is >> number;

        // vertex number is terminated by colon
        is >> std::ws;
        int colon;
        if ((colon = is.get ()) != ':')
            format_error ("expected colon, got %c", colon);

        // coords
        double x, y, z;
        is >> x >> y >> z;

        // rest of record (i.e. attributes) is ignored for now
        ignore_rest_of_line ();

        // z coordinate is ignored
        int vert_id = b->insert_vertex (vec_t (x, y));

        // saved for later reference
        if ((int)vertex_map.size () <= number) {
            vertex_map.resize (number + 1, Boundary::INVALID_VERTEX);
        }
        vertex_map[number] = vert_id;
    }

    int lookup_vertex (int v_id) {
        if ((int)vertex_map.size () > v_id && v_id >= 0) {
            int ret = vertex_map[v_id];
            if (ret != Boundary::INVALID_VERTEX)
                return ret;
        }
        format_error ("file refers to vertex %i which I do not know (vertex_map size = %i)",
                      v_id, (int)vertex_map.size ());
    }

    void read_poly () {
        int number;
        is >> number;

        // poly number is terminated by colon
        is >> std::ws;
        int colon;
        if ((colon = is.get ()) != ':')
            format_error ("PolyFileReader::read_poly: expected colon, got %c", colon);

        // first two vertices
        int initial_vertex;
        is >> initial_vertex;
        // translate POLY vertices into Boundary vertices
        initial_vertex = lookup_vertex (initial_vertex);
        int prev_vertex;
        is >> prev_vertex >> std::ws;
        prev_vertex = lookup_vertex (prev_vertex);

        int initial_edge = b->insert_edge (Boundary::INVALID_EDGE,
                                    initial_vertex, prev_vertex,
                                    Boundary::INVALID_EDGE);

        int prev_edge = initial_edge;
        while (is_digit (is.peek ())) {
            int vertex;
            is >> vertex >> std::ws;
            // translate POLY vertices into Boundary vertices
            vertex = lookup_vertex (vertex);

            int edge = b->insert_edge (prev_edge, prev_vertex,
                                       vertex, Boundary::INVALID_EDGE);
            prev_vertex = vertex;
            prev_edge = edge;
        }

        is >> std::ws;
        int closed_indic = is.peek ();
        if (closed_indic == '<') {
            // extract '<' character
            is.get ();
            // close contour
            b->insert_edge (prev_edge, prev_vertex, initial_vertex, initial_edge);
        } else {
            std::cerr << "WARNING: non-closed contour " << number << ". this is probably not what you want.\n";
        }

        ignore_rest_of_line ();
    }
};

void load_poly (Boundary *b, const string &polyfilename) {
    PolyFileReader r (b, polyfilename);
    r.extract_header ("POINTS");
    while (int p = r.is.peek ()) {
        if (r.is_digit (p)) {
            // incoming number
            r.read_vertex ();
        } else {
            // something else
            goto no_more_points;
        }
    }
no_more_points:
    r.extract_header ("POLYS");
    while (int p = r.is.peek ()) {
        if (r.is_digit (p)) {
            // incoming number
            r.read_poly ();
        } else {
            // something else
            goto no_more_polys;
        }
    }
no_more_polys:
    r.extract_header ("END");
    // rest of file is ignored.
    r.is.close ();
}

