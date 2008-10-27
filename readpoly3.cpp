// vim: et:sw=4:ts=4
// read some data out of .poly files

#include "readpoly.h"
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace {

    // object keeping the state during the process of reading a .poly file
    struct PolyFileReader {
        PolyFileReader (std::istream &is_, PolyFileSink *sink_)
            : is (is_), sink (sink_)
        {
            if (is) {
                throw_error ("Stream is not readable.");
            }
        }

        std::istream &is;
        PolyFileSink *sink;
        // mapping .poly vertex indices -> target vertex indices
        std::vector <int> vertex_map;
        typedef PolyFileSink::prop_list_t prop_list_t;

        static void throw_error (const std::string &message, ...) {
            char buffer[200];
            va_list va;
            va_start (va, message);
            vsnprintf (buffer, 199, message.c_str (), va);
            throw std::runtime_error (buffer);
        }

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

        static std::string printable_int (int c) {
            if (c == EOF)
                return "EOF";
            else
                return std::string (1, char (c));
        }

        void extract_header (const std::string &tag) {
            std::string tmp;
            std::getline (is, tmp);
            if (tmp != tag) {
                throw_error ("format error in .poly file: expected %s, got %s",
                    tag.c_str (), tmp.c_str ());
            }
        }

        void ignore_rest_of_line () {
            std::streamsize mx = std::numeric_limits <std::streamsize>::max ();
            is.ignore (mx, '\n');
        }

        int lookup_vertex (int v_id) {
            if ((int)vertex_map.size () > v_id && v_id >= 0) {
                int ret = vertex_map[v_id];
                if (ret >= 0)
                    return ret;
            }
            throw_error ("file refers to vertex %i which I do not know (vertex_map size = %i)",
                         v_id, (int)vertex_map.size ());
            abort ();  // silence warning
        }

        void read_vertex () {
            int number;
            is >> number;

            // vertex number is terminated by colon
            int colon = (is >> std::ws).get ();
            if (colon != ':')
                throw_error ("Expected colon, got %s.", printable_int (colon).c_str ());

            // coords
            double x, y, z;
            is >> x >> y >> z;
            if (!is)
                throw_error ("Cannot read coordinates.");

            // FIXME read props
            ignore_rest_of_line ();

            int vert_id = sink->insert_vertex (x, y, z, prop_list_t ());

            // saved for later reference
            if ((int)vertex_map.size () <= number) {
                vertex_map.resize (number+1, -1);
            }
            vertex_map[number] = vert_id;
        }

        void read_facet () {
            std::string line;
            std::getline (is, line);
            std::istringstream iss (line);

            int number;
            iss >> number;

            // vertex number is terminated by colon
            int colon = (iss >> std::ws).get ();
            if (colon != ':')
                throw_error ("Expected colon, got %s.", printable_int (colon).c_str ());

            iss >> std::ws;
            std::vector <int> vl;
            vl.reserve (10);
            while (is_digit (iss.peek ())) {
                int v;
                iss >> v >> std::ws;
                v = lookup_vertex (v);
                vl.push_back (v);
            }

            prop_list_t pl;
            std::string word;
            while (iss >> word)
                pl.push_back (word);

            sink->insert_facet (vl, pl);
        }

        void read_file () {
            extract_header ("POINTS");
            while (is_digit (is.peek ()))
                read_vertex ();
            extract_header ("POLYS");
            while (is_digit (is.peek ()))
                read_facet ();
            extract_header ("END");
            // rest of file is ignored.
        }
    };

}

void parse_poly_file (PolyFileSink *sink, std::istream &is) {
    PolyFileReader r (is, sink);
    r.read_file ();
}
