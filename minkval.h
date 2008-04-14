#ifndef MINKOWSKIVALUATIONS_H_INCLUDED 
#define MINKOWSKIVALUATIONS_H_INCLUDED 

#include "util.h"
#include <map>

class SurfaceMinkowskiValuator {
public:
    typedef Boundary::edge_iterator edge_iterator;
    typedef std::string string;

    SurfaceMinkowskiValuator ();
    void set_reference (vec_t);
    ~SurfaceMinkowskiValuator ();

    void add_contourseg (Boundary *b, edge_iterator begin, edge_iterator end);

    double scalar (const string &);
    vec_t  vector (const string &);
    mat_t  matrix (const string &);

    /*
    double W100 (const Boundary &);
    double W200 (const Boundary &);
    vec_t W110 (const Boundary &);
    vec_t W101 (const Boundary &);
    */

private:
    vec_t edge_normal (edge_iterator);
    vec_t edge_vertex0 (edge_iterator);
    vec_t edge_vertex1 (edge_iterator);
    vec_t my_ref;
    Boundary *my_b;
    std::map <string, vec_t> vacc;
    std::map <string, double> sacc;
    std::map <string, mat_t> macc;
};

#endif /* MINKOWSKIVALUATIONS_H_INCLUDED */
