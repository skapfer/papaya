#ifndef READ_POLY_H_INCLUDED 
#define READ_POLY_H_INCLUDED 
    
#include <string>
#include <vector>
#include <istream>

class PolyFileSink {
public:
    typedef std::vector <int> vertex_id_list_t;
    typedef std::vector <std::string> prop_list_t;

    virtual ~PolyFileSink () {}
    virtual int  insert_vertex (double, double, double,
                                const prop_list_t &) = 0;
    virtual void insert_facet  (const vertex_id_list_t &,
                                const prop_list_t &) = 0;
};

void parse_poly_file (PolyFileSink *sink, std::istream &);

#endif // READ_POLY_H_INCLUDED
