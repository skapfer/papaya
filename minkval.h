#ifndef MINKOWSKIVALUATIONS_H_INCLUDED 
#define MINKOWSKIVALUATIONS_H_INCLUDED 

#include "util.h"
#include <map>


template <typename VALUE_TYPE>
class GenericMinkowskiFunctional {
public:
    typedef VALUE_TYPE value_t;
    typedef int label_t;
    enum { MAX_LABELS = 20000 };

public:
    GenericMinkowskiFunctional (const std::string &name);

    const std::string &name () const;

    void dump (std::ostream &) const;

protected:
    value_t &acc (label_t);
    void reszacc (label_t);
private:
    std::vector <value_t> my_acc;
    std::string my_name;
};

template <typename VALUE_TYPE>
std::ostream &operator<< (std::ostream &os, GenericMinkowskiFunctional <VALUE_TYPE> &f) {
    f.dump (os);
    return os;
}

typedef GenericMinkowskiFunctional <double> ScalarMinkowskiFunctional;
typedef GenericMinkowskiFunctional <vec_t>  VectorMinkowskiFunctional;
typedef GenericMinkowskiFunctional <mat_t>  MatrixMinkowskiFunctional;

class SurfaceIntegral {
public:
    typedef Boundary::edge_iterator edge_iterator;

    SurfaceIntegral ();

    virtual void add_contour (const Boundary &, edge_iterator begin, edge_iterator end) = 0;

    void ref_vertex (const vec_t &);

protected:
    const vec_t &ref_vertex () const;
private:
    vec_t my_ref_vertex;
};

class VolumeIntegral {
public:
    typedef Boundary::edge_iterator edge_iterator;
    // FIXME
};


void calculate_all_surface_integrals (const Boundary &b);

//
// inline implementation
//
template <typename VALUE_TYPE>
inline GenericMinkowskiFunctional<VALUE_TYPE>::GenericMinkowskiFunctional (const std::string &name)
    : my_name (name) { }

template <typename VALUE_TYPE>
inline VALUE_TYPE &GenericMinkowskiFunctional<VALUE_TYPE>::acc (int label) {
    if (int (my_acc.size ()) > int (label)) {
        return my_acc[int (label)];
    } else {
        reszacc (label);
        return this->acc (label);
    }
}

template <typename VALUE_TYPE>
inline const std::string &GenericMinkowskiFunctional <VALUE_TYPE>::name () const {
    return my_name;
}

// allocate new labels (i.e. accumulators for more labels)
template <typename VALUE_TYPE>
void GenericMinkowskiFunctional<VALUE_TYPE>::reszacc (label_t l) {
    if (l < MAX_LABELS) {
        int currsize = (int)my_acc.size ();
        my_acc.resize (int (l) + 1);
        memset (&*(my_acc.begin () + currsize), 0, sizeof (value_t) * (int (l) + 1 - currsize));
    } else {
        die ("GenericMinkowskiFunctional: too many labels");
    }
}

// short debug output to std::ostream
template <typename VALUE_TYPE>
void GenericMinkowskiFunctional<VALUE_TYPE>::dump (std::ostream &os) const {
    os << "# " << this->name ()
       << "\n# label v\n";
    for (int i = 0; i != (int)my_acc.size (); ++i) {
        os << i << " " << my_acc[i] << "\n";
    }
}

inline SurfaceIntegral::SurfaceIntegral () {
    my_ref_vertex = vec_t (0., 0.);
}

inline void SurfaceIntegral::ref_vertex (const vec_t &m) {
    my_ref_vertex = m;
}

inline const vec_t &SurfaceIntegral::ref_vertex () const {
    return my_ref_vertex;
}


#if 0
class BasicMinkowskiValuator {
public:
    typedef Boundary::edge_iterator edge_iterator;
    typedef std::string string;

    BasicMinkowskiValuator ();
    void set_reference (vec_t);
    ~BasicMinkowskiValuator ();

    double scalar (const string &);
    vec_t  vector (const string &);
    mat_t  matrix (const string &);

    /*
    double W100 (const Boundary &);
    double W200 (const Boundary &);
    vec_t W110 (const Boundary &);
    vec_t W101 (const Boundary &);
    */

protected:
    vec_t edge_normal (edge_iterator);
    vec_t edge_vertex0 (edge_iterator);
    vec_t edge_vertex1 (edge_iterator);
    vec_t my_ref;
    Boundary *my_b;
    std::map <string, vec_t> vacc;
    std::map <string, double> sacc;
    std::map <string, mat_t> macc;
};

class SurfaceMinkowskiValuator : public BasicMinkowskiValuator {
public:
    void add_contourseg (Boundary *b, edge_iterator begin, edge_iterator end);
};
#endif // 0

#endif /* MINKOWSKIVALUATIONS_H_INCLUDED */
