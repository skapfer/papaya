#ifndef MINKOWSKIVALUATIONS_H_INCLUDED 
#define MINKOWSKIVALUATIONS_H_INCLUDED 

#include "util.h"
#include <map>
#include <ostream>
#include <iomanip>


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
    static void dump_accu (std::ostream &os, double);
    static void dump_accu (std::ostream &os, const vec_t &);
    static void dump_accu (std::ostream &os, const mat_t &);
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

    // reference point for calculation of Minkowski tensors
    // with r != 0
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

// return accumulator for label
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
       << "\n#--label-----------------v\n";
    for (int i = 0; i != (int)my_acc.size (); ++i) {
        os << std::setw (8) << i << " ";
        this->dump_accu (os, my_acc[i]);
        os << "\n";
    }
}

template <typename VALUE_TYPE>
inline void GenericMinkowskiFunctional<VALUE_TYPE>::dump_accu (std::ostream &os, double v) {
    os << std::setw (17) << v;
}

template <typename VALUE_TYPE>
inline void GenericMinkowskiFunctional<VALUE_TYPE>::dump_accu (std::ostream &os, const vec_t &v) {
    os << "(" << std::setw (16) << v[0] << " " << std::setw (16) << v[1] << ")";
}

template <typename VALUE_TYPE>
inline void GenericMinkowskiFunctional<VALUE_TYPE>::dump_accu (std::ostream &os, const mat_t &v) {
    os << "((" << std::setw (15) << v(0,0) << " " << std::setw (16) << v(0,1) << ") ";
    os << "("  << std::setw (16) << v(1,0) << " " << std::setw (15) << v(1,1) << "))";
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

#endif /* MINKOWSKIVALUATIONS_H_INCLUDED */
