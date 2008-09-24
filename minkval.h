#ifndef MINKOWSKIVALUATIONS_H_INCLUDED 
#define MINKOWSKIVALUATIONS_H_INCLUDED 

#include "util.h"
#include <map>
#include <ostream>
#include <iomanip>


class AbstractMinkowskiFunctional {
public:
    typedef Boundary::edge_iterator edge_iterator;
    typedef int label_t;
    // can be increased; this value is this low to catch runaway processes.
    enum { MAX_LABELS = 20000 };

    AbstractMinkowskiFunctional ();
    virtual ~AbstractMinkowskiFunctional () { }

    virtual void add_contour (const Boundary &, edge_iterator begin, edge_iterator end) = 0;

    // reference point for calculation of Minkowski tensors
    // with r != 0
    void ref_vertex (label_t, const vec_t &);
    void global_ref_vertex (const vec_t &);

protected:
    const vec_t &ref_vertex (label_t) const;
private:
    vec_t global_ref_vertex_;
    void reszref (label_t) const;
    mutable std::vector <vec_t> my_ref_vertex_;
};

template <typename VALUE_TYPE>
class GenericMinkowskiFunctional : public AbstractMinkowskiFunctional {
public:
    typedef VALUE_TYPE value_t;
    typedef GenericMinkowskiFunctional <value_t> this_t;

public:
    GenericMinkowskiFunctional (const std::string &name);

    const std::string &name () const;

    void dump (std::ostream &) const;

    const value_t &value (label_t) const;

protected:
    value_t &acc (label_t);
    void reszacc (label_t);
    static void dump_accu (std::ostream &os, double);
    static void dump_accu (std::ostream &os, const vec_t &);
    static void dump_accu (std::ostream &os, const mat_t &);
private:
    static value_t dummy_acc;
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


void calculate_all_surface_integrals (const Boundary &b);

ScalarMinkowskiFunctional *create_w000 ();
ScalarMinkowskiFunctional *create_w100 ();
ScalarMinkowskiFunctional *create_w200 ();
VectorMinkowskiFunctional *create_w010 ();
VectorMinkowskiFunctional *create_w110 ();
VectorMinkowskiFunctional *create_w210 ();
MatrixMinkowskiFunctional *create_w020 ();
MatrixMinkowskiFunctional *create_w120 ();
MatrixMinkowskiFunctional *create_w102 ();
MatrixMinkowskiFunctional *create_w220 ();
MatrixMinkowskiFunctional *create_w211 ();


//
// inline implementation
//

inline AbstractMinkowskiFunctional::AbstractMinkowskiFunctional () {
    global_ref_vertex_ = vec_t (0., 0.);
}

// allocate new labels (i.e. reference vertices for more labels)
inline void AbstractMinkowskiFunctional::reszref (label_t l) const {
    if (l < MAX_LABELS) {
        unsigned new_size = std::max (my_ref_vertex_.size (), l+1u);
        my_ref_vertex_.resize (new_size, vec_t (0., 0.));
    } else {
        die ("AbstractMinkowskiFunctional: too many labels");
    }
}

inline void AbstractMinkowskiFunctional::ref_vertex (
        AbstractMinkowskiFunctional::label_t l, const vec_t &m) {
    assert (l != Boundary::NO_LABEL);
    reszref (l);
    my_ref_vertex_.at (l) = m;
}

inline void AbstractMinkowskiFunctional::global_ref_vertex (
        const vec_t &ref) {
    my_ref_vertex_.clear ();
    global_ref_vertex_ = ref;
}

inline const vec_t &AbstractMinkowskiFunctional::ref_vertex (
        AbstractMinkowskiFunctional::label_t l) const {
    if (my_ref_vertex_.size ()) {
        if (l == Boundary::NO_LABEL) {
            return global_ref_vertex_;
        } else {
            reszref (l);
            return my_ref_vertex_.at (l);
        }
    } else {
        return global_ref_vertex_;
    }
}

template <typename VALUE_TYPE>
inline GenericMinkowskiFunctional<VALUE_TYPE>::GenericMinkowskiFunctional (const std::string &name)
    : my_name (name) {
        assert_not_nan (this_t::dummy_acc);
}

// return accumulator for label
template <typename VALUE_TYPE>
inline VALUE_TYPE &GenericMinkowskiFunctional<VALUE_TYPE>::acc (int label) {
    if (label == Boundary::NO_LABEL)
        return dummy_acc;
    reszacc (label);
    return my_acc[label];
}

template <typename VALUE_TYPE>
inline const VALUE_TYPE &GenericMinkowskiFunctional<VALUE_TYPE>::value (int label) const {
    assert (label >= 0);
    return const_cast <GenericMinkowskiFunctional<VALUE_TYPE> *> (this)->acc (label);
}

template <typename VALUE_TYPE>
inline const std::string &GenericMinkowskiFunctional <VALUE_TYPE>::name () const {
    return my_name;
}

// allocate new labels (i.e. accumulators for more labels)
template <typename VALUE_TYPE>
void GenericMinkowskiFunctional<VALUE_TYPE>::reszacc (label_t l) {
    if (l < MAX_LABELS) {
        unsigned new_size = std::max (my_acc.size (), l+1u);
        value_t zero;
        memset (&zero, 0, sizeof (zero));
        my_acc.resize (new_size, zero);
    } else {
        die ("GenericMinkowskiFunctional: too many labels");
    }
}

// short debug output to std::ostream
template <typename VALUE_TYPE>
void GenericMinkowskiFunctional<VALUE_TYPE>::dump (std::ostream &os) const {
    os << "# " << this->name ()
       << "\n#--label-----------------v-----------------v"
          "-----------------v-----------------v\n";
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

#endif /* MINKOWSKIVALUATIONS_H_INCLUDED */
