#ifndef MINK3_H_INCLUDED 
#define MINK3_H_INCLUDED 
#ifndef PAPAYA2
#error "This header is PAPAYA2 code and cannot be mixed with PAPAYA1."
#endif

#include "util.h"
#include <stdio.h> // snprintf


//
// C++ SCAFFOLDING
//

template <int RANK>
class TensorValuation {
public:
    typedef typename TensorRank <RANK>::value_type value_type;

    const value_type &value (int label) const {
        assert (label >= 0);
        return *value_reg (label);
    }

    const value_type global_value () const {
        value_type ret;
        // Eigen classes don't initialize themselves
        memset (&ret, 0, sizeof (value_type));
        for (int i = 1; i != (int)stor_.size (); ++i) {
            ret += stor_[i];
        }
        return ret;
    }

protected:
    void add_to_value_reg (int label, double prefactor, const value_type &contribution) {
        value_type *labelled_ = value_reg (label);
        *labelled_ += prefactor * contribution;
    }

    value_type *value_reg (int label) {
        ++label;
        assert (label >= 0);    // -1 is the "rubbish dump" label
        if (stor_.size () <= label)
            resize_stor_ (label+1);
        return &stor_[label];
    }

    void resize_stor_ (int size) {
        int old_size = stor_.size ()
        assert (size > old_size);
        stor_.resize (size);
        // Eigen classes don't initialize themselves
        memset (&stor_[old_size], 0, sizeof (value_type) * (size-old_size));
    }

    std::vector <value_type> stor_;
};

// class MinkowskiValuation contains the common scaffolding for the individual
// valuations.  it is not useful on its own.
template <int NU, int R, int S>
class MinkowskiValuation : public TensorValuation <R+S> {
public:
    typedef typename TensorValuation <R+S>::value_type value_type;

    const std::string name () const {
        return this->fct_name_;
    }

protected:
    MinkowskiValuation () {
        snprintf (this->fct_name_, sizeof (this->fct_name_),
                  "W%i%i%i", NU, R, S);
#if (DIMENSION == 2)
        if (NU == 0)
            normal_ = 1.;
        else
            normal_ = .5;
#elif (DIMENSION == 3)
        if (NU == 0)
            normal_ = 1.;
        else
            normal_ = 1./3;
#else
        #error "DIMENSION is not set sensibly."
#endif
        assert (NU >= 0);
        assert (NU < DIMENSION+1);
        assert (R >= 0);
        assert (S >= 0);
    }

    template <typename BODY>
    void add (const BODY &body, int facet, const value_type &contribution) {
        int label = body.facet_label (facet);
        this->add_to_value_reg (label, normal_, contribution);
    }

private:
    double normal_;
    char fct_name_[20];
};


//
// FUNCTIONALS PROPER
//

struct W_100 : public MinkowskiValuation <1, 0, 0> {
    template <typename BODY>
    void operator() (const BODY &body, int facet) {
        add (body, facet, body.facet_volume (facet));
    }
};

struct W_200 : public MinkowskiValuation <2, 0, 0> {
    template <typename BODY>
    void operator() (const BODY &body, int facet) {
        double contrib = 0;
        add (body, facet, contrib);
    }
};


#endif // MINK3_H_INCLUDED
