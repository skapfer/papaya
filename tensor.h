#ifndef PAPAYA_TENSOR_H_INCLUDED 
#define PAPAYA_TENSOR_H_INCLUDED 

// classes from Eigen
#include "vector.h"
#include "matrix.h"

// FIXME should be replaced by a symmetric matrix type.
typedef Eigen::Matrix2d mat2_t;
typedef Eigen::Matrix3d mat3_t;
typedef Eigen::Vector2d vec2_t;
typedef Eigen::Vector3d vec3_t;

typedef Eigen::Vector2d vec_t;
typedef Eigen::Matrix2d mat_t;

// TensorRank <N>::value_type is the class which stores
// rank-N tensors.
template <int RANK>
struct TensorRank;



void assert_not_nan (double x);
void assert_not_nan (const vec2_t &);
void assert_not_nan (const vec3_t &);
void assert_not_nan (const mat2_t &m);

namespace Papaya2 {
    bool is_finite (double);

    inline bool is_finite (const vec2_t &v) {
        return is_finite (v[0]) && is_finite (v[1]);
    }

    inline bool is_finite (const vec3_t &v) {
        return is_finite (v[0]) && is_finite (v[1])
               && is_finite (v[2]);
    }
}


//
// INLINES AND TEMPLATES
//

template <>
struct TensorRank<2> {
#ifdef DIMENSION
    typedef mat_t value_type;
#endif
};

template <>
struct TensorRank<1> {
#ifdef DIMENSION
    typedef vec_t value_type;
#endif
};

template <>
struct TensorRank<0> {
#ifdef DIMENSION
    typedef double value_type;
#endif
};

// FIXME move tensor products etc. here.

// seems we can't rely on isnan being available everywhere.
inline bool not_nan (double x) {
    return x==x;
}

inline bool not_nan (const vec2_t &v) {
    return not_nan (v[0]) && not_nan (v[1]);
}

inline void assert_not_nan (double x) {
#ifndef NDEBUG
    assert (not_nan (x));
#endif // NDEBUG
}

inline void assert_not_nan (const vec2_t &v) {
#ifndef NDEBUG
    assert (not_nan (v));
#endif // NDEBUG
}

inline void assert_not_nan (const vec3_t &v) {
#ifndef NDEBUG
    assert_not_nan (v(0));
    assert_not_nan (v(1));
    assert_not_nan (v(2));
#endif // NDEBUG
}

inline void assert_not_nan (const mat2_t &m) {
#ifndef NDEBUG
    assert_not_nan (m(0,0));
    assert_not_nan (m(1,0));
    assert_not_nan (m(0,1));
    assert_not_nan (m(1,1));
#endif // NDEBUG
}

inline void assert_not_nan (const mat3_t &m) {
#ifndef NDEBUG
    for (int i = 0; i != 3; ++i)
    for (int j = 0; j != 3; ++j) {
        assert_not_nan (m(i,j));
    }
#endif // NDEBUG
}

#endif /* PAPAYA_TENSOR_H_INCLUDED */
