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

#ifndef PAPAYA2
typedef Eigen::Vector2d vec_t;
typedef Eigen::Matrix2d mat_t;
#else
#ifdef DIMENSION
#if (DIMENSION==2)
typedef vec2_t vec_t;
typedef mat2_t mat_t;
#elif (DIMENSION==3)
typedef mat3_t mat_t;
typedef vec3_t vec_t;
#endif // DIMENSION == 2, 3
#endif // DIMENSION
#endif

// TensorRank <N>::value_type is the class which stores
// rank-N tensors.
template <int RANK>
struct TensorRank;



void assert_not_nan (double x);
void assert_not_nan (const vec2_t &);
void assert_non_nan (const vec3_t &);
void assert_not_nan (const mat2_t &m);


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

inline void assert_not_nan (double x) {
#ifndef NDEBUG
    assert (!isnan (x));
#endif // NDEBUG
}

inline void assert_not_nan (const vec2_t &v) {
#ifndef NDEBUG
    assert_not_nan (v(0));
    assert_not_nan (v(1));
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