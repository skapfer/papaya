// this requires a C99 feature and may have to be replaced
//
// for some reason, g++-4.2.4-1ubuntu3 does not expose isfinite from <math.h>
// but only from <cmath>, even though it is C99.  go figure.

#include <cmath>
#include "tensor.h"
using namespace std;

namespace Papaya2 {
    bool is_finite (double x) {
        return !!isfinite (x);
    }
}
