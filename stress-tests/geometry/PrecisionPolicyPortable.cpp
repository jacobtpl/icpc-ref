// PrecisionPolicy with the plain C++ __float128 math that Point.h uses
// where glibc's *f128 functions are missing (e.g. macOS).
#define GEO_OWN_F128
#include "PrecisionPolicy.cpp"
