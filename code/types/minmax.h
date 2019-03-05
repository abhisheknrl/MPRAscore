#ifndef __MINMAX_H__
#define __MINMAX_H__

// The __min and __max macros should be defined in stdlib.h
#include <stdlib.h>

// However, sometimes they are not (e.g., with gcc), in which case we define them post hoc
#ifndef __min
#define __min(a,b) ((a)<(b) ? (a) : (b))
#endif
#ifndef __max
#define __max(a,b) ((a)>(b) ? (a) : (b))
#endif
// #ifndef __clamp
#define __clamp(x,xmin,xmax) (__max(xmin, __min(x, xmax)))
// #endif
#ifndef __abs
#define __abs(a) ((a)>=0 ? (a) : -(a))
#endif


#endif
