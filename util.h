#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdint.h>

typedef uint32_t index;
#define ASSERT_CONCAT_(a, b) a##b
#define ASSERT_CONCAT(a, b) ASSERT_CONCAT_(a, b)
#define CT_ASSERT(e) enum { ASSERT_CONCAT(assert_line_, __LINE__) = 1/(!!(e)) }

/* TODO test */
#define align_power(x, y) (((x)+(y)-1) & (~((y)-1)))
#define upper_div(x, y) (((x)+(y)-1)/(y))
#define align_num(x, y) (upper_div(x,y)*y)

typedef double REAL;

#endif
