#ifndef __COMMON_CUH__
#define __COMMON_CUH__

/* Coalised read size (bytes) */
#define COALISED_SIZE (64)
#define COALISED_NUM(datatype) (COALISED_SIZE/(sizeof(datatype)))

typedef uint32_t copy_unit;

#endif
