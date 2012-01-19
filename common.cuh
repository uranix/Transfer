#ifndef __COMMON_CUH__
#define __COMMON_CUH__

/* Coalised read size (bytes) */
#define COALESCED_SIZE (64)
#define COALESCED_NUM(datatype) (COALESCED_SIZE/(sizeof(datatype)))

typedef uint32_t copy_unit;

#endif
