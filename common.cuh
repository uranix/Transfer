#ifndef __COMMON_CUH__
#define __COMMON_CUH__

#ifdef _MSC_VER
//#pragma warning (push, 4710)
#endif

/* Coalesced read size (bytes) */
#define COALESCED_SIZE (64)
#define COALESCED_NUM(datatype) (COALESCED_SIZE/(sizeof(datatype)))

typedef uint32_t copy_unit;

#endif
