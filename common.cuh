#ifndef __COMMON_CUH__
#define __COMMON_CUH__

/* Coalesced read size (bytes) */
#define COALESCED_SIZE (64)
#define COALESCED_NUM(datatype) (COALESCED_SIZE/(sizeof(datatype)))
#define ASLM_MAX (128)

#endif
