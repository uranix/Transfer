#ifndef __TETRAHEDRON_H__
#define __TETRAHEDRON_H__

#include "util.h"

#pragma pack(push, 1)
typedef struct ALIGN(8) _tetrahedron {
	/* 0 */
	idx p[4];
	/* 4x4 = 16 */
	REAL kappa_volume;
	REAL I_p;
	/* 4x4 + 2x8 = 32 */
	REAL s[4][3]; // S_i \times mathbf{n_i}
	/* 4x4 + 2x8 + 12x8 = 128 b */
} tetrahedron;
#pragma pack(pop)

CT_ASSERT(sizeof(tetrahedron) == 128);


#endif
