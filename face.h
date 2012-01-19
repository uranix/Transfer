#ifndef __FACE_H__
#define __FACE_H__

#include "util.h"

#pragma pack(push, 1)
struct __align__(8) face {
	/* 0 */
	index p[4]; // padd 3->4
	/* 4x4 = 16 */
	REAL s[4]; // S_i \times mathbf{n_i}
	/* 4x4 + 4x8 = 48 */
	REAL _padd[2];
};
#pragma pack(pop)

CT_ASSERT(sizeof(face) == 64);

#endif
