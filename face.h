#ifndef __FACE_H__
#define __FACE_H__

#include "util.h"

#pragma pack(push, 1)
typedef
struct ALIGN(8) _face {
	/* 0 */
	idx p[4]; // padd 3->4
	/* 4x4 = 16 */
	REAL s[3]; // S_i \times mathbf{n_i}
	/* 4x4 + 4x8 = 48 */
	REAL _padd[sizeof(REAL) == 8 ? 3 : 9];
} face;
#pragma pack(pop)

CT_ASSERT(sizeof(face) == 64);

#endif
