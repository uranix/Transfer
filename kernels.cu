#include "util.h"
#include "tetrahedron.h"
#include "face.h"

#include "common.cuh"

/*
Computes r = (e_i, f) + (1/kappa nabla e_i, 1/kappa nabla f) - (e_i, I_p)
	nP						: total vertices number
	start[nP+1]				: start[i+1] - start[i] = number of tetrahedrons incidental to vertex i
	idx[start[nP]]			: corresponding tetrahedron index
	pos[start[nP]]			: local vertex index in tetradedron
	mesh[nT]				: mesh
	slm						: anglar harmonics total number. aslm = align_power(slm, COALISED_NUM(REAL))
	omega[3*aslm*aslm]		: values part of <Omega_i Omega_j>_lm^lms. Symmetrical to i <-> j, lm <-> lms.
							omega[i][j][lm][lms] = omega[3*aslm*aslm + aslm*lms + lm] when j = omega_pos[i][lms][lm], otherwise 0
							- coalesced when read thread <-> lm
	omega_pos[3*aslm*aslm]	: coordinate part of <Omega_i Omega_j>_lm^lms
							omega_pos[i][lms][lm] = omega_pos[3*aslm*aslm + aslm*lms + lm]
							- coalesced when read thread <-> lm
	f[aslm*nP]				: degrees of freedom 
							f[point][lm]=f[aslm*point + lm]. 
							- coalesced when reading thread <-> lm
	r[aslm*nP]				: result

Assumed:
	blockDim.x = aslm, should be less or equal ASLM_MAX.
	blockDim.y = 1
	blockDim.z = 1
	gridDim.x*gridDim.y = nP
	gridSize.z = 1

	shmem per block = 128b + 32b * blockDim.x + ? [__syncthreads()]
		~4 Kb  blockDim.x = 128
		~8 Kb blockDim.x = 256
		~32 Kb blockDim.x = 1024
   */

#define ASLM_MAX (256)

__global__ void volumePart(	index nP, index *start, index *idx, index *pos, tetrahedron *mesh, 
							index slm, REAL *omega, index *omega_pos, REAL *f, REAL *r) {
	__shared__ tetrahedron tetas;
	__shared__ REAL sums_j[4*ASLM_MAX]; /* 3 -> 4 for align*/
	index aslm = align_power(slm, COALISED_NUM(REAL));
	index vertex = blockIdx.x + blockIdx.y * gridDim.x;
	index lm = threadIdx.x;
	REAL sum_i[3];
	REAL *sum_j = &sums_j[4*lm];
	index lo, hi;
	REAL fl[4];
	REAL fc;

	REAL sum = 0;
	lo = start[vertex];
	hi = start[vertex+1];
	for (int j = lo; j < hi; j++) {
		__syncthreads();
		tetrahedron *tet = &tetas;
		copy_unit *dst = (copy_unit *)tet;
		copy_unit *src = (copy_unit *)(mesh + idx[j]);
		index local = pos[j];
		
		index copy_incr = blockDim.x;
		for (int s=0, smax = sizeof(tetrahedron); s < smax; s += sizeof(copy_unit)*copy_incr, dst += copy_incr, src += copy_incr)
			if (s + sizeof(copy_unit) * threadIdx.x < smax)
				dst[threadIdx.x] = src[threadIdx.x];
		__syncthreads();
		fc = 0;
		#pragma unroll
		for (int s = 0; s < 4; s++) {
			REAL tmp = f[aslm*tet->p[s] + lm];
			fl[s] = tmp; 
			fc += tmp;
		}
		sum += tet->kappa_volume * (fl[local] + fc) * (1. / 20.);
		if (lm == 0)
			sum -= tet->kappa_volume * tet->I_p * (1. / 4.);
		#pragma unroll
		for (int si = 0; si < 3; si++)
			sum_i[si] = tet->s[local][si] / tet->kappa_volume * (1. / 9.); 
		#pragma unroll
		for (int sj = 0; sj < 3; sj++) {
			sum_j[sj] = 0;
			#pragma unroll
			for (int k = 0; k<4; k++)
				sum_j[sj] += fl[k] * tet->s[k][sj];
		}
		__syncthreads();
		REAL rowsum;
		index v = lm;
		#pragma unroll
		for (int row = 0; row < 3; row++, v += (aslm-slm)*aslm) {
			rowsum = 0;
			for (int lms = 0; lms < slm; lms++, v += aslm) {
				rowsum += omega[v] * sums_j[4*lms + omega_pos[v]];
			}
			sum += rowsum * sum_i[row];
		}
	}
	r[aslm*vertex + lm] = sum;
}

/*
WORKS ONLY IF NORMAL IS (+/-1,0,0), (0,+/-1,0) or (0,0,+/-1).
Computes r += int_{dG x 4pi} |Omega n(x)| e_i f d Omega dS
	nP						: total vertices number
	start[nP+1]				: start[i+1] - start[i] = number of faces incidental to vertex i
	idx[start[nP]]			: corresponding face index
	pos[start[nP]]			: local vertex index in face
	bnd[nT]					: boundary faces
	slm						: anglar harmonics total number. aslm = align_power(slm, COALISED_NUM(REAL))
	Ox,Oy,Oz[aslm*aslm]		: <|Omega_x|>, <|Omega_y|>, <|Omega_z|>.
	f[aslm*nP]				: degrees of freedom 
							f[point][lm]=f[aslm*point + lm]. 
							- coalesced when reading thread <-> lm
	r[aslm*nP]				: result

Assumed:
	blockDim.x = aslm, should be less or equal ASLM_MAX.
	blockDim.y = 1
	blockDim.z = 1
	gridDim.x*gridDim.y = nP
	gridSize.z = 1

	shmem per block = 32*ASLM_MAX * blockDim.x + ? [__syncthreads()]
   */
__global__ void surfacePart( index nP, index *start, index *idx, index *pos, face *bnd, 
							 index slm, REAL *Ox, REAL *Oy, REAL *Oz, REAL *f, REAL *r) 
{
	__shared__ REAL fv[4*ASLM_MAX]; /* f1,f2,f3,fc */
	REAL *On;
	face triangle;

	index aslm = align_power(slm, COALISED_NUM(REAL));	
	index vertex = blockIdx.x + blockIdx.y * gridDim.x;
	index lm = threadIdx.x;
	index lo, hi;

	REAL sum = 0;
	lo = start[vertex];
	hi = start[vertex+1];
	for (int j = lo; j < hi; j++) {
		__syncthreads();
		face *tr = &triangle;
		copy_unit *dst = (copy_unit *)tr;
		copy_unit *src = (copy_unit *)(bnd + idx[j]);
		index local = pos[j];
		
		index copy_incr = blockDim.x;
		for (int s=0, smax = sizeof(face); s < smax; s += sizeof(copy_unit)*copy_incr, dst += copy_incr, src += copy_incr)
			if (s + sizeof(copy_unit) * threadIdx.x < smax)
				dst[threadIdx.x] = src[threadIdx.x];
		__syncthreads();
		REAL surf = abs(tr->s[0] + tr->s[1] + tr->s[2]);
		if (2*abs(tr->s[0]) > surf)
			On = Ox;
		else if (2*abs(tr->s[1]) > surf)
			On = Oy;
		else
			On = Oz;
		fv[3*aslm+lm] = 0;
		#pragma unroll
		for (int s = 0; s<3; s++) {
			REAL tmp = f[aslm*tr->p[s]+lm];
			fv[s*aslm+lm] = tmp;
			fv[3*aslm+lm] += tmp;
		}
		__syncthreads();
		for (int lms = 0; lms < slm; lms ++)
			sum += surf * On[aslm * lms + lm] * (fv[local*aslm + slm] + fv[3*aslm + slm]) * (1. / 12.);
	}
	r[aslm*vertex + lm] += sum;
}
