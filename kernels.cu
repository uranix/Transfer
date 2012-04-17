#include "util.h"

#include "MeshData.h"
#include "AngularData.h"

#include "common.cuh"

/*
   sz must be multiple of sizeof(copy_unit);
   dst and src should be aligned of sizeof(copy_unit) boundary
 */
__device__ void dmemcpy(void *_dst, const void *_src, idx sz) {
	copy_unit *dst = (copy_unit *)_dst;
	copy_unit *src = (copy_unit *)_src;
	for (idx s = threadIdx.x, smax = sz / sizeof(copy_unit); s < smax; s++)
		dst[s] = src[s];
}

/*
Computes right-hand size of SLAE r = (e_i, I_p)
 */

__global__ void rightHandSide(	DeviceMeshDataRaw md,
								REAL *r)
{
	__shared__ tetrahedron stet;
	tetrahedron *tet = &stet;

	idx *start = md.tetstart;
	idx *tetidx = md.tetidx;
	idx *pos = md.tetpos;
	tetrahedron *mesh = md.mesh; 

	idx vertex = blockIdx.x + blockIdx.y * gridDim.x;
	idx lm = threadIdx.x;
	idx lo, hi;

	REAL sum = 0;
	lo = start[vertex];
	hi = start[vertex+1];
	for (int j = lo; j < hi; j++) {
		__syncthreads();
		idx local = pos[j];
		
		dmemcpy(tet, mesh+tetidx[j], sizeof(tetrahedron));
		__syncthreads();
		if (lm == 0)
			sum += tet->kappa_volume * tet->I_p * (1. / 4.);
	}
	if (vertex < md.nP)
		r[aslm*vertex + lm] = sum;
}

/*
Computes r = (e_i, f) + (1/kappa nabla e_i, 1/kappa nabla f)

	nP						: total vertices number
	start[nP+1]				: start[i+1] - start[i] = number of tetrahedrons incidental to vertex i
	idx[start[nP]]			: corresponding tetrahedron idx
	pos[start[nP]]			: local vertex idx in tetradedron
	mesh[nT]				: mesh
	slm						: angular harmonics total number. aslm = align_power(slm, COALESCED_NUM(REAL))
	omega[3*aslm*aslm]		: values part of <Omega_i Omega_j>_lm^lms. Symmetrical to i <-> j, lm <-> lms.
							omega[i][j][lm][lms] = omega[i*aslm*aslm + aslm*lms + lm] when j = omega_pos[i][lms][lm], otherwise 0
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

	shmem per block = 128b + 24b * blockDim.x + ? [__syncthreads()]
		~4 Kb  blockDim.x = 128
		~8 Kb blockDim.x = 256
		~32 Kb blockDim.x = 1024
   */

#define ASLM_MAX (256)

__global__ void volumePart(	DeviceMeshDataRaw md,
							DeviceAngularDataRaw ad,
							REAL *f,
							REAL *r) 
{
	__shared__ tetrahedron stet;
	__shared__ REAL sums_j[3*ASLM_MAX]; 
	tetrahedron *tet = &stet;

	idx slm = ad.slm;
	idx aslm = ad.aslm;
	REAL *omega = ad.omega;
	idx *omega_pos = ad.omega_pos;

	idx *start = md.tetstart;
	idx *tetidx = md.tetidx;
	idx *pos = md.tetpos;
	tetrahedron *mesh = md.mesh; 

	idx vertex = blockIdx.x + blockIdx.y * gridDim.x;
	idx lm = threadIdx.x;
	REAL sum_i[3];
	idx lo, hi;
	REAL fl[4];
	REAL fsum;

	REAL sum = 0;
	lo = start[vertex];
	hi = start[vertex+1];
	for (int j = lo; j < hi; j++) {
		__syncthreads();
		idx local = pos[j];
		
		dmemcpy(tet, mesh+tetidx[j], sizeof(tetrahedron));
		__syncthreads();
		fc = 0;
		#pragma unroll
		for (int s = 0; s < 4; s++) {
			REAL tmp = f[aslm*tet->p[s] + lm];
			fl[s] = tmp; 
			fsum += tmp;
		}
		sum += tet->kappa_volume * (fl[local] + fsum) * (1. / 20.);
		#pragma unroll
		for (int si = 0; si < 3; si++)
			sum_i[si] = tet->s[local][si] / tet->kappa_volume * (1. / 9.); 
		#pragma unroll
		for (idx sj = 0, v = lm; sj < 3; sj++, v += aslm) {
			sums_j[v] = 0;
			#pragma unroll
			for (int k = 0; k<4; k++)
				sums_j[v] += fl[k] * tet->s[k][sj];
		}
		__syncthreads();
		REAL rowsum;
		/* nvcc can't unroll this.
		idx v = lm;
		#pragma unroll
		for (int row = 0; row < 3; row++) {
			rowsum = 0;
			for (int lms = 0; lms < slm; lms++, v += aslm) {
				rowsum += omega[v] * sums_j[omega_pos[v]*aslm + lms];
			}
			v += (aslm-slm)*aslm;
			sum += rowsum * sum_i[row];
		}
		*/
		rowsum = 0;
		for (idx lms = 0, v = lm; lms < slm; lms++, v += aslm)
			rowsum += omega[v] * sums_j[omega_pos[v]*aslm + lms]; /* TODO: resolve bank conflict (x16) */
		sum += rowsum * sum_i[0];

		rowsum = 0;
		for (idx lms = 0, v = lm + aslm*aslm; lms < slm; lms++, v += aslm)
			rowsum += omega[v] * sums_j[omega_pos[v]*aslm + lms]; /* TODO: resolve bank conflict (x16) */
		sum += rowsum * sum_i[1];

		rowsum = 0;
		for (idx lms = 0, v = lm + aslm*aslm*2; lms < slm; lms++, v += aslm)
			rowsum += omega[v] * sums_j[omega_pos[v]*aslm + lms]; /* TODO: resolve bank conflict (x16) */
		sum += rowsum * sum_i[2];
	}
	if (vertex < md.nP)
		r[aslm*vertex + lm] = sum;
}

/*
WORKS ONLY IF NORMAL IS (+/-1,0,0), (0,+/-1,0) or (0,0,+/-1). Issue #15
Computes r += int_{dG x 4pi} |Omega n(x)| e_i f d Omega dS
	nP						: total vertices number
	start[nP+1]				: start[i+1] - start[i] = number of faces incidental to vertex i
	idx[start[nP]]			: corresponding face idx
	pos[start[nP]]			: local vertex idx in face
	bnd[nF]					: boundary faces
	slm						: angular harmonics total number. aslm = align_power(slm, COALESCED_NUM(REAL))
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
__global__ void surfacePart( DeviceMeshDataRaw md,
							 DeviceAngularDataRaw ad,
							 REAL *f, 
							 REAL *r) 
{
	__shared__ REAL fv[4*ASLM_MAX]; /* f1,f2,f3,fsum */
	__shared__ face stri;
	face *tr = &stri;
	REAL *On;

	idx *start = md.facestart;
	idx *faceidx = md.faceidx;
	idx *pos = md.facepos;
	face *bnd = md.bnd;

	idx slm = ad.slm; 
	idx aslm = ad.aslm;
	REAL *Ox = ad.Ox;
	REAL *Oy = ad.Oy;
	REAL *Oz = ad.Oz;

	idx vertex = blockIdx.x + blockIdx.y * gridDim.x;
	idx lm = threadIdx.x;
	idx lo, hi;

	REAL sum = 0;
	lo = start[vertex];
	hi = start[vertex+1];
	for (int j = lo; j < hi; j++) {
		__syncthreads();
		idx local = pos[j];	
		dmemcpy(tr, bnd + faceidx[j], sizeof(tr));
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
	if (vertex < md.nP)
		r[aslm*vertex + lm] += sum;
}
