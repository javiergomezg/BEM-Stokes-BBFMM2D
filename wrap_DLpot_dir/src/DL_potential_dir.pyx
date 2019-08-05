import cython
import numpy as np
cimport numpy as np

cdef extern from "DL_potential_dir.hpp":

	ctypedef double REAL

	void direct_DL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *ux, int n_ux,
			REAL *uy, int n_uy,
			REAL *nx, int n_nx,
			REAL *ny, int n_ny,
			REAL *b, int n_b);


def direct_DL_potential(np.ndarray[REAL, ndim = 1] xsrc,
			np.ndarray[REAL, ndim = 1] ysrc,
			np.ndarray[REAL, ndim = 1] xtar,
			np.ndarray[REAL, ndim = 1] ytar,
			np.ndarray[REAL, ndim = 1] Rp,
			np.ndarray[REAL, ndim = 1] wc,
			np.ndarray[REAL, ndim = 1] ux,
			np.ndarray[REAL, ndim = 1] uy,
			np.ndarray[REAL, ndim = 1] nx,
			np.ndarray[REAL, ndim = 1] ny,
			np.ndarray[REAL, ndim = 1, mode = "c"] b):

	cdef np.int32_t xsrcSize = len(xsrc)
	cdef np.int32_t ysrcSize = len(ysrc)
	cdef np.int32_t xtarSize = len(xtar)
	cdef np.int32_t ytarSize = len(ytar)
	cdef np.int32_t n_pan = len(Rp)
	cdef np.int32_t n_wc = len(wc)
	cdef np.int32_t n_ux = len(ux)
	cdef np.int32_t n_uy = len(uy)
	cdef np.int32_t n_nx = len(nx)
	cdef np.int32_t n_ny = len(ny)
	cdef np.int32_t n_b = len(b)

	direct_DL_potential_cpp(<REAL*> &xsrc[0], <int> xsrcSize,
			<REAL*> &ysrc[0], <int> ysrcSize,
			<REAL*> &xtar[0], <int> xtarSize,
			<REAL*> &ytar[0], <int> ytarSize,
			<REAL*> &Rp[0], <int> n_pan,
			<REAL*> &wc[0], <int> n_wc,
			<REAL*> &ux[0], <int> n_ux,
			<REAL*> &uy[0], <int> n_uy,
			<REAL*> &nx[0], <int> n_nx,
			<REAL*> &ny[0], <int> n_ny,
			<REAL*> &b[0], <int> n_b)
	
