import cython
import numpy as np
cimport numpy as np

cdef extern from "SL_potential_dir.hpp":

	ctypedef double REAL

	void direct_SL_potential_cpp(REAL *xsrc, int xsrcSize,
			REAL *ysrc, int ysrcSize,
			REAL *xtar, int xtarSize,
			REAL *ytar, int ytarSize,
			REAL *Rp, int n_pan,
			REAL *wc, int n_wc,
			REAL *fx, int n_fx,
			REAL *fy, int n_fy,
			REAL Re,
			REAL *Ax, int n_Ax);


def direct_SL_potential(np.ndarray[REAL, ndim = 1] xsrc,
			np.ndarray[REAL, ndim = 1] ysrc,
			np.ndarray[REAL, ndim = 1] xtar,
			np.ndarray[REAL, ndim = 1] ytar,
			np.ndarray[REAL, ndim = 1] Rp,
			np.ndarray[REAL, ndim = 1] wc,
			np.ndarray[REAL, ndim = 1] fx,
			np.ndarray[REAL, ndim = 1] fy,
			Re,
			np.ndarray[REAL, ndim = 1, mode = "c"] Ax):

	cdef np.int32_t xsrcSize = len(xsrc)
	cdef np.int32_t ysrcSize = len(ysrc)
	cdef np.int32_t xtarSize = len(xtar)
	cdef np.int32_t ytarSize = len(ytar)
	cdef np.int32_t n_pan = len(Rp)
	cdef np.int32_t n_wc = len(wc)
	cdef np.int32_t n_fx = len(fx)
	cdef np.int32_t n_fy = len(fy)
	cdef np.int32_t n_Ax = len(Ax)

	direct_SL_potential_cpp(<REAL*> &xsrc[0], <int> xsrcSize,
			<REAL*> &ysrc[0], <int> ysrcSize,
			<REAL*> &xtar[0], <int> xtarSize,
			<REAL*> &ytar[0], <int> ytarSize,
			<REAL*> &Rp[0], <int> n_pan,
			<REAL*> &wc[0], <int> n_wc,
			<REAL*> &fx[0], <int> n_fx,
			<REAL*> &fy[0], <int> n_fy,
			<REAL> Re,
			<REAL*> &Ax[0], <int> n_Ax)
	
