import numpy as np
import Classnode
import cuadratura
#import fs_gf
import time
from time import time
from numpy import pi, cos, sin
from DL_potential import bbfmm_DL_potential
from DL_potential_dir import direct_DL_potential
from SL_potential import bbfmm_SL_potential
from SL_potential_dir import direct_SL_potential
from ufield_dir import direct_ufield
from pfield_dir import direct_pfield
from ufield_fmm import fmm_ufield
from pfield_fmm import fmm_pfield


def targetvector(Ncol,Xcol,Ycol):
    target = Target(Ncol)
    target.z[:].real = Xcol[:]    
    target.z[:].imag = Ycol[:]
    return target
    
def sourcevector(Xpan, Ypan, t):
    N_sources = (len(Xpan)-1)*len(t)
    sources = Sources(N_sources)
    #print len(Xpan)
    #print len(t)    
    
    #print sources.z
    #print len(sources.z)
    
    i = 0
    for p in range(len(Xpan)-1):
        Xg, Yg = nodos_gauss_panel(t, Xpan[p], Ypan[p], Xpan[p+1], Ypan[p+1])
        #print Xpan[p], Ypan[p]
        #print Xg
        #print Yg
        #print Xpan[p+1], Ypan[p+1]
        
        for j in range(len(t)):
            sources.z[i] = Xg[len(t)-1-j] + 1j*Yg[len(t)-1-j]
            sources.w[i] = (pi/len(t))*(1.0-(t[len(t)-1-j])**2)**(0.5)
            #sources.z[i] = Xg[j] + 1j*Yg[j]
            #sources.w[i] = (1.0-t[j]**2)**(0.5)
            i+=1
            
    #print sources.z
    #print
    #print t
    #print
    #print sources.w
    
    return sources
    
def vector_b_p2p(T,S,nc, ux, uy, Rp,nx, ny):
    print "Entrando a calculo directo de double layer potential"
    Nq = len(T.z)
    Ns = len(S.z)

    T_load = time()    
    
    #print Nq, Ns
    
    bx,by = np.zeros(Nq, dtype=np.float64), np.zeros(Nq, dtype=np.float64)
    bx[:],by[:] = ux[:], uy[:]

    #print bx, ux
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag  
       
    #print X
    #print Y    
    
    for q in range(Nq):
        X0 = T.z[q].real
        Y0 = T.z[q].imag
        Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy =gf_stress_fs(X,Y,X0,Y0)

        if (q % 1000 == 0):
            print "Se esta calculando el punto de colocacion:", q
            print "tiempo transcurrido: ", time() - T_load, "[s]" 
        
        for ic in range(Ns):
            if (ic//nc != q):
                bx[q] = bx[q] + (-1./(2.*pi))*(Rp[ic//nc]/2.)*S.w[ic]*(ux[ic//nc]*
                (Txxx[ic]*nx[ic//nc]+Txxy[ic]*ny[ic//nc]) + uy[ic//nc]*
                (Tyxx[ic]*nx[ic//nc]+Tyxy[ic]*ny[ic//nc]))
                
                by[q] = by[q] + (-1./(2.*pi))*(Rp[ic//nc]/2.)*S.w[ic]*(ux[ic//nc]*
                (Txyx[ic]*nx[ic//nc]+Txyy[ic]*ny[ic//nc]) + uy[ic//nc]*
                (Tyyx[ic]*nx[ic//nc]+Tyyy[ic]*ny[ic//nc]))
                
    b = np.zeros([2*Nq], dtype = np.float64)

    b[:Nq], b[Nq:] = bx[:],by[:]
    
    
    return b

def vector_b_p2p_cy(T,S,nc, ux, uy, Rp,nx, ny):
    print "Entrando a calculo directo de double layer potential"
    Nq = len(T.z)
    Ns = len(S.z)

    T_load = time()    
    
    #print Nq, Ns
    #print bx, ux
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)
    wc = np.zeros(Ns)
    XT = np.zeros(Nq)
    YT = np.zeros(Nq)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag
    
    wc[:] = S.w[:]
    
    XT[:] = T.z[:].real
    YT[:] = T.z[:].imag

    b = np.zeros([2*Nq], dtype = np.float64)
       
    direct_DL_potential(X,
			Y,
			XT,
			YT,
			Rp,
			wc,
			ux,
			uy,
			nx,
			ny,
			b)

    del X, Y, wc, XT, YT                
    
    return b
    

def vector_b_fmm(T,S,nc, ux, uy, Rp,nx, ny, nchebpol):
    Nq = len(T.z)
    Ns = len(S.z)
    
    #print Nq, Ns
       
    #print bx, ux
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)
    XT = np.zeros(Nq)
    YT = np.zeros(Nq)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag  
    XT[:] = T.z[:].real
    YT[:] = T.z[:].imag
    
    
    #print X
    #print Y    

    bfmm = np.zeros([2*Nq], dtype = np.float64)
    
    bbfmm_DL_potential(X,
			Y,
			XT,
			YT,
			Rp,
			S.w[:],
			ux,
			uy,
			nx,
			ny,
			bfmm,
                nchebpol)

    del X, Y, XT, YT   
    
    return bfmm
    
def vector_Ax_p2p(T,S,nc,f,Rp,Xp,Yp,mu):
    print "entrando a calculo de directo de sigle layer potential"
    Nq, Ns = len(T.z), len(S.z)
    T_load = time()

    fx, fy = np.zeros(Nq, dtype = np.float64), np.zeros(Nq, dtype = np.float64)
    fx[:], fy[:] = f[:Nq], f[Nq:]
    
    Ax = np.zeros(2*Nq, dtype = np.float64)

        
    for q in range(Nq):
        Ax[q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
                int_alpha_xx(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
                int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])
               
        Ax[Nq + q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
                int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
                int_alpha_yy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])
                
    #print Ax
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)

    X[:] = S.z[:].real
    Y[:] = S.z[:].imag    
    
    for q in range(Nq):
        X0 = T.z[q].real
        Y0 = T.z[q].imag        
        Gxx, Gxy, Gyx, Gyy = gf_vel_fs(X,Y,X0,Y0)

        if (q % 1000 == 0):
            print "Se esta calculando el punto de colocacion:", q
            print "tiempo transcurrido: ", time() - T_load, "[s]"         
        
        for ic in range(Ns):
            if (ic//nc != q):
                Ax[q] = Ax[q] + (1./(2*pi*mu))*(Rp[ic//nc]/2.)*S.w[ic]*(
                        Gxx[ic]*fx[ic//nc] + Gxy[ic]*fy[ic//nc])
                        
                Ax[q + Nq] = Ax[q + Nq] + (1./(2.*pi*mu))*(
                             Rp[ic//nc]/2.)*S.w[ic]*(
                             Gyx[ic]*fx[ic//nc] + Gyy[ic]*fy[ic//nc])

    del X, Y, Nq, Ns, fx, fy   
    
    return Ax
    
    
def vector_Ax_fmm(T,S,nc,f,Rp,Xp,Yp,Re,nchebpol):
    Nq, Ns = len(T.z), len(S.z)

    fx, fy = np.zeros(Nq, dtype = np.float64), np.zeros(Nq, dtype = np.float64)
    fx[:], fy[:] = f[:Nq], f[Nq:]
    
    AxS = np.zeros(2*Nq, dtype = np.float64)
    
    #print
    #print "Inicio calculo de singular en Sigle-layer-potential" 
    
    T_s_i = time()
    #for q in range(Nq):
    #    AxS[q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
    #            int_alpha_xx(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
    #            int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])
                
    #    AxS[Nq + q] = (1./(2*pi*mu))*(  # (Rp[q]/2.)*(
    #            int_alpha_xy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fx[q] + 
    #            int_alpha_yy(Xp[q],Xp[q+1],Yp[q],Yp[q+1])*fy[q])

    AxS[:Nq] = (Re/(2*pi))*(  # (Rp[q]/2.)*(
                int_alpha_xx(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fx[:] + 
                int_alpha_xy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fy[:])

    AxS[Nq:] = (Re/(2*pi))*(  # (Rp[q]/2.)*(
            int_alpha_xy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fx[:] + 
            int_alpha_yy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fy[:])
    
    T_s_f = time()
    #print "Tiempo calculo singularidad: ", T_s_f - T_s_i, "[s]"   
    #print "Fin calculo de singular en Sigle-layer-potential"
    #print
    #print Ax
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)
    XT = np.zeros(Nq)
    YT = np.zeros(Nq)

    Axfmm = np.zeros(2*Nq, dtype = np.float64)


    X[:] = S.z[:].real
    Y[:] = S.z[:].imag  
    XT[:] = T.z[:].real
    YT[:] = T.z[:].imag
    

    bbfmm_SL_potential(X,
        			Y,
        			XT,
        			YT,
        			Rp,
        			S.w[:],
        			fx,
        			fy,
        			Re,
        			Axfmm,
        			nchebpol)
           
    Axfmm[:] = Axfmm[:] + AxS[:]
    
    del X, Y, XT, YT, AxS, fx, fy
    del Nq, Ns
 
    return Axfmm   
    
def vector_Ax_p2p_cy(T,S,nc,f,Rp,Xp,Yp,Re):
    Nq, Ns = len(T.z), len(S.z)

    fx, fy = np.zeros(Nq, dtype = np.float64), np.zeros(Nq, dtype = np.float64)
    fx[:], fy[:] = f[:Nq], f[Nq:]
    
    AxS = np.zeros(2*Nq, dtype = np.float64)
    
    #print
    #print "Inicio calculo de singular en Sigle-layer-potential" 
    
    T_s_i = time()

    AxS[:Nq] = (Re/(2*pi))*(  # (Rp[q]/2.)*(
                int_alpha_xx(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fx[:] + 
                int_alpha_xy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fy[:])

    AxS[Nq:] = (Re/(2*pi))*(  # (Rp[q]/2.)*(
            int_alpha_xy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fx[:] + 
            int_alpha_yy(Xp[:Nq],Xp[1:],Yp[:Nq],Yp[1:])*fy[:])
    
    T_s_f = time()
    #print "Tiempo calculo singularidad: ", T_s_f - T_s_i, "[s]"   
    #print "Fin calculo de singular en Sigle-layer-potential"
    #print
    #print Ax
    
    X = np.zeros(Ns)
    Y = np.zeros(Ns)
    XT = np.zeros(Nq)
    YT = np.zeros(Nq)

    Ax = np.zeros(2*Nq, dtype = np.float64)


    X[:] = S.z[:].real
    Y[:] = S.z[:].imag  
    XT[:] = T.z[:].real
    YT[:] = T.z[:].imag
    

    direct_SL_potential(X,
        			Y,
        			XT,
        			YT,
        			Rp,
        			S.w[:],
        			fx,
        			fy,
        			Re,
        			Ax)
           
    Ax[:] = Ax[:] + AxS[:]
    
    del X, Y, XT, YT, AxS, fx, fy
    del Nq, Ns
 
    return Ax 

def u_field_dir_cy(S, T, Rp, u_0, nx, ny, fx, fy, Re, theta_rot):
    
    Ns, Nt, n_pan = len(S.z), len(T.z), len(Rp)
    
    xsrc = np.zeros(Ns, dtype = np.float64)
    ysrc = np.zeros(Ns, dtype = np.float64)    

    wc = np.zeros(Ns, dtype = np.float64)    

    xtar = np.zeros(Nt, dtype = np.float64)
    ytar = np.zeros(Nt, dtype = np.float64)    
    
    xsrc[:] = S.z[:].real
    ysrc[:] = S.z[:].imag
    
    xtar[:] = T.z[:].real
    ytar[:] = T.z[:].imag
    
    wc[:] = S.w[:]
    
    
    ux, uy = np.zeros(n_pan, dtype=np.float64), np.zeros(n_pan, dtype=np.float64)
    
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0
    
    ux_f, uy_f = np.zeros(Nt, dtype = np.float64), np.zeros(Nt, dtype = np.float64)
    
    ux_f[:] = cos(theta_rot) #*u_0/u_0
    uy_f[:] = sin(theta_rot) #*u_0/u_0
    
    direct_ufield(xsrc,
			ysrc,
			xtar,
			ytar,
			Rp,
			wc,
			ux,
			uy,
			nx,
			ny,
			fx,
			fy,
			Re,
			ux_f,
			uy_f)
    
    del xsrc, ysrc, xtar, ytar, ux, uy
    
    return ux_f, uy_f

def u_field_fmm(S, T, Rp, u_0, nx, ny, fx, fy, Re, theta_rot, ncheb):
    
    Ns, Nt, n_pan = len(S.z), len(T.z), len(Rp)
    
    xsrc = np.zeros(Ns, dtype = np.float64)
    ysrc = np.zeros(Ns, dtype = np.float64)    

    wc = np.zeros(Ns, dtype = np.float64)    

    xtar = np.zeros(Nt, dtype = np.float64)
    ytar = np.zeros(Nt, dtype = np.float64)    
    
    xsrc[:] = S.z[:].real
    ysrc[:] = S.z[:].imag
    
    xtar[:] = T.z[:].real
    ytar[:] = T.z[:].imag
    
    wc[:] = S.w[:]
    
    
    ux, uy = np.zeros(n_pan, dtype=np.float64), np.zeros(n_pan, dtype=np.float64)
    
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0
    
    ux_f, uy_f = np.zeros(Nt, dtype = np.float64), np.zeros(Nt, dtype = np.float64)
    
    ux_f[:] = cos(theta_rot) #*u_0/u_0
    uy_f[:] = sin(theta_rot) #*u_0/u_0
    
    fmm_ufield(xsrc,
			ysrc,
			xtar,
			ytar,
			Rp,
			wc,
			ux,
			uy,
			nx,
			ny,
			fx,
			fy,
			Re,
			ux_f,
			uy_f,
                ncheb)
    
    del xsrc, ysrc, xtar, ytar, ux, uy   
    
    return ux_f, uy_f


def p_field_dir_cy(S, T, Rp, u_0, nx, ny, fx, fy, Re, theta_rot):
    
    Ns, Nt, n_pan = len(S.z), len(T.z), len(Rp)
    
    xsrc = np.zeros(Ns, dtype = np.float64)
    ysrc = np.zeros(Ns, dtype = np.float64)    

    wc = np.zeros(Ns, dtype = np.float64)    

    xtar = np.zeros(Nt, dtype = np.float64)
    ytar = np.zeros(Nt, dtype = np.float64)    
    
    xsrc[:] = S.z[:].real
    ysrc[:] = S.z[:].imag
    
    xtar[:] = T.z[:].real
    ytar[:] = T.z[:].imag
    
    wc[:] = S.w[:]
    
    
    ux, uy = np.zeros(n_pan, dtype=np.float64), np.zeros(n_pan, dtype=np.float64)
    
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0
    
    p_f = np.zeros(Nt, dtype = np.float64)
        
    direct_pfield(xsrc,
			ysrc,
			xtar,
			ytar,
			Rp,
			wc,
			ux,
			uy,
			nx,
			ny,
			fx,
			fy,
			Re,
			p_f)
       
    del xsrc, ysrc, xtar, ytar, ux, uy
    
    
    return p_f
    
def p_field_fmm(S, T, Rp, u_0, nx, ny, fx, fy, Re, theta_rot, ncheb):
    
    Ns, Nt, n_pan = len(S.z), len(T.z), len(Rp)
    
    xsrc = np.zeros(Ns, dtype = np.float64)
    ysrc = np.zeros(Ns, dtype = np.float64)    

    wc = np.zeros(Ns, dtype = np.float64)    

    xtar = np.zeros(Nt, dtype = np.float64)
    ytar = np.zeros(Nt, dtype = np.float64)    
    
    xsrc[:] = S.z[:].real
    ysrc[:] = S.z[:].imag
    
    xtar[:] = T.z[:].real
    ytar[:] = T.z[:].imag
    
    wc[:] = S.w[:]
    
    
    ux, uy = np.zeros(n_pan, dtype=np.float64), np.zeros(n_pan, dtype=np.float64)
    
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0
    
    P_f = np.zeros(Nt, dtype = np.float64)
        
    fmm_pfield(xsrc,
			ysrc,
			xtar,
			ytar,
			Rp,
			wc,
			ux,
			uy,
			nx,
			ny,
			fx,
			fy,
			Re,
			P_f,
                ncheb)
    
    del xsrc, ysrc, xtar, ytar, ux, uy
    
    return P_f

