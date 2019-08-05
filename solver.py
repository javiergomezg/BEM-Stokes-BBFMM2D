import numpy as np
import gmres
import cilmesh
import matrixfreefunction
from numpy import pi, cos, sin

def f_traction(A,B, ux, uy):
    n = len(ux)
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    
    b = np.matmul(B,u)
    print b
    f = np.linalg.solve(A, b)
        
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n,0], f[n:,0]
        
    return fx, fy
    
    
def f_traction_freeM(A,T, S, nc, ux, uy,Rp,nx,ny):
    n = len(ux)
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    
    b = vector_b_p2p(T,S,nc, ux, uy, Rp, nx, ny)    
    
    #b = np.matmul(B,u)
    #print b
    f = np.linalg.solve(A, b)
    
    #print f
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n,0], f[n:,0]
    
    #print fx
    #print
    #print fy
    
    return fx, fy
    
def f_traction_gmres(A,B, ux, uy):
    n = len(ux)
    nn = 2*n
    u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0E-16
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    b = np.matmul(B,u)
   
    
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)
    
    #print b
    #print b2
    f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-6)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy

def f_traction_gmres_Mfree(T, S, nc,Xp, Yp, ux, uy, Rp, nx, ny, Re):
    n = len(ux)
    nn = 2*n
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([2*n,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 1.0
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    #b = np.matmul(B,u)
    b = vector_b_p2p_cy(T,S,nc, ux, uy, Rp, nx, ny)
   
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)

    
    #print b

    #f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-14)
    f, iteration = gmres_mgs(T,S,nc,f,b2,Rp,Xp,Yp,Re, 200, 200, 1.0E-6)
    
    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy, iteration
    
def f_traction_gmres_Mfree_fmm(T, S, nc,Xp, Yp, ux, uy, Rp, nx, ny, Re,n_cheb):
    n = len(ux)
    nn = 2*n
    #u = np.zeros([2*n,1], dtype = np.float64)
    f = np.zeros([nn,1], dtype = np.float64)
    
    #u[:n,0], u[n:,0] = ux[:], uy[:]
    f = np.zeros(nn, dtype = np.float64)
    f[:] = 0.0
    
    #b = prueba_cont(T,S,nc, ux, uy, Rp,nx, ny)
    #b = np.matmul(B,u)
    b = vector_b_fmm(T,S,nc, ux, uy, Rp, nx, ny, n_cheb)
   
    #D = np.zeros([nn,nn], dtype=np.float64)
    #for i in range(nn):
    #    D[i,i] = A[i,i]**(-1)


    #AA = np.matmul(D,A)
    #bb = np.matmul(D,b)
    b2 = np.transpose(b)

    
    #print b

    #f,_ = gmres_mgs(A,f,b2, nn, 500, 1.0E-14)
    print "Entrando al GMRES"
    print
    f,iteration, time_iteration = gmres_mgs_fmm(T,S,nc,f,b2,Rp,Xp,Yp,Re, n_cheb, 200, 200, 1.0E-6)

    fx = np.zeros(n, dtype = np.float64) 
    fy = np.zeros(n, dtype = np.float64)
    
    fx[:], fy[:] = f[:n], f[n:]
    
    del f
    
    #print fx
    #print
    #print fy
    
    return fx, fy, iteration, time_iteration
    
    
    
def tangential_of_ft(fx,fy, tx, ty, h):
    ftn = np.zeros(len(fx), dtype = np.float64)
    dDrag = np.zeros(len(fx), dtype = np.float64)    
    
    ftn[:] = fx[:]*tx[:] + fy[:]*ty[:]
    
    return ftn
    
def calc_ndrag(fx,h):
    dDrag = np.zeros(len(fx), dtype = np.float64)
    dDrag[:] = fx[:]*h[:]
    
    return sum(dDrag)
    
    
def test_DRAG_dir(n_pan, nq, R_cil, X_c, Y_c, u_0, rho, Re, L):   
    h_theta = 2.0*pi/(n_pan) 
    X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)
    
    X_col, Y_col, X_pan, Y_pan, R_pan = X_col/L, Y_col/L, X_pan/L, Y_pan/L, R_pan/L
    
    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    target = targetvector(n_pan,X_col,Y_col)
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
        
    ux[:] = u_0/u_0
    
##################################################################################    
    ## Calculo directo del DRAG (Sin speed up)
            
    dirT_i = time()
    fx, fy, iteration_dir = f_traction_gmres_Mfree(target, sources, nq, X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, Re)
                    
    NDRAG_direct = calc_ndrag(fx, R_pan)*rho*L*u_0**2

    dirT_f = time()

    dirT = dirT_f - dirT_i
    
    del X_col, Y_col, X_pan, Y_pan, R_pan
    del t, tx_pan, ty_pan, nx_pan, ny_pan
    del target
    del sources
    del uy, ux
    
    return NDRAG_direct, fx, fy, iteration_dir, dirT
    
def test_DRAG_fmm(n_pan, nq, n_cheb, R_cil, X_c, Y_c, u_0, rho, Re, L):  
    
    h_theta = 2.0*pi/(n_pan) 
    X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)
    
    #X_col, Y_col, X_pan, Y_pan, R_pan = X_col/L, Y_col/L, X_pan/L, Y_pan/L, R_pan/L
    
    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    target = targetvector(n_pan,X_col,Y_col)
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
        
    ux[:] = u_0#/u_0
    
##################################################################################    
    ## Calculo directo del DRAG (Sin speed up)
            
    fmmT_i = time()
    fx, fy, iteration_fmm, time_iteration = f_traction_gmres_Mfree_fmm(target, sources, 
                                                                   nq, 
                                                                   X_pan, 
                                                                   Y_pan, 
                                                                   ux, 
                                                                   uy, 
                                                                   R_pan, 
                                                                   nx_pan, 
                                                                   ny_pan, 
                                                                   Re,
                                                                   n_cheb)
                    
    NDRAG_fmm = calc_ndrag(fx, R_pan)

    fmmT_f = time()

    fmmT = fmmT_f - fmmT_i
    
    del X_col, Y_col, X_pan, Y_pan, R_pan
    del t, tx_pan, ty_pan, nx_pan, ny_pan
    del target
    del sources
    del uy, ux
    
    return NDRAG_fmm, fx, fy, iteration_fmm, fmmT
    
def DRAG_test_RBC_dir(d,h,b,r, nq, n_pan_fc, X_c, Y_c, theta_rot, L, u_0, rho, Re):

    X_pan, Y_pan, X_col, Y_col, R_pan, n_pan = RBC_mesh(d,h,b,r, 
                                                n_pan_fc, 
                                                X_c, Y_c)
                                                
    X_pan, Y_pan, X_col, Y_col, R_pan = X_pan/L, Y_pan/L, X_col/L, Y_col/L, R_pan/L                                            
                                                
    #X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    target = targetvector(n_pan,X_col,Y_col)
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
    
    #se que es lo mismo que dejar solo seno y coseno de theta_rot, solo es para esquematizar
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0
    
    
    ##################################################################################    
    ## Calculo directo del DRAG (Sin speed up)
    
    dirT_i = time()
    fx, fy, it_dir = f_traction_gmres_Mfree(target, sources, 
        nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, Re)
    dirT_f = time()
    
    dirT = dirT_f-dirT_i
    
    NDRAG_direct = calc_ndrag(fx, R_pan)*rho*L*u_0**2
    del ux, uy #por ahora voy a borrar, quizas en un futuro hare que las retorne al main para usarlas
    
    return fx, fy, it_dir, NDRAG_direct, n_pan, X_col, Y_col, tx_pan, ty_pan, R_pan, nx_pan, ny_pan, sources, dirT   
    
def DRAG_test_RBC_fmm(d,h,b,r, nq, n_pan_fc, X_c, Y_c, theta_rot, L, u_0, rho, Re, n_cheb):

    X_pan, Y_pan, X_col, Y_col, R_pan, n_pan = RBC_mesh(d,h,b,r, 
                                                n_pan_fc, 
                                                X_c, Y_c)
                                                
    X_pan, Y_pan, X_col, Y_col, R_pan = X_pan/L, Y_pan/L, X_col/L, Y_col/L, R_pan/L                                            
                                                
    #X_col, Y_col, X_pan, Y_pan, R_pan = geometria_cil(R_cil, n_pan, h_theta, X_c, Y_c)    
    t = nodos_gauss_chebyshev(nq)   
    tx_pan, ty_pan, nx_pan, ny_pan = unit_vector1(X_pan, Y_pan, R_pan)
    
    target = targetvector(n_pan,X_col,Y_col)
    
    sources = sourcevector(X_pan, Y_pan, t)
    
    ux = np.zeros(n_pan, dtype = np.float64)
    uy = np.zeros(n_pan, dtype = np.float64)
    
    #se que es lo mismo que dejar solo seno y coseno de theta_rot, solo es para esquematizar
    ux[:] = cos(theta_rot) #*u_0/u_0
    uy[:] = sin(theta_rot) #*u_0/u_0

    ##################################################################################    
    ## Calculo directo del DRAG (Sin speed up)
    
    fmmT_i = time()
    fx, fy, it_fmm, T_it = f_traction_gmres_Mfree_fmm(target, sources, 
        nq,X_pan, Y_pan, ux, uy, R_pan, nx_pan, ny_pan, Re, n_cheb)
    fmmT_f = time()

    fmmT = fmmT_f-fmmT_i    
    
    NDRAG_fmm = calc_ndrag(fx, R_pan)

    del ux, uy #por ahora voy a borrar, quizas en un futuro hare que las retorne al main para usarlas
    
    return fx, fy, it_fmm, NDRAG_fmm, n_pan, X_col, Y_col , tx_pan, ty_pan, R_pan, nx_pan, ny_pan, sources, fmmT   