import cilmesh
import cuadratura
#import matrix_gen
import matrixfreefunction
import Classnode
import solver
import numpy as np
import time
from time import time
from numpy import pi
from numpy import log, sqrt

R_cil = 0.0005   # radio del cilindo [m]
u_0 = 1.    # velocidad [m/s]

L = 2.*R_cil        # Diametro cilindro [m] (longitud caracteristica)

n_pan = [1000]#, 20000]#, 200, 400, 1000, 2000]#, 4000, 10000, 20000, 40000] 

nq = 8       # nodos cuadratura gaussiana
mu = 1.0          
rho = 1.

ni = mu/rho

Re = u_0*(L)/ni

nc_pol = 7 #cantidad de nodos que se utilizan para aproximar el kernel

#R_cil = 5.E-7          # velocidad e-6 [m/s] 
print
print "###########################################################"
print "############PRUEBA: ECUACION DIMENSIONAL##################"
print "###########################################################"
print
print 'Re = ', Re

X_c, Y_c = 0.0 , 0.0 # Centro del cilindro

print

################################################################
########## Formulas Drag de Lamb################################
################################################################

eul_const = 0.577215664901532860606
eps = (0.5-eul_const - log(Re/4.))**(-1)
DRAGL = 4.*pi*mu*u_0*eps

Del1 = (log(4./Re)-eul_const-0.5)**(-1)
Cd_teo = (4.*pi/Re)*(Del1 -0.87*(Del1**3))

#print DRAGL 
#print Cd_teo

#NDRAG_direct = np.zeros(len(n_pan), dtype = np.float64)
#iteration_dir = np.zeros(len(n_pan), dtype = np.int64)
#dirT = np.zeros(len(n_pan), dtype = np.float64)

NDRAG_fmm = np.zeros(len(n_pan), dtype = np.float64)
iteration_fmm = np.zeros(len(n_pan), dtype = np.int64)
fmmT = np.zeros(len(n_pan), dtype = np.float64)

    
for i in range(len(n_pan)):
    
    #NDRAG_direct[i], fx, fy, iteration_dir[i], dirT[i] = test_DRAG_dir(n_pan[i],
    #                                                        nq, R_cil, 
    #                                                        X_c, Y_c, 
    #                                                        u_0, 
    #                                                        rho, Re, L)    
    
    NDRAG_fmm[i], fx_fmm, fy_fmm, iteration_fmm[i], fmmT[i] = test_DRAG_fmm(n_pan[i], 
                                                                            nq, nc_pol, R_cil, 
                                                                            X_c, Y_c, 
                                                                            u_0, rho,
                                                                            (1./mu), L)
    
    del fx_fmm, fy_fmm
    #del fx, fy    
    #del fx, fy, fx_fmm, fy_fmm 

for i in range(len(n_pan)):
    print
    print 'numero de paneles: ', n_pan[i]
    print 'numero de particulas total: ', n_pan[i]*nq + n_pan[i]
 #   print '****Calculo directo******************'
 #   print 'Drag [N/m]= ', NDRAG_direct[i]
 #   print 'Iteraciones = ', iteration_dir[i]
 #   print 'Tiempo [s] = ', dirT[i]
    print '****Calculo con Speed-up******************'
    print 'Drag [N/m]= ', NDRAG_fmm[i]*L/(rho*u_0**2)
    print 'Iteraciones = ', iteration_fmm[i]
    print 'Tiempo [s] = ', fmmT[i]
    
#del NDRAG_direct, iteration_dir, dirT
del NDRAG_fmm, iteration_fmm, fmmT

