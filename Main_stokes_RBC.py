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
from matplotlib import pyplot, cm
from matplotlib import rcParams

""""
Parametros:

d: Diametro principal del globulo rojo
h: Altura mayor
b: altura menor
r: radio interno
"""""""""

d = 7.82
h = 2.58
b = 0.81
r = d/3
a = 1.

L = h

N_cont = 21

#d = 10.
#h = 4.
#b = 2.
#r = d/3

n_pan_fc = 3200#, 80, 160]

#centro del eritrocito

X_c = 0.0
Y_c = 0.0

theta_rot = 0.0
theta_rot = theta_rot*pi/180

u_0 = 0.5E6

#n_pan = [100, 200, 400, 1000, 2000, 4000, 10000]#, 20000, 40000] 

nq = 9        # nodos cuadratura gaussiana
mu = 1.5E-9         
rho = 1.0540E-15

ni = mu/rho

Re = u_0*(L)/ni

nc_pol = 8 #cantidad de nodos que se utilizan para aproximar el kernel

#R_cil = 5.E-7         # radio e-6 [m]
#u_0 = 1.E-6            # velocidad e-6 [m/s] 
print
print "###########################################################"
print "############PRUEBA: ECUACION ADIMENSIONAL##################"
print "###########################################################"
print
print 'Re = ', Re


    
#fx, fy, it_dir, NDRAG_direct, n_pan, X_col, Y_col, tx, ty, Rp, nx_pan, ny_pan, sources, dirT = DRAG_test_RBC_dir(d,h,b,r, 
#                                                                                           nq, n_pan_fc, 
#                                                                                           X_c, Y_c, 
#                                                                                           theta_rot, 
#                                                                                           L, u_0, 
#                                                                                           rho, Re)
                                                                                        

fx_fmm, fy_fmm, it_fmm, NDRAG_fmm, n_pan, X_col, Y_col, tx, ty, Rp, nx_pan, ny_pan, sources, fmmT = DRAG_test_RBC_fmm(d,h,b,r, 
                                                                                              nq, n_pan_fc, 
                                                                                              X_c, Y_c, 
                                                                                              theta_rot, 
                                                                                              L, u_0, 
                                                                                              rho, Re, nc_pol)
                                                                      
x_cont, y_cont = contourn_RBC(d, a, h, b, r, N_cont)

x_cont, y_cont = x_cont/L, y_cont/L


N_cont_p = 81

x_cont_p, y_cont_p = contourn_RBC(d, a, h, b, r, N_cont_p)

x_cont_p, y_cont_p = x_cont_p/L, y_cont_p/L


target = targetvector(len(x_cont), x_cont, y_cont) 

target_p = targetvector(len(x_cont_p), x_cont_p, y_cont_p)


#ux_f, uy_f = u_field_dir_cy(sources, target, Rp, u_0, nx_pan, ny_pan, fx, fy, Re, L, theta_rot)
#ux_f_fmm, uy_f_fmm = u_field_fmm(sources, target, Rp, u_0, nx_pan, ny_pan, fx_fmm, fy_fmm, Re, L, theta_rot,nc_pol)

p_f = p_field_dir_cy(sources, target_p, Rp, u_0, nx_pan, ny_pan, fx, fy, Re, L, theta_rot)
p_f_fmm = p_field_fmm(sources, target_p, Rp, u_0, nx_pan, ny_pan, fx_fmm, fy_fmm, Re, L, theta_rot, nc_pol)

#print target.z[:]

print 'Numero de Paneles: ', n_pan
#print 'iteraciones directas: ', it_dir
print 'iteraciones fmm: ', it_fmm
#print 'Drag directo= ',-NDRAG_direct
#print 'tiempo directo [s]= ', dirT
print 'Drag FMM= ',-NDRAG_fmm
print 'tiempo fmm [s]=', fmmT
#print 'Desviacion= ', abs(NDRAG_direct-NDRAG_fmm)


#ft_direct = tangential_of_ft(fx,fy, tx, ty, Rp)*(rho*u_0**2)*1E6
ft_fmm = tangential_of_ft(fx_fmm,fy_fmm, tx, ty, Rp)*(rho*u_0**2)*1E6

X_col, Y_col = X_col*L, Y_col*L

x_cont, y_cont = x_cont*L, y_cont*L 

#ux_f, uy_f = -ux_f*1E-1, -uy_f*1E-1
ux_f_fmm, uy_f_fmm = -ux_f_fmm*1E-1, -uy_f_fmm*1E-1


#C = np.zeros(len(ux_f), dtype=np.float64)

#C[:] = (ux_f[:]**2 + uy_f[:]**2)**0.5

C_fmm = np.zeros(len(ux_f_fmm), dtype=np.float64)

C_fmm[:] = (ux_f_fmm[:]**2 + uy_f_fmm[:]**2)**0.5


x_cont_p, y_cont_p = x_cont_p*L, y_cont_p*L 
#p_f = (-p_f*rho*u_0**2.)
p_f_fmm = (-p_f_fmm*rho*u_0**2.)


#fig1 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
#pyplot.title("Esfuersos de corte en Eritrocito (calculo directo)", fontsize = 10)
#pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_direct[:]), alpha = 0.5,
#       edgecolors = "face")
#pyplot.quiver(x_cont, y_cont, ux_f, uy_f)
#pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
#pyplot.xlabel("x [$\\mu$m]")
#pyplot.ylabel("y [$\\mu$m]")
#pyplot.grid(True)
#pyplot.axis([-1.05*d/2 + X_c,1.05*d/2 + X_c,-1.05*d/2 + Y_c,1.05*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('cell_direct_calc.png')
#pyplot.show()    

 


fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
pyplot.title("Esfuerzos de corte en Eritrocito calculo con speed up", fontsize = 10)
pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_fmm[:]), marker = ",", alpha = 0.5,
       edgecolors = "face")
#pyplot.quiver(x_cont, y_cont, ux_f, uy_f)
pyplot.xlabel("x [$\\mu$m]")
pyplot.ylabel("y [$\\mu$m]")     
pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
pyplot.grid(True)
#pyplot.axis([-1.1*d/2 + X_c,1.1*d/2 + X_c,-1.1*d/2 + Y_c,1.1*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
pyplot.savefig('cell_fmm_calc.png')
pyplot.show()      



#fig3 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
#pyplot.title("campo velocidad calculo directo", fontsize = 10)
#pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_fmm[:]), marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.quiver(x_cont, y_cont, ux_f, uy_f, C)
#pyplot.scatter(x_cont_p[:], y_cont_p[:], c=p_f[:], marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.contour(x_cont, y_cont, p_f)
#pyplot.colorbar().set_label('velocidad [m/s]')
#pyplot.xlabel("x [$\\mu$m]")
#pyplot.ylabel("y [$\\mu$m]")    
#pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
#pyplot.grid(True)
#pyplot.axis([-1.1*d/2 + X_c,1.1*d/2 + X_c,-1.1*d/2 + Y_c,1.1*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('camp_vel_calc_dir.png')
#pyplot.show()     



fig4 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
pyplot.title("campo velocidad calculo con speed up", fontsize = 10)
#pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_fmm[:]), marker = ",", alpha = 0.5,
#       edgecolors = "face")
pyplot.quiver(x_cont, y_cont, ux_f_fmm, uy_f, C_fmm)
#pyplot.scatter(x_cont_p[:], y_cont_p[:], c=p_f[:], marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.contour(x_cont, y_cont, p_f)
pyplot.colorbar().set_label('velocidad fmm [m/s]')
pyplot.xlabel("x [$\\mu$m]")
pyplot.ylabel("y [$\\mu$m]")    
#pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
pyplot.grid(True)
#pyplot.axis([-1.1*d/2 + X_c,1.1*d/2 + X_c,-1.1*d/2 + Y_c,1.1*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
pyplot.savefig('camp_vel_calc_fmm.png')
pyplot.show()     



#fig5 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
#         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
#pyplot.title("campo presiones calculo calculo directo", fontsize = 10)
#pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_fmm[:]), marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.quiver(x_cont, y_cont, ux_f, uy_f)
#pyplot.scatter(x_cont_p[:], y_cont_p[:], c=p_f_fmm[:], marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.contour(x_cont, y_cont, p_f)
#pyplot.colorbar().set_label('presion [N/$m^2$]')
#pyplot.xlabel("x [$\\mu$m]")
#pyplot.ylabel("y [$\\mu$m]")    
#pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
#pyplot.grid(True)
#pyplot.axis([-1.1*d/2 + X_c,1.1*d/2 + X_c,-1.1*d/2 + Y_c,1.1*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#pyplot.savefig('camp_press_calc_dir.png')
#pyplot.show()     


fig6 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
         facecolor = 'w', edgecolor = 'k')
#pyplot.plot(X_r,Y_r, 'go--', linewidth = 3)
#pyplot.plot(X_c,Y_c, 'go', linewidth = 3)
pyplot.title("campo presiones calculo con speed up", fontsize = 10)
#pyplot.scatter(X_col[:], Y_col[:], c=abs(ft_fmm[:]), marker = ",", alpha = 0.5,
#       edgecolors = "face")
#pyplot.quiver(x_cont, y_cont, ux_f, uy_f)
pyplot.scatter(x_cont_p[:], y_cont_p[:], c=p_f_fmm[:], marker = ",", alpha = 0.5,
       edgecolors = "face")
#pyplot.contour(x_cont, y_cont, p_f)
pyplot.colorbar().set_label('presion fmm [N/$m^2$]')
pyplot.xlabel("x [$\\mu$m]")
pyplot.ylabel("y [$\\mu$m]")    
#pyplot.colorbar().set_label('Esfuerzos de corte [N/$m^2$]')
pyplot.grid(True)
#pyplot.axis([-1.1*d/2 + X_c,1.1*d/2 + X_c,-1.1*d/2 + Y_c,1.1*d/2 + Y_c])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
pyplot.savefig('camp_press_calc_fmm.png')
pyplot.show()     
##################################################################################  

