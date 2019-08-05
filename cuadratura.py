import numpy as np
#import gauss_legendre_terms
#import fs_gf 
from numpy import pi, sin, cos, log
from matplotlib import pyplot, cm
from matplotlib import rcParams

#### nodos de cuadratura en t -1 a 1
#### Cuadratura GAUSS CHEBYSHEV
def nodos_gauss_chebyshev(n_gsch):
    i = np.arange(1., n_gsch + 1., dtype=np.float64)
    t = np.zeros(n_gsch, dtype = np.float64)
    t[:] = cos(pi*(2.*i[:] -1 )/(2.*n_gsch))
    return t

#### Cuadratura GAUSS Legendre  
def nodos_gauss_legendre(NQ):
    w, t = gauss_leg_coef(NQ)
    return w, t

#### evaluacion de t en X e Y
def nodos_gauss_panel(t, X1, Y1, X2, Y2):
    Xg = np.zeros(len(t), dtype = np.float64)
    Yg = np.zeros(len(t), dtype = np.float64)
    
    Xg[:] = 0.5*(X2 + X1) + 0.5*(X2 - X1)*t[:]
    Yg[:] = 0.5*(Y2 + Y1) + 0.5*(Y2 - Y1)*t[:]
    
    return Xg, Yg

### Integral gauss chebichev para las componentes de la matris B
    
def beta_comp_gs(t,Xg,Yg, Xq, Yq, nx, ny, R):
    Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy = gf_stress_fs(Xg,Yg,Xq,Yq)
    
    bxxaux = np.zeros(len(t), dtype=np.float64)
    bxyaux = np.zeros(len(t), dtype=np.float64)
    byxaux = np.zeros(len(t), dtype=np.float64)
    byyaux = np.zeros(len(t), dtype=np.float64)

    bxxaux[:] = sqrt(1.-t[:]**2)*(Txxx[:]*nx + Txxy[:]*ny)
    bxyaux[:] = sqrt(1.-t[:]**2)*(Txyx[:]*nx + Txyy[:]*ny)
    byxaux[:] = sqrt(1.-t[:]**2)*(Tyxx[:]*nx + Tyxy[:]*ny)
    byyaux[:] = sqrt(1.-t[:]**2)*(Tyyx[:]*nx + Tyyy[:]*ny)

    bxx = 0.5*pi*R*sum(bxxaux)/len(t)
    bxy = 0.5*pi*R*sum(bxyaux)/len(t)
    byx = 0.5*pi*R*sum(byxaux)/len(t)
    byy = 0.5*pi*R*sum(byyaux)/len(t)
    
    return bxx, bxy, byx, byy

### Integral gauss Legendre para las componentes de la matris B
def beta_comp_gl(w, t,Xg,Yg, Xq, Yq, nx, ny, R):
    Txxx, Txxy, Tyxx, Tyxy, Txyx, Txyy, Tyyx, Tyyy = gf_stress_fs(Xg,Yg,Xq,Yq)
    
    bxxaux = np.zeros(len(t), dtype=np.float64)
    bxyaux = np.zeros(len(t), dtype=np.float64)
    byxaux = np.zeros(len(t), dtype=np.float64)
    byyaux = np.zeros(len(t), dtype=np.float64)

    bxxaux[:] = w[:]*(Txxx[:]*nx + Txxy[:]*ny)
    bxyaux[:] = w[:]*(Txyx[:]*nx + Txyy[:]*ny)
    byxaux[:] = w[:]*(Tyxx[:]*nx + Tyxy[:]*ny)
    byyaux[:] = w[:]*(Tyyx[:]*nx + Tyyy[:]*ny)

    bxx = 0.5*R*sum(bxxaux)
    bxy = 0.5*R*sum(bxyaux)
    byx = 0.5*R*sum(byxaux)
    byy = 0.5*R*sum(byyaux)
    
    return bxx, bxy, byx, byy
    
### integrales en paneles sigulares de matris A
    
def int_alpha_xx(X1,X2,Y1,Y2):
    l2 = (X2-X1)**2  + (Y2-Y1)**2
    ll = sqrt(l2)    
    
    I1 = - ll*(log(ll/2.)-1.)
    I2 = ((X2-X1)**2)/ll
    return I1 + I2 
    
def int_alpha_xy(X1,X2,Y1,Y2):
    l2 = (X2-X1)**2  + (Y2-Y1)**2
    ll = sqrt(l2)
    return ((X2-X1)*(Y2-Y1))/ll

def int_alpha_yy(X1,X2,Y1,Y2):
    l2 = (X2-X1)**2  + (Y2-Y1)**2
    ll = sqrt(l2)    
    
    I1 = - ll*(log(ll/2.)-1.)
    I2 = ((Y2-Y1)**2)/ll

    return I1 + I2
    
### integrales en paneles no sigulares de matris A
### Gauss CHEBYSHEV
def alpha_comp_gs(t,Xg,Yg, Xq, Yq, R):
    Gxx, Gxy, Gyx, Gyy = gf_vel_fs(Xg,Yg,Xq,Yq)
    
    axxaux = np.zeros(len(t), dtype=np.float64)
    axyaux = np.zeros(len(t), dtype=np.float64)
    ayxaux = np.zeros(len(t), dtype=np.float64)
    ayyaux = np.zeros(len(t), dtype=np.float64)
    
    axxaux[:] = sqrt(1.-t[:]**2)*Gxx[:]
    axyaux[:] = sqrt(1.-t[:]**2)*Gxy[:]
    ayxaux[:] = sqrt(1.-t[:]**2)*Gyx[:]
    ayyaux[:] = sqrt(1.-t[:]**2)*Gyy[:]
    
    axx = 0.5*pi*R*sum(axxaux)/len(t)
    axy = 0.5*pi*R*sum(axyaux)/len(t)
    ayx = 0.5*pi*R*sum(ayxaux)/len(t)
    ayy = 0.5*pi*R*sum(ayyaux)/len(t)
    
    return axx, axy, ayx, ayy      

### Gauss Legendre    
def alpha_comp_gl(w,t,Xg,Yg, Xq, Yq, R):
    Gxx, Gxy, Gyx, Gyy = gf_vel_fs(Xg,Yg,Xq,Yq)
    
    axxaux = np.zeros(len(t), dtype=np.float64)
    axyaux = np.zeros(len(t), dtype=np.float64)
    ayxaux = np.zeros(len(t), dtype=np.float64)
    ayyaux = np.zeros(len(t), dtype=np.float64)
    
    axxaux[:] = w[:]*Gxx[:]
    axyaux[:] = w[:]*Gxy[:]
    ayxaux[:] = w[:]*Gyx[:]
    ayyaux[:] = w[:]*Gyy[:]
    
    axx = 0.5*R*sum(axxaux)
    axy = 0.5*R*sum(axyaux)
    ayx = 0.5*R*sum(ayxaux)
    ayy = 0.5*R*sum(ayyaux)
    
    return axx, axy, ayx, ayy      

