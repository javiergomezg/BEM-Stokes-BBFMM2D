import numpy as np 
from numpy import pi, sin, cos, arctan2
from math import atan2
from matplotlib import pyplot, cm
from matplotlib import rcParams

a = 4.
b = 2.

n_pan_sup = 41

n_pan = 2*n_pan_sup

x_in = -a
x_out = a

x_sup = np.zeros(n_pan_sup + 1, dtype = np.float64)
y_sup = np.zeros(n_pan_sup + 1, dtype = np.float64)

x_sup[0] = x_in
y_sup[0] = b*((1. - (x_sup[0]**2/a**2))**0.5)

dx = 2*a/n_pan_sup

print 0, x_sup[0], y_sup[0]

for i in range(1,n_pan_sup):
    x_sup[i] = x_sup[i-1]+dx
    y_sup[i] = b*((1. - (x_sup[i]**2/a**2))**0.5)
    print i, x_sup[i], y_sup[i]
    
x_sup[-1] = -x_sup[0]    
y_sup[-1] = y_sup[0]
    
#print x_sup



X_pan = np.zeros(n_pan+1, dtype = np.float64)
Y_pan = np.zeros(n_pan+1, dtype = np.float64)

X_pan[:n_pan_sup+1] = x_sup[:]
Y_pan[:n_pan_sup+1] = y_sup[:]

#print X_pan

X_pan[n_pan_sup+1:-1] = -x_sup[1:-1]
Y_pan[n_pan_sup+1:-1] = -y_sup[1:-1]

X_pan[-1] = x_sup[0]
Y_pan[-1] = y_sup[0]

X_col = np.zeros(n_pan, dtype = np.float64)
Y_col = np.zeros(n_pan, dtype = np.float64)

X_col[:] = 0.5*(X_pan[:-1] + X_pan[1:])

Y_col[:] = 0.5*(Y_pan[:-1] + Y_pan[1:])


#print X_pan
#print
#print Y_pan


#print y_sup



fig2 = pyplot.figure(num = None, figsize = (8, 6), dpi = 80, 
                 facecolor = 'w', edgecolor = 'k')
pyplot.plot(X_pan,Y_pan, 'go--', linewidth = 3)
pyplot.plot(X_col,Y_col, 'ro', linewidth = 3)
#    pyplot.plot(XX,YY, 'bo-', linewidth = 3)
pyplot.grid(True)
pyplot.axis([-5,5,-5,5])
#    pyplot.grid(color = '0.3', linestyle = '--', linewidth = 0.3)
#    pyplot.savefig('circulito.png')
pyplot.show()