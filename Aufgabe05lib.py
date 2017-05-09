# -*- coding: utf-8 -*-
"""
Gleichungssystem für Blatt 5, Aufgabe 5
"""
from pylab import ones, linspace, meshgrid, sqrt, figure, cm, colorbar
from mpl_toolkits.mplot3d import Axes3D

from scipy.sparse import spdiags


def system(m):
    n = m*m
    
    e = ones(n)
    l = ones(n)
    l[m-1::m] = 0.0
    r = ones(n)
    r[::m] = 0.0
    
    A = spdiags([-e, -l, 4.0*e, -r, -e], [-m, -1, 0, 1, m], n, n, format='csr')
    
    b = -e / float(n)
    
    return A, b


def plotxk(xk):
    n = len(xk)
    m = int(sqrt(n))
    
    print m
    
    h = linspace(0, 1, m)
    yy,xx = meshgrid(h,h)
    
    fig = figure('xk, m = {0}'.format(m))
    ax  = fig.gca(projection='3d')
    
    surf = ax.plot_surface(xx, yy, xk.reshape(m,m), cmap = cm.jet, rstride = 5, cstride = 5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("Hoehe")
    colorbar(surf)
    

