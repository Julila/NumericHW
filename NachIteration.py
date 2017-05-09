# -*- coding: utf-8 -*-
"""
Created on Sun May 07 19:45:46 2017

@author: yuliia
"""
from pylab  import *
from scipy.linalg import lu
from numpy import eye, mat, tril_indices_from, ones, diag_indices_from, zeros,\
     arange, argmax, triu, tril
from numpy.linalg.linalg import norm, inv

def zerlegungPivot(A):
    A = array(A)
    dim = A.shape[0]
    p = zeros(dim-1, int)
    
    for i in range(0,dim-1):
        # Suchen max element
        maxZeile = i;
        maxElem = abs(A[i][i])
        ###############
        #laufen jede Zeile
        for k in range( i+1, dim):
            if( abs(A[k][i]) > maxElem):
                maxElem = abs(A[k][i])
                maxZeile = k
        # Die größte Zeil gefunden
        
        A[[i,maxZeile],:] =A[[maxZeile,i],:]

        p[i]=maxZeile
        ########
        for zeile in range(i+1,dim): #erste  Zeil bleibt unveränderlich
            koeff = A[zeile][i]/float(A[i][i])
            
            substZeile = multiply(koeff, A[i,i+1:])  # Koeff mal neue Zeile
            A[zeile,i+1:] = subtract(A[zeile,i+1:],substZeile) # Sustrairen von der alte Zeile die neue zeile
            A[zeile,i]=float(koeff)

             
    return [A,p]    
     
def permutation(p,x):
    x = array(x)
    dim = p.shape[0]
    for i in range(0,dim):
        x[[i,p[i]]]=x[[p[i],i]]
    return x

def rueckwaerts(LU,x):
    dim = LU.shape[0]
    y=zeros(dim)
    y[dim-1]=float(x[dim-1])/LU[dim-1,dim-1]
    for i in range(2,dim+1):
        z=dim-i;
        h = dot(y[z+1:dim], LU[z,z+1:dim])
        y[z] = (x[z]-h)/float(LU[z,z])
    return y

def vorwaerts(LU,x):
    dim = LU.shape[0]
    y=zeros(dim)
    y[0]=float(x[0])
    for i in range(1,dim):
        h = dot(y[:i], LU[i,:i])
        y[i] = x[i]-h        
    return y     
        
        
############################################################################++
#Start der  Aufgabe 4
############################################################################
        
        
 
def nachiteration(L, U, b, xk,p):
    pk=2*xk# just to make it start
    while (norm(pk,2)/norm(xk,2)>10**-6):
        Ak=L.dot(U)
        b=permutation(p, b)
        rk=b-Ak.dot(xk)
        pk=(inv(U).dot(inv(L))).dot(rk)
        xk+=pk
    return xk

def erstelleA(n):
    A=mat(eye(n))
    inds=tril_indices_from(A)
    A[inds]=ones((n**2+n)/2)*(-1)
    A[diag_indices_from(A)]=ones(n)
    A.T[n-1]=ones(n)   
    return A  

def erstelleB(n):
    b=arange(2.,2.-n,-1.)
    b[n-1]-=1
    bt = b.T
    return bt

def erstelleA2(n):
    A=mat(eye(n))
    for j in range(n):
        for i in range(j+1,n):
            A[i,j]=i+j+2
    return A

def erstelleB2(n):
    b=zeros(n)
    b[0]=1
    return b

############################################################################
# TESTEN
###########################################################################
 

print " Teil a.)"
for n  in [40,50]: # ( 60 hinzufügen) Die Berechnung dauert zu  lange 
    A= erstelleA(n)
    b = erstelleB(n)
    # Einfache LU Zerlegung mit Spalte Pivotis.

    (LU,p)=zerlegungPivot(A.copy())
    Lb = permutation(p, b)
    Ux = vorwaerts(LU, Lb)
    x = rueckwaerts(LU, Ux)
 
    print "LU Zerlegung mit Spatl.piv. n=",n
    print x
    print "-------------------------------------------------------"
    
    #xk=ones(len(A))*2000
    L=tril(LU).copy()
    L[diag_indices_from(L)]=ones(len(LU))
    U=triu(LU)
    xNach=nachiteration(L, U, b, x,p)
    print "Nachiteration verfahren n=",n
    
    print xNach
    print "-------------------------------------------------------"


 
    
 
print " Teil b.)"
for n in [5,7,10,12]: # Man sollte 15 auch testen, aber es braucht zu lange um durchzulaufen.!!! 
    A=erstelleA2(n)
    b= erstelleB2(n)
    (LU,p)=zerlegungPivot(A.copy())
    Lb = permutation(p, b)
    Ux = vorwaerts(LU, Lb)
    xk = rueckwaerts(LU, Ux)
 
    print "Lösung der LU Zerlegung mit Spatl.piv. n=",n
    print xk
    print "---------------------------------------"
    L=tril(LU).copy()
    L[diag_indices_from(L)]=ones(len(LU))
    U=triu(LU)
    xk=nachiteration(L, U, b, xk,p)
    print "Nachiteration verfahren n=",n
    print xk
    print "---------------------------------------"

       
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
