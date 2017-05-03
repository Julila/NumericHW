from pylab  import *
import numpy as np
import sys
import copy
import numpy as np
from numpy.linalg import cond
from scipy.linalg import norm, inv
#########LU################################

def zerlegung(A):
    A = array(A)
    dim = A.shape[0]
    p = zeros(dim, int)
    
    for i in range(0,dim):
        h = i;
        if(A[i, i]==0): # Falls Element auf Position [i,i] = 0 ist, dann tauschen wir die Zeilen
            h = find(A[i:,i]!=0)[0]+piv;
            A[[i,h],:] =A[[h,i],:] #Zeilen tauschen
            
        p[i]=h
        for zeile in range(i+1,dim):
            koeff=float(A[zeile,i])/A[i, i] # Finden L11
            substZeile = multiply(koeff, A[i,i+1:]) # Koeff mal neue Zeile
            A[zeile,i+1:] = subtract(A[zeile,i+1:],substZeile)# Sustrairen von der alte Zeile die neue zeile
            A[zeile,i]=koeff # Speichern L11 in A matrix
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

#erstelle matrix A für die Aufgabe 5
def MatrixA( sigma):
    A = array([[3,2,1], [2, 2*10**sigma, 2*10**sigma ], [1, 2*10**sigma, -10**sigma ]])
 
    return A

#erstelle vektor b für die Aufgabe 5 
def VektorB(sigma):
    b = array([3+ 3*(10**sigma),6*(10**sigma),2*(10**sigma)])
    return b        
#############################################


 
print "Aufgabe 5 " 
print "Satz von Prager und Oetti"
print "------------------"

Beta = 10
n = [-8,-10,-12]

for val in n:
    exacteLoesung =   array([ 10**val,1,1])
    
    print " n = " , 10**val
    print "Exacte Loesung : ", exacteLoesung

    A = MatrixA(val)
    B = VektorB(val)

    #Suchen xstange
    LU = zerlegung(A)
    pb = permutation(LU[1],B)
    y = vorwaerts(LU[0],pb)
    xstange = rueckwaerts(LU[0],y)
    
    sstange = (np.abs(A)).dot(np.abs(xstange)) + np.abs(B)
    #Satz Routine 
    rstange = B - A.dot(xstange)
    eps =  sys.float_info.epsilon  # aus sys bibliothek 
    epsilon = np.max(np.divide(np.abs(rstange), sstange))
    k = np.linalg.cond(A)
    print "Kondition A (lilang bib)  : K(a)*eps:" , eps * k

    print "Kondition A   : K(a)*eps: ", eps * norm(A, np.inf) * norm(inv(A), np.inf),
    print  '\n'


###
# Lösung: Obwohl sigma immer kleiner wird, vergrößern sich Kondition und epsilon Zahl
###







