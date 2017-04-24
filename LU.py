

from pylab  import *


def zerlegung(A):
    A = array(A)
    dim = A.shape[0]
    p = zeros(dim, int)
   
    for piv in range(0,dim):
        h = piv;
        if(A[piv, piv]==0):
            h = find(A[piv:,piv]!=0)[0]+piv;
            A[[piv,h],:] =A[[h,piv],:]
            
        p[piv]=h
        for zeile in range(piv+1,dim):
            div=float(A[zeile,piv])/A[piv, piv]
            zh = multiply(div, A[piv,piv+1:])
            A[zeile,piv+1:] = subtract(A[zeile,piv+1:],zh)
            A[zeile,piv]=div
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
        
def erstelleA(n):
    A = zeros([n,n])
    for i in range(1,n+1):
        for j in range(1,n+1):
            A[i-1,j-1]=1.0/(i+j-1)
    return A

def erstelleB(n):
    b = zeros(n)
    for i in range(1,n+1):
        b[i-1]=1.0/(i+1)
    return b





print "LU Zerlegung"          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])


LUP = zerlegung(A)
print LUP[0]
print erstelleA(3)
x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)
print x1
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2

