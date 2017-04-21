# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 13:23:39 2017

@author: waldemar
"""


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


# LU Zerlegung mit Pivotisierung

def LU(A, b):
    
    n = len(A) # Give us total of lines
    # (1) Permutation Vector
    V = [0 for i in range(n)]
    print V



    # (2) Fill L matrix and its diagonal with 1
    L = [[0 for i in range(n)] for i in range(n)]
    for i in range(0,n):
        L[i][i] = 1

 
    # (3) Fill U matrix
    U = [[0 for i in range(0,n)] for i in range(n)]
    for i in range(0,n):
        for j in range(0,n):
            U[i][j] = A[i][j]


    n = len(U)

    # (4) Find both U and L matrices
    for i in range(0,n): # for i in [0,1,2,..,n]
                
        maxElem = abs(U[i][i]) # Wir nehmen erstes Element, als wäre er der größte
        maxRow = i
        print "U before swapping "
        print U
        
        for k in range(i+1, n): # Interacting over the next line
            if(abs(U[k][i]) > maxElem):
                maxElem = abs(U[k][i]) # Next line on the diagonal
                maxRow = k
        print "V"
        V[i] = maxRow
        print V
        # (4.2) Swap the rows pivoting the maxRow, i is the current row
        for k in range(i, n): # Interacting column by column
            tmp=U[maxRow][k]
            U[maxRow][k]=U[i][k]
            U[i][k]=tmp

        # (4.3) Subtract lines
        for k in range(i+1,n):
            c = -U[k][i]/float(U[i][i])
            L[k][i] = c # (4.4) Store the multiplier
            for j in range(i, n):
                U[k][j] += c*U[i][j] # Multiply with the pivot line and subtract
        print "L "
        print L
        # (4.5) Make the rows bellow this one zero in the current column
        for k in range(i+1, n):
            U[k][i]=0

    n = len(L)

    # (5) Perform substitutioan Ly=b
    y = [0 for i in range(n)]
    for i in range(0,n,1):
        y[i] = b[i]/float(L[i][i])
        for k in range(0,i,1):
            y[i] -= y[k]*L[i][k]

    n = len(U)

    # (6) Perform substitution Ux=y
    x = [0 in range(n)]
    """for i in range(n-1,-1,-1):
        x[i] = y[i]/float(U[i][i])
        for k in range (i-1,-1,-1):
            U[i] -= x[i]*U[i][k]"""

    return x

    






print "LU Zerlegung"          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])


LUP = zerlegung(A)
print LUP[0]
"""
x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)
print x1
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2"""

print "Der Unterschied liegt daran, dass man ..."


print "LU Zerlegung mit Spalten-Pivot-Suche"          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])
print "Start matrix"
print A
LUP = LU(A,b1)
print LUP
"""x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)
print x1
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2
"""



"""

for n in range (5,21,5):
    
    bn = erstelleB(n)
    An = erstelleA(n)
    LUPn = zerlegung(An)
    x = permutation(LUPn[1],bn)
    x = vorwaerts(LUPn[0],x)
    x = rueckwaerts(LUPn[0],x)
    print "n={}".format(n)
    print "x={}".format(x)
 
    
u = array([0,1,2,3])
v = array([0,0,0,1])

As = A + outer(u,v)
xd = shermanMorrison(LUP,b1,u,v)
print As
print xd
print dot(As, xd)


"""
    
    




            

                
        
