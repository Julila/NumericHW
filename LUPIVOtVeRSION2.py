from pylab  import *
import copy


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
            
            zh = multiply(koeff, A[i,i+1:])
            A[zeile,i+1:] = subtract(A[zeile,i+1:],zh)
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
        
 

print "LU Zerlegung mit Pivotisierung "          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]], dtype=float16)
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])


LUP = zerlegungPivot(A)

 
 
 
x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)

"""
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2
"""

 

#Matrix A
def MatrixA(n, Beta):
    A = [n*[0] for i in range(n)]
    for j in xrange(0, n):
        if j != n-1:
            A[j][j] = 1
        else:
            A[j][j] = 0
            
        A[j][j-1] = -Beta
    A[0][n - 1] = Beta
    return A

#Vektor B
def VektorB(n, Beta):
    b = zeros(n, float)

    for j in xrange(0, n):
        b[j] = float(1 - Beta ) 

    b[0] = 1 + Beta
    b[n - 1] = -Beta
    return b



 

################
print "Aufgabe 4 " 
Beta = 10
n = [10, 15, 20]

for val in n:
    print " n = " , val 
    A = MatrixA(val, Beta)
    B = VektorB(val, Beta)
 
    LU = zerlegungPivot(A)
 
    y = vorwaerts(LU[0],permutation(LU[1],B))
    x = rueckwaerts(LU[0],y)
    print x

