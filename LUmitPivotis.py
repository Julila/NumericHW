from pylab  import *
import copy
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

#############################################
 
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
############################################
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
        


#EINFACHES BEISPIEL aus dem Skript
"""
print "LU Zerlegung mit Pivotisierung "          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]], dtype=float16)
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])
 
LUP = zerlegungPivot(A)
 
x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)

 
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2 # Die Lösungen stimmen überein
"""

 

#erstelle matrix A für die Aufgabe 4 
def MatrixA(n, beta):
    A = [n*[0] for i in range(n)]
    for j in xrange(0, n):
        if j != n-1:
            A[j][j] = 1
        else:
            A[j][j] = 0
            
        A[j][j-1] = -beta
    A[0][n - 1] = beta
    return A

#erstelle vektor b für die Aufgabe 4 
def VektorB(n, beta):
    b = zeros(n, float)

    for j in xrange(0, n):
        b[j] = float(1 - beta ) 

    b[0] = 1 + beta
    b[n - 1] = -beta
    return b

 

#EINFACHES BEISPIEL
"""
print "LU Zerlegung mit Pivotisierung "          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]], dtype=float16)
b1 = array([3,5,4,5])
b2 = array([4,10,12,11])
 
LUP = zerlegungPivot(A)
 
x1 = permutation(LUP[1],b1)
x1 = vorwaerts(LUP[0],x1)
x1 = rueckwaerts(LUP[0],x1)

 
x2 = permutation(LUP[1],b2)
x2 = vorwaerts(LUP[0],x2)
x2 = rueckwaerts(LUP[0],x2)
print x2
""" 

################
print "Aufgabe 4 " 
Beta = 10
n = [10, 15, 20]

for val in n:
    exacteLoesung =  [1 for i in range(val)]
    print " n = " , val 
    A = MatrixA(val, Beta)
    B = VektorB(val, Beta)
  
    #LU Pivotisierung
    LUP = zerlegungPivot(A)
    yP = vorwaerts(LUP[0],permutation(LUP[1],B))
    xP = rueckwaerts(LUP[0],yP)
    # LU (ohne Pivotisierung)
    LU  = zerlegungPivot(A)
    y = vorwaerts(LU[0],permutation(LU[1],B))
    x = rueckwaerts(LU[0],y)

    print "Die Lösung der LGS mit LU mit Pivotisierung"
    print xP
    print "Die Lösung der LGS mit LU ohne Pivotisierung"
    print x
    print "Exacte Lösung"
    print exacteLoesung
    print "--------------------------------------------"

    # Die Lösungen stimmen ja überein.
    
# Man kann behaupten, dass die Pivot-Verfahren  numerischen Fehler reduzieren. Aber jedoch beim Schriftlichen Lösung der Aufgaben ergibt sich mehr Aufwand. 
#Beim Pivot-Verfahren sucht man  in der  ersten Spalte das betragsmäßig maximale Element, tauscht die Zeilen und richtet danach die Elimination aus.

