# -*- coding: utf-8 -*-
"""
Created on Tue May 02 11:17:34 2017

@author: yuliia
"""

#encoding: utf-8
import numpy

from pylab  import *
import numpy as np

#############################################
#LU mit Pivot.  
#############################################
def zerlegungPivot(AA):
    A = array(AA)

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
         
###################################################
         
        
###################################################
# Hausholder
###################################################
def householderQR(A):
    # Aus der Aufgabe 3 folgt:
    # d = - sign(a11) norm(a1)
    
    #A transpose
    Atransp = A.transpose() # Tipp aus der Aufgabe 4
    dim = len(A[0])
    hauptDiagonal = numpy.full(dim, 0.)
    
    for i in numpy.arange(dim-1):  
        v = Atransp[i][i:] # Erste Spalte
        
        hauptDiagonal[i] = -numpy.sign(v[0]) * numpy.linalg.norm(v) # ||a||2 Euklidische Norm
 
        v[0] -= hauptDiagonal[i]
        v2 = -2*v[0]*hauptDiagonal[i] #||v2||2 
        
        for j in numpy.arange(i+1,dim):  
            Qi = Atransp[j][i:]  
            Qi -= 2*Qi.dot(v)/v2*v

         
    hauptDiagonal[dim-1] = Atransp[dim-1][dim-1]
    
    return hauptDiagonal

# Berechnen y
def berechneY(A, diag, B):
    Atrans = A.transpose()
    dim = len(B)
    y = B.copy()
     
    for i in numpy.arange(dim-1) :
        v = numpy.zeros(dim)
        v[i:] = Atrans[i][i:] 
        v2 = -2*v[i]*diag[i]  
        y -= 2*y.dot(v)/v2*v  
    return y



def rueckwaerts_HH(A, diag, Y):
    dim =A.shape[0]
    x = zeros(dim)
    
    for i in numpy.arange(len(Y)-1, 0-1, -1):
        
        x[i] = float(Y[i])/float(diag[i]) # Koeffizient 
        for j in numpy.arange(len(Y)-1, i, -1): 
            x[i] -= A[i][j]/float(diag[i])*x[j]

    return x
 



#Erstelle Matrix A für die Audfabe 5
def erstelleMatrixA(n):
    A =  numpy.eye(n)
    for zeile in range(0, n):
        for spalte in range(0, n):
            if(zeile==spalte):
                A[zeile][spalte]=1
                
            elif (spalte==(n-1)):
                A[zeile][spalte]=1
                
            elif (spalte<zeile):
                A[zeile][spalte]=(-1)
                
            elif (zeile>spalte):
                A[zeile][spalte]=0;
                
    return A;


 
#Erstelle Vektor B für die Audfabe 5
def erstelleB(n):
    b = numpy.full(n, 2.)
    for i in numpy.arange(n-1):
        b[i] -= i
    b[n-1] = 2-n
    return b






##############################################

"""
# Beispiel aus Folien
A = numpy.array([[20,18,44],[0,40,45],[-15,24,-108]], dtype= float)
B= numpy.array([-4,-45,78], dtype= float)
QR = householderQR(A)

#Berechnen zuerst Y und dann X
y = berechneY(A, QR, B) # y = QTb
X =   householder_rueckwaerts(A, QR, y)
print X
 
#Beispiel aus der Aufgabe 3

print "Aufgabe 3"
#A1 = np.array([[1, -5, -20], [-4, 11, -1], [8, -4, 2]])
A = numpy.array([[20,18,44],[0,40,45],[-15,24,-108]], dtype= float)
QR = householderQR(A) 
print around(QR)
print "-----"
 

"""


print "Aufgabe 5"
n = [40, 50, 60]
for val in n:
    exacteLoesung =  [1 for i in range(val)]
    print "n = " , val 
    A = erstelleMatrixA(val)
    B = erstelleB(val)
    
    HauptDiagonale = householderQR(A) # d 
    Y= berechneY(A, HauptDiagonale, B)
    print "Householder Verfahren"
    print rueckwaerts_HH(A, HauptDiagonale, Y)
 
    
    A = erstelleMatrixA(val)
    LUP = zerlegungPivot(A)
    yP = vorwaerts(LUP[0],permutation(LUP[1],B))
    xP = rueckwaerts(LUP[0],yP)
    
    print "LU-Zerlegung mit Spalten-Pivot"  
    print  around(xP)
    
    print "Exacte Loesung"
    print exacteLoesung
    print "--------------------------------------------"


# Die Ergebnisse beim n = 60 unterscheiden sich. Das kann daran liegen----

# Vereinfachung des Schrittes durch speicherung des Hauptdiagonales speicher der Arbeitaufwand

 
