from pylab  import *


def LUPivot( A):
    A = array(A)
    print A
    dim = A.shape[0]
    p = zeros(dim-1, int)
    # erstellen L und setze 1 auf Diagonal
 
    #erstelle U
    U = [[0 for i in range(0,dim)] for i in range(dim)]
    for i in range(0,dim):
        for j in range(0,dim):
            U[i][j] = A[i][j]

    for i in range( 0, dim-1): # Zeilenweiße

        # Suchen max element
        maxElem = abs(U[i][i])
        maxZeile = i

        #laufen jede Zeile
        for k in range( i+1, dim):
            if( abs(U[k][i]) > maxElem):
                maxElem = abs(U[k][i])
                maxZeile = k
        # Die größte Zeil gefunden
        p[i] = maxZeile

        #Tauschen die Zeile

        for k in range(i, dim) :
            tmp=U[maxZeile][k]
            U[maxZeile][k]=U[i][k]
            U[i][k]=tmp
 
 
        
        for k in range(i+1,dim):#erste  Zeil bleibt unveränderlich
            koeff = U[k][i]/float(U[i][i])
        
           
            c = -U[k][i]/float(U[i][i])
        
            for j in range(i, dim):
                U[k][j] += c*U[i][j] # Multiply with the pivot line and subtract


    
  
    print U
    print p
        
    




print "LU Zerlegung"          
A = array([[0, 0, 0, 1], [2, 1, 2, 0], [4, 4, 0, 0], [2, 3, 1, 0]])
print LUPivot(A)
