 
from pylab import *

from Aufgabe05lib import *
#Konjugierte Gradienten
def CG(A,b,x0):
    xo=array(x0);
    # p0 = r0 = b - Ax0
    r0=b-A.dot(x0);
    p0=array(r0)        
    rAnfang=array(r0)
    arr = []
    
    while True:
        # <rk, rk>
        rknorm = r0.dot(r0)
        # alfa =  <rk, rk> / < pk , Apk> 
        alfa = rknorm/p0.dot(A.dot(p0))
        # xk+1 = xk + alfa*pk
        x1=x0+alfa*p0
        # rk+1 == r1
        # rk+1 = rk - alfa A * pk
        r1=r0-alfa*A.dot(p0)
        # norm beta
        # |beta| = < rk+1 , rk+1> / <rk , rk>
        beta=r1.dot(r1)/rknorm
        # pk+1 = rk+1 + bk+1
        p1=r1+beta*p0
        
        div = norm(r1)/norm(rAnfang)# Hier vergleichen wir aktuellen Wert mit Anfangswert
        arr.append(div)
        if(div<10**(-6)): # Falls Naeherung eps. wert klein ist, dann malen wir schon Grafik 
            break
        else: # sonst wiederhollen wir die Berechnung
            x0=x1;
            r0=r1;
            p0=p1;
            
    return x1, arr
# Aufgabe  5 Implementierung
for i in  [50,100,200]:
    x0 = zeros(i*i)  
    Ab = system(i)
    res = CG(Ab[0],Ab[1],x0)
    figure()
    semilogy(res[1])
    plotxk(res[0])

 