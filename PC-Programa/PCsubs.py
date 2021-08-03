#Subrrutinas para calculo de cristales fotónicos (PC) en 1-D 
#Por: Fidel Alejandro Navarro Salazar
#-----------------------------------------------------------------------------
#Este archivo contiene las funciones necesarias para la definición de matrices y cálculo de los PC

#Usamos Numpy para generar arreglos de los datos
import numpy as np
#Usamos Scipy para realizar un ajuste de los datos experimentales
from scipy.interpolate import interp1d as fit
#usamos Random para generar números pseudo-aleatorios
import random

def Beta(n,l, theta): 
    return n * 2 * np.pi *np.cos(theta) / l

#matriz transmicion
def Tr(n1,n2):
    T_n = np.zeros( ( 2 , 2 ) , dtype=np.complex )
    
    T_n[ 0 , 0 ] = n1 + n2
    T_n[ 0 , 1 ] = n1 - n2
    T_n[ 1 , 0 ] = n1 - n2
    T_n[ 1 , 1 ] = n1 + n2
    
    return T_n / (n1 * 2)
    
def P(Beta1 ,a) :
    #matriz propagacion
    P = np.zeros( ( 2 , 2 ) , dtype=np.complex ) #propagación material

    P[ 0 , 0 ] = np.exp( -1j * Beta1 * a )
    P[ 0 , 1 ] = 0
    P[ 1 , 0 ] = 0
    P[ 1 , 1 ] = np.exp(1j * Beta1 * a )
    
    return P

#calcula el producto de las matrices
def Trans (n0, n1, n2, Beta1 , Beta2, wavelength , m, d) :
    #Los valores igualados a cero corresponden al defecto
    #si no se declaran el programa los tomará como 0 y se calcularán las matrices de forma regular
    
    To = [] #se crea una lista, esta lista contiene las matrices de transmisión de manera ordenada
    for i in range(len(d)+1):   #para N capas se tienen N+1 matrices de transmición
        if i==0:                #definimos los valores de las mattrices: cuando i==0 estaremos en la primera interfase (aire-material1)
            To.append(Tr(n0,n1))
        elif i==len(d) :        #cuando i==N estaremos en la última interfase (material2-aire)
            To.append(Tr(n2,n0))
        elif i%2==0. :          #cuando i sea número par estaremos en la interfase material2-material1
            To.append(Tr(n2,n1))
        else:                   #para todos los demás casos (es decir los impares) estaremos en la interfase material1-material2
            To.append(Tr(n1,n2))
    
    #Similarmente definimos las matrices de propagación
    Po = [] 
    for i in range(len(d)):             #para N capas tendremos N matrices de propagación
        if i%2==0. :                    #para i pares P correspondera al del material 1
            Po.append( P(Beta1, d[i]) ) 
        else:                           #para cualquier otra i (es decir impares) P correspondera al material 2
            Po.append( P(Beta2, d[i]) ) 
    
    M = np.dot(To[0], Po[0])    #podemos comenzar multiplicando las primeras 2 matrices To[0]=T_aire_mat1 y Po=P_mat1
    
    #Para multiplicar las matrices faltantes podemos hacer uso del siguiente contador
    for i in range(len(Po)-1):                  #Como ya multiplicamos una matriz Po y una To nos quedan N-1 matrices Po y N matrices To
        M = np.dot(M, np.dot(To[i+1],Po[i+1]) ) #Multiplicamos el valor previo de M por el producto de To[i] y Po[i] : M x (To x Po)
        
    M = np.dot(M, To[-1]) #previo a esta linea se utilizaron las primeras N matrices To, solo queda la última (material 2 - aire)
    

    return M #regresamos la matriz de tranferencia
#ley de snell
def snell(n1,n2,theta):
    theta2 = np.arcsin( ( n1/n2 )  * np.sin(theta) )
    return theta2

#Función que organiza la información y llama a las funciones anteriores
#def cristal(n0, n1, n2, theta, l, a, b, N, nD=0, c=0, M=0):  #versión para valores constantes de n y k
def cristal(n0, n_1, k_1, n_2, k_2, theta, l, d, M): #vérsión para valores experimentales de n y k
    r = [ ]
    t = [ ]
    for i in l:  
        n1 = n_1(i) + 1j* k_1(i)  #comentar en el caso de valores constantes
        n2 = n_2(i) + 1j* k_2(i)  #comentar en el caso de valores constantes
        
        theta1 = snell(n0, n1.real, theta) #ángulo material 1 
        theta2 = snell(n0, n2.real, theta) #ángulo material 2


        Beta1 = Beta(n1,i,theta1) #material 1
        Beta2 = Beta(n2,i,theta2) #material 2
            
        #Trans (n0, n1 , n2, Beta1 , Beta2, wavelength , a, b, n, nD=0, BetaD=0, cD=0,m=0 )
        Final_M = Trans(n0, n1, n2, Beta1 , Beta2, i , M, d)
        #calculo de transmitancia
        to = 1/Final_M[0,0]
        T = (((to)*to.conjugate())).real
        #calculo de reflectancia
        ro = Final_M[1,0] / Final_M[0,0]
        R = (((ro)*ro.conjugate())).real
        t.append(T)
        r.append(R)
    
    return r

def PC(M,mat1,mat2, l, espesor1, espesor2, d=0):
    # M : número de capaz
    # mat1 : dataframe con valores (longitud de onda, n, k) del material 1
    # mat2 : dataframe con valores (longitud de onda, n, k) del material 2
    # d : lista de los espesores de las capaz, si no se especifica el valor se tomara como d=0 y se asignaran valores aleatorios
    #dd = 1 #ancho de paso

    #Obtenemos la longitud de onda incial, para ello se compara la longitud de onda mínima de los tres materiales
    #lini = [mat1[0].min(), mat2[0].min()]
    #lini = max(lini)
    #Obtenemos la longitud de onda final, para ello se compara la longitud de onda máxima de los tres materiales
    #lfin = [mat1[0].max(), mat2[0].max()]
    #lfin = min(lfin)

    #Arreglo con longitudes de onda
    #l = np.linspace(lini,lfin, int((lfin-lini)/dd) +1)

    #indices de refración del aire
    n0 = 1. #aire


    #espesor material 
    if d==0:
        d=[]
        for i in range(M):  #si no se asignaron valores a d entonces se calculará una lista de longitud M con valores entre 10-1000 nm
            d.append(random.uniform(espesor1, espesor2))


    #angulo de incidencia
    theta = np.pi*0.

    #Ajustes
    #material 1
    n_1=fit(mat1[0],mat1[1],kind='cubic') #ajuste de indice de refraccion 
    k_1=fit(mat1[0],mat1[2],kind='cubic') #ajuste de índice de extincion 
    #material 2
    n_2=fit(mat2[0],mat2[1],kind='cubic') #ajuste de indice de refraccion 
    k_2=fit(mat2[0],mat2[2],kind='cubic') #ajuste de índice de extincion

    r = cristal(n0, n_1, k_1, n_2, k_2, theta, l, d, M)
    
    return(r,d)
