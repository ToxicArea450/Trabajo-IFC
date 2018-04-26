#Experimento de Rutherford
# -*- coding: utf-8 -*-
from scipy.integrate import odeint
from numpy import *
import matplotlib.pyplot as plt
from math import pi, atan

#Parámetros del problema
R, zs, qe, u, e0= 7.5e-14, 1.0e-21, 1.6e-19, 1.66e-27, 8.08542e-12

#Datos de las cargas y masas
Q=79*qe         #Carga del núclero de oro
q=2*qe          #Carga de la partícula alfa (núcleo de Helio)
m=4*u          #Masa de la partícula alfa 

cte=(q*Q)/(4*pi*e0*m) #Constante que multiplica a las ecuaciones de las trayectorias

#Definimos la función derivadas que queremos integrar
def derivadas(var,t):
	x, y, vx, vy= var[0], var[1], var[2], var[3]
	ax=cte*(x/((x**2+y**2)**(3./2)))
	ay=cte*(y/((x**2+y**2)**(3./2)))
	return array([vx, vy, ax, ay])

#Definimos la función recta para poder calcular el ángulo de dispersión
def recta(x, a, b):
	return x*a+b
	
#Creamos el vector de tiempos e introducimos las condiciones iniciales
tiempos=linspace(0,100*zs,1000)
x, vx, vy= 100*R, -(2*R)/zs, 0.0
b=random.normal(0,0.5,(50,)) #Parámetro de impacto, es la distancia entre la dirección de la particula incidente y el centro de fuerzas, es decir, la y.
n=b.size

#Como b van a ser varios números, acemos un bucle for en el que para cada número se llame a la función odeint

phit=zeros_like(b) #Array de 0 con el mismo tamaño que el array b al que le iremos metiendo los datos de cada ángulo

for i in range(n-1):
	y=b[i]*R*5
	trayectorias=odeint(derivadas,array([x, y, vx, vy]), tiempos)
	xt, yt, vxt, vyt = trayectorias[:,0], trayectorias[:,1], trayectorias[:,2], trayectorias[:,3]
	
	#Dibujo de las trayectorias
	plt.subplot(2,1,1)
	plt.plot(xt, yt)
	
	#Ángulo de dispersión
	k=1
	E=4
	phi=2*atan(k/(2*E*b[i]))*180/pi
	phit[i]=phi
	
plt.plot(0,0,'ro')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(2,1,2)
plt.plot(b,phit,'bo')
plt.xlabel('b')
plt.ylabel('phit')
plt.show()

