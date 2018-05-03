#Experimento de Rutherford (Modo Euler)
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
def derivadas(x, y, vx, vy, t):
	return cte*(x/((x**2+y**2)**(3./2))), cte*(y/((x**2+y**2)**(3./2)))
	
#Definimos el método de Euler, que es el que usaremos para integrar la función
def euler(funcion, x, y, vx, vy, tiempos):
	n=tiempos.size
	xt, yt, vxt, vyt= zeros_like(tiempos), zeros_like(tiempos), zeros_like(tiempos), zeros_like(tiempos)
	xt[0], yt[0], vxt[0], vyt[0]= x, y, vx, vy
	print '\n    x    |    y     |    vx    |    vy    '
	print '-----------------------------------------'
	for k in range(1,n):
		t=tiempos[k-1]
		Dt=(tiempos[k]-tiempos[k-1])
		ax, ay= funcion(x, y, vx, vy, t)
		t, vx, vy, x, y= t+Dt, vx+ax*Dt, vy+ay*Dt, x+vx*Dt, y+vy*Dt
		xt[k], yt[k], vxt[k], vyt[k]= x, y, vx, vy
		print '%.2e | %.2e | %.2e | %.2e '%(x,y,vx,vy)
	return xt, yt, vxt, vyt
	
#Creamos el vector de tiempos e introducimos las condiciones iniciales
tiempos=linspace(0,100*zs,100)
x, vx, vy= 100*R, -(2*R)/zs, 0.0
b=arange(0,1,0.01) #Parámetro de impacto, es la distancia entre la dirección de la particula incidente y el centro de fuerzas, es decir, la y.

#
for k in b:
	y=k*R
	xt, yt, vxt, vyt=euler(derivadas, x, y, vx, vy, tiempos)
	plt.plot(xt, yt)
	
plt.plot(0,0,'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
