#Experimento de Rutherford
# -*- coding: utf-8 -*-
from scipy.integrate import odeint
from numpy import *
import matplotlib.pyplot as plt
from math import pi, atan

#Parámetros del problema
R=7.5e10-14
zs=1.0e-21
qe=1.6e-19
u=1.66e-27 
e0=8.08542e-12

#Datos de las cargas y masas
Q=79*qe         #Carga del núclero de oro
q=2*qe          #Carga de la partícula alfa (núcleo de Helio)
M=4*u          #Masa de la partícula alfa 

cte=(q*Q)/(M*4*pi*e0)
#Plantear las derivadas
def derivadas(var,t): 
	x, y, vx, vy= var[0], var[1], var[2], var[3] #Los valores iniciales vendrán en el array var
	ax=cte*(x/(x**2+y**2)**(3./2))
	ay=cte*(y/(x**2+y**2)**(3./2))
	return array([vx, vy, ax, ay])
	
tiempos=linspace(0,100*zs,1000)

#Condiciones iniciales
x, y, vx, vy= 10*R, R, -2.e-7, 0.
b=random.randint(0,10,(5,))

#Bucle para meter todos los b
'''
for y in b:
	trayectorias= odeint(derivadas, array([x, y, vx, vy]), tiempos)
	xt, yt, vxt, vyt = trayectorias[:,0], trayectorias[:,1], trayectorias[:,2], trayectorias[:,3]
	print xt
	print yt
	
	#Dibujo de trayectorias
	plt.plot(xt, yt,'b-')
	
plt.plot(0,0,'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.axis()
plt.show()
'''
trayectorias= odeint(derivadas, array([x, y, vx, vy]), tiempos)
xt, yt, vxt, vyt = trayectorias[:,0], trayectorias[:,1], trayectorias[:,2], trayectorias[:,3]
print xt
print "----------------------------"
print yt

plt.plot(xt, yt,'b-')	
plt.plot(0,0,'ro')
plt.xlabel('x')
plt.ylabel('y')
plt.axis()
plt.show()

#Ángulo de dispersión
k=1./4*pi*e0
E=0.5*M*vx
phi=atan(float(k)/(-82*E*b))*180/pi
