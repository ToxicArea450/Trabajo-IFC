#Experimento de Rutherford
# -*- coding: utf-8 -*-
from scipy.integrate import odeint
from numpy import *
import matplotlib.pyplot as plt
from math import pi

#Parámetros del problema
R=7.5e10-14
zs=1.0e-21
qe=1.6e-19
u=1.66e-27 
e0=808542e-12

#Datos de las cargas y masas
Q=79         #Carga del núclero de oro
q=2          #Carga de la partícula alfa (núclero de Helio)
M=4          #Masa de la partícula alfa 

cte=(qe**2*zs**2*q*Q)/(u*R**3*M*4*pi*e0)

#Plantear las derivadas
def derivadas(var,t): 
	x, y, vx, vy= var[0], var[1], var[2], var[3] #Los valores iniciales vendrán en el array var
	ax=cte*(x/(x**2+y**2)**(3/2))
	ay=cte*(y/(x**2+y**2)**(3/2))
	return array([vx, vy, ax, ay])
	
tiempos=linspace(0,100,1000)

#Condiciones iniciales
x, y, vx, vy= -100, 10., 2, 0

trayectorias= odeint(derivadas, array([x, y, vy, vy]), tiempos)
xt, yt, vxt, vyt = trayectorias[:,0], trayectorias[:,1], trayectorias[:,2], trayectorias[:,3]

#Dibujo de trayectorias
plt.plot(tiempos, xt, tiempos, yt, tiempos, vxt, tiempos, vyt)
plt.legend(["Posicion en x", "Posicion en y", "Velocidad en x", "Velocidad en y"], loc="best")
plt.show()
