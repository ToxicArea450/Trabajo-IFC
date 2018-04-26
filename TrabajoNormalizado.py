# -*- coding: utf-8 -*-
#Experimento de Rutherford
#En este caso normalizaremos las ecuaciones para no trabajar con números muy grandes o pequeños.
from scipy.integrate import odeint
from scipy.optimize import curve_fit
from numpy import *
import matplotlib.pyplot as plt
from math import pi, atan

#Parámetros del problema
R, zs, qe, u, e0= 7.5e-14, 1.0e-21, 1.6e-19, 1.66e-27, 8.08542e-12

#Datos de las cargas y masas
Q=79        #Carga del núclero de oro
q=2          #Carga de la partícula alfa (núcleo de Helio)
m=4          #Masa de la partícula alfa 

cte=((qe**2)*(zs**2)*q*Q)/(u*(R**3)*m*4*pi*e0) #Constante que multiplica a las ecuaciones del movimiento normalizadas

#Definimos la función derivadas que queremos integrar
def derivadas(var,t):
	'''
	Esta función devuelve las derivadas de x, y, vx y vy
	Es la que usará odeint para poder integrar estas derivadas y obtener las trayectorias de la partícula alfa
	'''
	x, y, vx, vy= var[0], var[1], var[2], var[3]
	ax=cte*(x/((x**2+y**2)**(3./2)))
	ay=cte*(y/((x**2+y**2)**(3./2)))
	return array([vx, vy, ax, ay])

#Definimos la función recta para poder calcular el ángulo de dispersión
def recta(x, a, b):
	'''
	Esta función expresa la funcion y=ax+b, es decir, es la función de una recta.
	La usaremos para hallar los ángulos de dispersión de las particulas a partir de la pendiente de su recta de ajuste
	'''
	return x*a+b
	
#Creamos el vector de tiempos e introducimos las condiciones iniciales
tiempos=linspace(0,100,1000)
x, vx, vy= -100, 2, 0.0
b=random.normal(0,0.1,(50,)) #Genera un número aleatorio de números que siguen una distribución normal N(0,0.08), cada número será el parámetro de impacto de cada partícula que vamos a lanzar
n=b.size

#Como b va a tener varios valores, ya que tenemos mas de una partícula, hacemos un bucle for en el que para cada valor se llame a la función odeint
phit=zeros_like(b) #Array de 0 con el mismo tamaño que el array b al que le iremos metiendo los datos de cada ángulo para luego poder representarlos

for i in range(n-1):
	y=b[i] #Parámetro de impacto de cada particula
	trayectorias=odeint(derivadas,array([x, y, vx, vy]), tiempos)
	xt, yt, vxt, vyt = trayectorias[:,0], trayectorias[:,1], trayectorias[:,2], trayectorias[:,3]
	
	#Dibujo de las trayectorias
	plt.subplot(2,1,1)
	plt.plot(xt, yt)
	
	#Cálculo del ángulo de dispersión
	p, pcov=curve_fit(recta, xt[-50:], yt[-50:])
	pendiente=p[0]
	angulo=atan(pendiente)*(180/pi) #El problema que tiene este ángulo es que si la particula se desvia mas de 90º te da el angulo complemetario cambiado de signo, no el que nos interesa
	#Para solucionar ese problema hacemos lo siguiente:
	if y>0:                   #Si el parámetro de impacto es positivo, la pendiente será positiva por lo que el ángulo también, si el ángulo sale negativo significa que nos han dado el complementario.
		if pendiente >=0:
			phi=angulo
		elif pendiente <0:
			phi=(180-angulo)
	if y<0:                   #Si el parámetro de impacto es negativo, la pendiente será negativa por lo que el ángulo también, si el ángulo sale positivo significa que nos han dado el complementario.
		if pendiente <=0:
			phi=angulo
		elif pendiente >0:
			phi=(180-abs(angulo))*(-1)	
	phit[i]=phi
	
plt.plot(0,0,'ro')
plt.title('Trayectoria de la particula al chocar con otra')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(2,1,2)
plt.title('Relacion parametro de impacto - angulo de dispersion')
plt.plot(b,phit,'bo')
plt.xlabel('b')
plt.ylabel('phit')
plt.show()
