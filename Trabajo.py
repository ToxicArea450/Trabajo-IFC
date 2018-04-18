#Experimento de Rutherford
# -*- coding: utf-8 -*-
from scipy.integrate import odeint
import numpy as np

#Parámetros del problema
R=7.5e10-14
zs=1.0e-21
qe=1.6e-19
u=1.66e-27 

#Plantear las derivadas
def derivadas(var,t): 
	x, y, vx, vy= var[0], var[1], var[2], var[3] #Los valores iniciales vendrán en el array var
	return [ , , 
	
