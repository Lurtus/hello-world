#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:39:15 2019

@author: eric
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate
from scipy import optimize
import pdb
import Nbodyodes as nbo



#t0 = 0
#t1 = 100
#
#m1 = 1
#m2 = 1
#
#G = 5
#
#x1_0 = 0
#y1_0 =  0
#vx1_0 = 0.1
#vy1_0 = -0.1
#
#x2_0 = 5
#y2_0 =  0
#vx2_0 = 0
#vy2_0 = 1

#t0 = 0
#t1 = 100
#
#m1 = 0.05
#m2 = 1
#
#G = 10
#
#x1_0 = 0
#y1_0 =  0
#vx1_0 = 0
#vy1_0 = 0.1
#
#x2_0 = 5
#y2_0 =  1
#vx2_0 = -1
#vy2_0 = 1.2



#Earh moon system
#t0 = 0
#t1 = 24*3600*30*2
#G=6.67408e-11
#m1 =  5.972e24 
#m2 = 7.342e22*16
#x1_0 = 0
#y1_0 =  0
#vx1_0 = 0
#vy1_0 = 0
#
#x2_0 = 384400000
#y2_0 =  0
#vx2_0 = 0
#vy2_0 = 1022


###Earh Sun system
t0 = 0
t1 = 24*3600*365
G=6.67408e-11
m1 =  1.989e30 
m2 = 5.972e24 
x1_0 = 0
y1_0 =  0
vx1_0 = 0
vy1_0 = 0

x2_0 = 149.6e9
y2_0 =  0
vx2_0 = 0
vy2_0 = 107000/3.6


###Earh-Satalite system
#t0 = 0
#t1 = 3600*24*30
#G=6.67408e-11
#m1 =  5.972e24 
#m2 = 2e3
#x1_0 = 0
#y1_0 =  0
#vx1_0 = 0
#vy1_0 = 0
#
#x2_0 = 1e8
#y2_0 =  8e7
#vx2_0 = 1.5e3
#vy2_0 = -1e3

#107,000 km/h

Y0 = [x1_0, y1_0, vx1_0, vy1_0, x2_0, y2_0 , vx2_0, vy2_0]

nstep = 150
timeinterval = np.linspace(t0,t1,nstep)
dt_save = (t1-t0)/(nstep+1)


rtolode = 1e-6 #default values
atolode = 1e-6 #default values

#fun_wrap = lambda t, Y: nbo.odefun_2body_2d(t,Y,m1,m2,G) 
#sol = integrate.solve_ivp(fun_wrap,[t0,t1],Y0, method='RK45',rtol=rtolode, atol=atolode, t_eval = timeinterval)

dt = t1/1e7
fun_wrap = lambda Y: nbo.acc_2body_2d(Y,m1,m2,G) 
sol = nbo.leapfrogsolver_2d_2body(fun_wrap,[t0,t1],Y0, dt, dt_save)


time = sol.t
x1_all = sol.y[0,:]
y1_all = sol.y[1,:]

x1dot_all = sol.y[2,:]
y1dot_all = sol.y[3,:]

x2_all = sol.y[4,:]
y2_all = sol.y[5,:]

x2dot_all = sol.y[6,:]
y2dot_all = sol.y[7,:]

plt.figure(num=1)
   
for t,x1,y1,x2,y2,x1dot,y1dot,x2dot,y2dot in zip(time,x1_all,y1_all,x2_all,y2_all,x1dot_all,y1dot_all,x2dot_all,y2dot_all):
    plt.figure(1)
    plt.plot(x1,y1,'bo')
    plt.plot(x2,y2,'ro')
    plt.show()   
    plt.pause(0.2)


    #distance
    r12 = ((x2-x1)**2 + (y2-y1)**2)**(1/2)
    v1 = np.sqrt(x1dot**2 + y1dot**2)
    v2 = np.sqrt(x2dot**2 + y2dot**2)

    #Potential Energy
    U = -G*m1*m2/r12 
    #U = 0
    
    #Kinetic Energy
    T = (m1*v1**2)/2 + (m2*v2**2)/2
    
    E = U + T
    print(E)
    plt.figure(2)
    #plt.plot(t,U,'bo')
    #plt.plot(t,T,'ro')
    plt.plot(t,E,'go')
    plt.show()   
    plt.pause(0.2)
#jac_wrap = lambda z, Y: de.odejac(z,Y,n,thick,Epar,km,beta,dmc,Deltamc,es0)