#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 12:48:52 2019

@author: Alvaro Gonzalez Bilbao
@function: used to create a channel with two slopes with an ellipsoid as obstacle
@version: 1.02 (22/03/21)
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
import math

class vector:
        def __init__(self,x=0, y=0, z=0):
                self.x= float(x)
                self.y= float(y)
                self.z= float(z)
                
        def List(self):
                return [self.x,self.y,self.z]
                    
        def mag(self):
                return round(math.sqrt(self.x**2+self.y**2+self.z**2),6)
            
        def __add__(self,other):
                return vector(self.x+other.x,self.y+other.y,self.z+other.z)
            
        def __sub__(self,other):
                return vector(self.x-other.x,self.y-other.y,self.z-other.z)
            
        def __mul__(self,alpha):
                return vector(self.x*alpha,self.y*alpha,self.z*alpha)
            
        def __truediv__(self,alpha):
                return vector(self.x/alpha,self.y/alpha,self.z/alpha)

        def __eq__(self, other):
                if self.x == other.x and self.y == other.y and self.z == other.z:
                    return True
                else:
                    return False
          
        def angle(self):              
            if self.x != 0 and self.y != 0:
                if self.x > 0 and self.y > 0:
                    angle = math.atan(self.y/self.x)
                elif self.x < 0:
                    angle = math.pi+math.atan(self.y/self.x)
                else:
                    angle = 2*math.pi+math.atan(self.y/self.x)
            elif self.x == 0:
                if self.y > 0:
                    angle = math.pi/2
                else:
                    angle = 3*math.pi/2
            else:
                if self.x > 0:
                    angle = 0
                else:
                    angle = math.pi                
            return angle
    
cellsize = 2.5
xllcorner = 0
yllcorner = 0
NODATA_value = -9999
filename = 'cuna_two_slopes_20y10.asc'

alpha1 = 10
alpha2 = 20
Lx = 200
Ly = 1000
ym = 250
A1 = 50
A2 = 20
rx = 15
ry = 100
exc_x = 1.05
exc_y = 1.4

ncols = int(Lx/cellsize+1)
nrows = int(Ly/cellsize+1)

x = [float(xllcorner+i*cellsize) for i in range(ncols)]
y = [float(yllcorner+i*cellsize) for i in range(nrows)] 
z = np.zeros((len(x),len(y)))

index = int((ym)/cellsize)
alpha1 = alpha1*np.pi/180
alpha2 = alpha2*np.pi/180

l = 1/2*ym/np.cos(alpha1)
L = ym/2*math.sqrt((1+np.cos(alpha2)/np.cos(alpha1))**2+np.tan(alpha1)**2*(1+np.sin(alpha2)/np.sin(alpha1))**2)
gamma = (alpha2-alpha1)/2
r = L/2*1/math.sin(gamma)
yc = ym/2-r*np.sin(alpha1)
zc = ym/2*np.tan(alpha1)+r*np.cos(alpha1)
yp1 = ym/2
yp3 = ym+l*np.cos(alpha2)

for j in range(nrows):
    if j > int(index):
        if (yp1 <= y[j] or y[j] <= yp3):
            for i in range(ncols):
                z[i,j] = zc-math.sqrt(r**2-(yc-y[j])**2)-A1*np.sin(i*np.pi/(ncols-1))*np.sin((y[j])/(Ly)*np.pi/2)                
        else:
            for i in range(ncols):
                z[i,j] = (y[j]-ym)*np.tan(alpha2)+ym*np.tan(alpha1)-A1*np.sin(i*np.pi/(ncols-1))*np.sin((y[j])/(Ly)*np.pi/2)                
    else:
        if (yp1 <= y[j] or y[j] <= yp3):
            for i in range(ncols):
                z[i,j] = zc-math.sqrt(r**2-(yc-y[j])**2)-A1*np.sin(i*np.pi/(ncols-1))*np.sin((y[j])/(Ly)*np.pi/2)                 
        else:
            for i in range(ncols):
                z[i,j] = (y[j])*np.tan(alpha1)-A1*np.sin(i*np.pi/(ncols-1))*np.sin((y[j])/(Ly)*np.pi/2) 

for j in range(nrows):
    for i in range(ncols):
        if (y[j]-Ly/2*exc_y)**2/(ry)**2+(x[i]-Lx/2*exc_x)**2/(rx)**2<=1:
            z[i,j] += A2*(1-((y[j]-Ly/2*exc_y)**2/(ry)**2+(x[i]-Lx/2*exc_x)**2/(rx)**2))

#figure1 = plt.figure()
#plt.plot(y,z[0,:])
#
#figure2 = plt.figure()
#plt.plot(x,z[:,0])
#
#figure3 = plt.figure()
#plt.plot(x,z[:,nrows-1])
#
#figure4 = plt.figure()
#plt.plot(x,z[:,index])

path = os.getcwd()
path += '/'+filename
fi = open(path, 'w')

fi.write("ncols       {}\n".format(ncols))
fi.write("nrows       {}\n".format(nrows))
fi.write("xllcorner       {}\n".format(xllcorner))
fi.write("yllcorner       {}\n".format(yllcorner))
fi.write("cellsize       {}\n".format(cellsize))
fi.write("NODATA_value       {}\n".format(NODATA_value))

for j in range(nrows):
    for i in range(ncols):
        fi.write(str(z[i,nrows-j-1])+" ")
    fi.write("\n")
fi.close()

#fig = plt.figure(figsize=(9, 3))
#plt.contourf(y,x,z,cmap=plt.get_cmap('hot'))
#plt.colorbar()
#plt.show()

rectangle = np.zeros((5,2))

rectangle[0,0] = 0 ; rectangle[0,1] = 0
rectangle[1,0] = Lx ; rectangle[1,1] = 0
rectangle[2,0] = Lx ; rectangle[2,1] = Ly
rectangle[3,0] = 0 ; rectangle[3,1] = Ly
rectangle[4,0] = rectangle[0,0] ; rectangle[4,1] = rectangle[0,1]

new_rx = 1.5*rx
new_ry = 1.1*ry

deltatheta = 2
ntheta = int(360/deltatheta+1)
theta = [float(0+i*deltatheta*math.pi/180) for i in range(ntheta)]

ellipse = np.zeros((len(theta),2))
for i in range(len(ellipse)):
    ellipse[i,0] = Lx/2*exc_x+rx*math.cos(theta[i])
    ellipse[i,1] = Ly/2*exc_y+ry*math.sin(theta[i])

new_ellipse = np.zeros((len(theta),2))
for i in range(len(new_ellipse)):
    new_ellipse[i,0] = Lx/2*exc_x+new_rx*math.cos(theta[i])
    new_ellipse[i,1] = Ly/2*exc_y+new_ry*math.sin(theta[i])

maxy = Ly/2*exc_y+new_ry*math.sin(math.pi/2)
miny = Ly/2*exc_y+new_ry*math.sin(-math.pi/2)

deltalp = 2
lp1 = np.zeros((int(Ly/deltalp+1),2))
lp2 = np.zeros((int(Ly/deltalp+1),2))

for i in range(len(lp1)):
    lp1[i,1] = Ly-i*deltalp
    lp2[i,1] = lp1[i,1]
    if lp1[i,1] > maxy or lp1[i,1] < miny:
        lp1[i,0] = Lx/2
        lp2[i,0] = Lx/2        
    else:
        theta_i = np.arcsin((lp1[i,1]-Ly/2*exc_y)/new_ry)
        lp1[i,0] = Lx/2*exc_x+new_rx*np.cos(math.pi-theta_i)
        lp2[i,0] = Lx/2*exc_x+new_rx*np.cos(theta_i)
        if lp1[i,0]>Lx/2:
            lp1[i,0] = Lx/2
        if lp2[i,0]<Lx/2:
            lp2[i,0] = Lx/2
            
dist1 = 60
dist2 = 120
dist3 = 80

tp1 = []
tp2 = []

for i in range(len(lp1)-2):
    v1 = vector(lp1[i,0], lp1[i,1], 0)
    v2 = vector(lp1[i+1,0], lp1[i+1,1], 0)
    v3 = vector(lp1[i+2,0], lp1[i+2,1], 0)
    d1 = v1-v2
    d2 = v3-v2
    theta = d2.angle()-d1.angle()
    if theta < 0: theta += 2*math.pi
    elif theta >= 2*math.pi: theta -= 2*math.pi
    theta = theta/2+d2.angle()

    if v2.y > maxy:
        vp1 = vector(math.cos(theta), math.sin(theta),0)*dist1*0.5+v2
        vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist1*0.5+v2 
    elif v2.y < miny:
        vp1 = vector(math.cos(theta), math.sin(theta),0)*dist2*0.5+v2
        vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist2*0.5+v2
    else:
        vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist3*0.5+v2 
        vp1x = Lx/2*exc_x
        vp1y = v2.y-(v2.y-vp2.y)/(v2.x-vp2.x)*(v2.x-vp1x)
        vp1 = vector(vp1x, vp1y,0)
  
    tp1.append([vp1, vp2])

for i in range(len(lp2)-2):
    v1 = vector(lp2[i,0], lp2[i,1], 0)
    v2 = vector(lp2[i+1,0], lp2[i+1,1], 0)
    v3 = vector(lp2[i+2,0], lp2[i+2,1], 0)
    d1 = v1-v2
    d2 = v3-v2
    theta = d2.angle()-d1.angle()
    if theta < 0: theta += 2*math.pi
    elif theta >= 2*math.pi: theta -= 2*math.pi
    theta = theta/2+d2.angle()

    if v2.y > maxy:
        vp1 = vector(math.cos(theta), math.sin(theta),0)*dist1*0.5+v2
        vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist1*0.5+v2 
    elif v2.y < miny:
        vp1 = vector(math.cos(theta), math.sin(theta),0)*dist2*0.5+v2
        vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist2*0.5+v2
    else:
        vp1 = vector(math.cos(theta), math.sin(theta),0)*dist3*0.5+v2
        vp2x = Lx/2*exc_x
        vp2y = v2.y-(v2.y-vp1.y)/(v2.x-vp1.x)*(v2.x-vp2x)
        vp2 = vector(vp2x, vp2y,0)
   
    tp2.append([vp1, vp2])

list_lp1 = [[lp1[0,0],lp1[0,1]]]
list_lp2 = [[lp2[0,0],lp2[0,1]]]
list_tp1 = []
list_tp2 = []

for i in range(len(tp1)):
    if lp1[i+1,1] > maxy+15 or (lp1[i+1,1] < maxy-15 and lp1[i+1,1] > miny+30) or lp1[i+1,1] < miny-20:
        list_lp1.append([lp1[i+1,0],lp1[i+1,1]])
        list_lp2.append([lp2[i+1,0],lp2[i+1,1]])
        list_tp1.append(tp1[i])
        list_tp2.append(tp2[i])        
        
list_lp1.append([lp1[-1,0],lp1[-1,1]])
list_lp2.append([lp2[-1,0],lp2[-1,1]])

lp1 = np.zeros((len(list_lp1),2))
lp2 = np.zeros((len(list_lp2),2))

for i in range(len(lp1)):
    lp1[i,0] = list_lp1[i][0]
    lp1[i,1] = list_lp1[i][1]
    lp2[i,0] = list_lp2[i][0]
    lp2[i,1] = list_lp2[i][1]


path = os.getcwd()
filename = 'lp_wedge_two_slopes_20y10.dat'
path += '/'+filename
fi = open(path, 'w')

fi.write("n_Channels"+" "+str(2)+"\n")
fi.write("\n")
fi.write("Channel1"+"\n")
fi.write("{"+"\n")
for i in range(2):
    for j in range(len(lp1)):
        fi.write(str(round(lp1[j,i],3))+" ")
    fi.write("\n")
fi.write("}"+"\n")
fi.write("Channel2"+"\n")
fi.write("{"+"\n")
for i in range(2):
    for j in range(len(lp2)):
        fi.write(str(round(lp2[j,i],3))+" ")
    fi.write("\n")
fi.write("}"+"\n")
fi.close()

path = os.getcwd()
filename = 'tp_wedge_two_slopes_20y10.dat'
path += '/'+filename
fi = open(path, 'w')

fi.write("n_Channels"+" "+str(2)+"\n")
fi.write("\n")
fi.write("Channel1"+"\n")
fi.write("{"+"\n")
fi.write("Left"+"\n")
fi.write("{"+"\n")
for j in range(len(list_tp1)):
    fi.write(str(round(list_tp1[j][0].x,3))+" ")
fi.write("\n")
for j in range(len(list_tp1)):
    fi.write(str(round(list_tp1[j][0].y,3))+" ")
fi.write("\n")
fi.write("}"+"\n")
fi.write("Right"+"\n")
fi.write("{"+"\n")
for j in range(len(list_tp1)):
    fi.write(str(round(list_tp1[j][1].x,3))+" ")
fi.write("\n")
for j in range(len(list_tp1)):
    fi.write(str(round(list_tp1[j][1].y,3))+" ")
fi.write("\n")
fi.write("}"+"\n")
fi.write("}"+"\n")
fi.write("Channel2"+"\n")
fi.write("{"+"\n")
fi.write("Left"+"\n")
fi.write("{"+"\n")
for j in range(len(list_tp2)):
    fi.write(str(round(list_tp2[j][0].x,3))+" ")
fi.write("\n")
for j in range(len(list_tp2)):
    fi.write(str(round(list_tp2[j][0].y,3))+" ")
fi.write("\n")
fi.write("}"+"\n")
fi.write("Right"+"\n")
fi.write("{"+"\n")
for j in range(len(list_tp2)):
    fi.write(str(round(list_tp2[j][1].x,3))+" ")
fi.write("\n")
for j in range(len(list_tp2)):
    fi.write(str(round(list_tp2[j][1].y,3))+" ")
fi.write("\n")
fi.write("}"+"\n")
fi.write("}"+"\n")
fi.close()

#fig = plt.figure(figsize=(8, 20))
#plt.plot(rectangle[:,0], rectangle[:,1])
#plt.plot(ellipse[:,0], ellipse[:,1])
#plt.plot(new_ellipse[:,0], new_ellipse[:,1])
#plt.plot(lp1[:,0], lp1[:,1])
#plt.plot(lp2[:,0], lp2[:,1])
#for i in range(len(list_tp1)):
#    vertexes = list_tp1[i]
#    plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])
#for i in range(len(list_tp2)):
#    vertexes = list_tp2[i]
#    plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])
#plt.grid()
#plt.show()