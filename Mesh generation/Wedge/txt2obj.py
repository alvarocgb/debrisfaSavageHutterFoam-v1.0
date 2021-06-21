#!/usr/bin/env python 

'''
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Code used to create an OBJ-file from a Raster-file
Author
    Matthias Rauter matthias.rauter@uibk.ac.at
Modified by
    Álvaro González Bilbao alvaro.gonzalez.bilbao@gmail.com
    
@version: 1.04 (21/06/21)
'''

# these 3 lines make python2 compatible with python3
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from scipy import interpolate
import re
import argparse
import sys
from time import gmtime, strftime
from operator import itemgetter
import math

##########################################  Definitions  ##############################################################

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
            if type(self) == type(alpha):
                return self.x*alpha.x+self.y*alpha.y+self.z*alpha.z   
            else:
                return vector(self.x*alpha,self.y*alpha,self.z*alpha)
        
    def __truediv__(self,alpha):
            return vector(self.x/alpha,self.y/alpha,self.z/alpha)

    def __eq__(self, other):
            if self.x == other.x and self.y == other.y and self.z == other.z:
                return True
            else:
                return False
          
def distxy(v1,v2):
    return round(math.sqrt((v1.x-v2.x)**2+(v1.y-v2.y)**2),6)

def N(xis):
    return 1/4.*(1+xis[0])*(1+xis[1])

def getX(xi, x0):
    xi0 = np.array([1, -1])
    xi1 = np.array([1, 1])
    xi2 = np.array([-1, 1])
    xi3 = np.array([-1, -1])
    
    xg = x0[0]*N(xi*xi0)+x0[1]*N(xi*xi1)+x0[2]*N(xi*xi2)+x0[3]*N(xi*xi3)
    return xg

pattern = r"""      
    \-*[0-9]+\.*[0-9]*
    """
#re module offers a set of functions that allow us to search a string for a match.
rx = re.compile(pattern, re.VERBOSE) 

def coords(s):
    try:        
        s = s[1:len(s)-1]
        x, y = map(float, s.split(','))
        return x, y
    except:
        raise argparse.ArgumentTypeError("Coordinates must be 'x,y'")

def get_alpha_list(alpha,name):
    if alpha is not None:
        alpha = list(alpha.split(","))
        new_list = [float(i) for i in alpha]
        alpha = new_list
        alpha.sort()
        if int(len(alpha)/2) != float(len(alpha)/2):
            raise Exception('{} should have an even number of values. The value of {} was: {}'.format(name, name,alpha))
        else:
            return alpha

def get_zone(zone, sense, offsetx, offsety):
    if zone is not None:
        zone = zone.split(";")
        for i in range(len(zone)):
            zone[i] = zone[i].replace('(','')
            zone[i] = zone[i].replace(')','')
            zone[i] = zone[i].split(",")
            zone[i][0] = float(zone[i][0])+offsetx
            zone[i][1] = float(zone[i][1])+offsety
        if sense == 'x':
            zone = sorted(zone,key=itemgetter(0))
        else:
            zone = sorted(zone,key=itemgetter(1))
        return zone
    else:
        return []

def check_zone(zone, minx, maxx, miny, maxy):
    if zone != []:
        for i in range(len(zone)):
            if zone[i][0] < minx:
                zone[i][0] = minx
            elif zone[i][0] > maxx:
                zone[i][0] = maxx
            if zone[i][1] < miny:
                zone[i][1] = miny
            elif zone[i][1] > maxy:
                zone[i][1] = maxy
            zone[i] = vector(zone[i][0], zone[i][1], 0)

def get_boundary_points(list, xyloc):
    new_zone = []
    if len(list) == 2:
        for i in range(len(xyloc)):
            new_zone.append(list[0]+(list[1]-list[0])*xyloc[i])
    else:
        dist = 0
        for i in range(len(list)-1):
            dist += distxy(list[i], list[i+1])
            
        index_list = [0]
        new_dist = 0
        for i in range(len(list)-1):
            new_dist += distxy(list[i], list[i+1])
            index_list.append(new_dist/dist)

        last_index = 0
        last_xloc = 0
        for i in range(len(xyloc)-1):
            coeff = 1/(index_list[last_index+1]-index_list[last_index])
            new_zone.append(list[last_index]+(list[last_index+1]-list[last_index])*(xyloc[i]-last_xloc)*coeff)
            if xyloc[i] < index_list[last_index+1] and xyloc[i+1] >= index_list[last_index+1]:
                last_xloc = xyloc[i+1]
                last_index += 1
                
        new_zone.append(list[-1])
        
    return new_zone

##########################################  Main function  #############################################################

def main(argv):
    
    parser = argparse.ArgumentParser(description='Taking a Raster-File and creating a OBJ-File')

    parser.add_argument('-i', type=str, help='input filename')
    parser.add_argument('-o', type=str, help='output filename')
    parser.add_argument('-xres', type=int, help='number of points in the x-axis')
    parser.add_argument('-yres', type=int, help='number of points in the y-axis')
    parser.add_argument('-p1', type=str, help='lower left edgepoint')
    parser.add_argument('-p2', type=str, help='upper left edgepoint')
    parser.add_argument('-p3', type=str, help='upper right edgepoint')
    parser.add_argument('-p4', type=str, help='lower right edgepoint')
    parser.add_argument('-fillup', help='fill up NaN values with nearest neighbor', action="store_true")
    parser.add_argument('-exactcopy', help='create a triangulation with the resolution of the dem', action="store_true")
    parser.add_argument('-offsetx', type=float, help='offset mesh by x m', default=0)
    parser.add_argument('-offsety', type=float, help='offset mesh by y m', default=0)
    parser.add_argument('-alphamaxY', type=str, help='values to divide the maxY patch for different boundary conditions, it must be an even number of values')
    parser.add_argument('-alphaminY', type=str, help='values to divide the minY patch for different boundary conditions, it must be an even number of values')
    parser.add_argument('-alphamaxX', type=str, help='values to divide the maxX patch for different boundary conditions, it must be an even number of values')
    parser.add_argument('-alphaminX', type=str, help='values to divide the minX patch for different boundary conditions, it must be an even number of values')    
    parser.add_argument('-maxY', type=str, help='points coordinates of the boundary patch')
    parser.add_argument('-minY', type=str, help='points coordinates of the boundary patch')
    parser.add_argument('-maxX', type=str, help='points coordinates of the boundary patch')
    parser.add_argument('-minX', type=str, help='points coordinates of the boundary patch') 
    parser.add_argument('-zmax', type=float, help='z maximum value used for the creation of the finite volume mesh', default=0)
    args = parser.parse_args()
        
    if args.xres is not None:
        xres = args.xres
    else:
        xres = 100
    if args.yres is not None:
        yres = args.yres
    else:
        yres = 100

    infile  = args.i
    outfile = args.o   
    p1in = args.p1
    p2in = args.p2 
    p3in = args.p3 
    p4in = args.p4 
    fillup = args.fillup
    exactcopy = args.exactcopy
    zmax = args.zmax

    #This is used to separate a patch in subpatches where different boundary conditions may be applied
    alphamaxY = get_alpha_list(args.alphamaxY, 'alphamaxY')
    alphaminY = get_alpha_list(args.alphaminY, 'alphaminY')
    alphamaxX = get_alpha_list(args.alphamaxX, 'alphamaxX')
    alphaminX = get_alpha_list(args.alphaminX, 'alphaminX')

    #This is used when the user wants to create a mesh with a shape different to the one of the original raster-file
    maxY = get_zone(args.maxY,'x',args.offsetx, args.offsety)
    minY = get_zone(args.minY,'x',args.offsetx, args.offsety)
    maxX = get_zone(args.maxX,'y',args.offsetx, args.offsety)
    minX = get_zone(args.minX,'y',args.offsetx, args.offsety)
    
    if p1in is not None:
        p1in = coords(args.p1)
        p1in = np.array(p1in)
    if p2in is not None:
        p2in = coords(args.p2)
        p2in = np.array(p2in)
    if p3in is not None:
        p3in = coords(args.p3)
        p3in = np.array(p3in)
    if p4in is not None:
        p4in = coords(args.p4)
        p4in = np.array(p4in)

    print("Reading from {}".format(infile))
    print("Writing to {}".format(outfile))    
    if exactcopy: print("creating an exact copy of the raster")
    if fillup:    print("filling up nodata values")
    
    print("moving by ({}, {})".format(args.offsetx, args.offsety))
    
    file  = open(infile, 'r')
    ncols = int(rx.findall(file.readline())[0]) 
    nrows = int(rx.findall(file.readline())[0])
    xllcenter = float(rx.findall(file.readline())[0])
    yllcenter = float(rx.findall(file.readline())[0])
    cs = float(rx.findall(file.readline())[0])

    xllcenter = xllcenter + args.offsetx 
    yllcenter = yllcenter + args.offsety
    
    x = [float(xllcenter+i*cs) for i in range(ncols)]
    y = [float(yllcenter+(nrows-i-1)*cs) for i in range(nrows)] 

    i = 0
    z = np.zeros((nrows, ncols))
    for line in file:
        st = rx.findall(line)
        if len(st) == ncols:
            z[i,:] = [float(s) for s in st]
            i = i+1

    if exactcopy:
        xres = ncols
        yres = nrows

    #It may be the case that some of the z-values in the raster are wrong, in which case they are replaced for valid values.
    if fillup:
        fn = True 
        k = 1
        while fn and k < 10:
            print("filling up ({}) ...".format(k)); k = k+1
            fn = False #if in an iteration there is no more filling up, then fn = False and the iteration stops
            fixed = 0
            unfixed = 0
            if k%2 == 0: #when k is an even number only the central values of the matrix are filled
                irun = [z.shape[0]-1-i for i in range(z.shape[0])]
                jrun = [z.shape[1]-1-j for j in range(z.shape[1])]
            else: #when k is an odd number the values at the borders are filled
                irun = range(z.shape[0]) 
                jrun = range(z.shape[1])
            jj = 0
            for i in irun:
                for j in jrun:
                    jj = jj+1
                    if z[i,j] < 0:
                        val = 0.; n = 0
                        if i > 0: 
                            if z[i-1,j] > 0: val = val+z[i-1,j]; n = n+1
                        if i < z.shape[0]-1: 
                            if z[i+1,j] > 0: val = val+z[i+1,j]; n = n+1
                        if j > 0: 
                            if z[i,j-1] > 0: val = val+z[i,j-1]; n = n+1
                        if j < z.shape[1]-1: 
                            if z[i,j+1] > 0: val = val+z[i,j+1]; n = n+1
                        if n>0:
                            z[i,j] = val/n; fixed = fixed+1
                        else:
                            unfixed = unfixed+1
                        fn = True
            print("...checked: {}, fixed: {}, unfixed: {}".format(jj,fixed,unfixed))       
       
    if exactcopy: #f function returns the z-coordinate of the pair (x,y)
        f = lambda x,y: z[int((yllcenter-y)/cs+nrows-1), int((x-xllcenter)/cs)]
    else: 
        f = interpolate.interp2d(x, y, z, kind='cubic')

    minx = min(x)
    maxx = max(x)
    miny = min(y)
    maxy = max(y)

    print("Limits of raster:")
    print("    x = [{},{}]".format(minx, maxx))
    print("    y = [{},{}]".format(miny, maxy))

    p1 = [minx, miny]
    p2 = [minx, maxy]
    p3 = [maxx, maxy]
    p4 = [maxx, miny]

    if p1in is not None:
        p1 = p1in+np.array([args.offsetx, args.offsety])
    if p2in is not None:
        p2 = p2in+np.array([args.offsetx, args.offsety])
    if p3in is not None:
        p3 = p3in+np.array([args.offsetx, args.offsety])
    if p4in is not None:
        p4 = p4in+np.array([args.offsetx, args.offsety])

    #It may be the case that the points introduced in maxY, minY, maxX or minX do not be completely contained in the original raster limits,
    #in which case they are modified to fulfill that condition.
    check_zone(maxY, minx, maxx, miny, maxy)
    check_zone(minY, minx, maxx, miny, maxy)
    check_zone(maxX, minx, maxx, miny, maxy)
    check_zone(minX, minx, maxx, miny, maxy)  

    pe = np.array([p1, p2, p3, p4])

    print("Meshing between:")
    print("    p1 = ({},{})".format(pe[0,0],pe[0,1]))
    print("    p2 = ({},{})".format(pe[1,0],pe[1,1]))
    print("    p3 = ({},{})".format(pe[2,0],pe[2,1]))
    print("    p4 = ({},{})".format(pe[3,0],pe[3,1]))
    s = 0.5*np.abs(np.dot(pe[:,0],np.roll(pe[:,1],1))-np.dot(pe[:,1],np.roll(pe[:,0],1)))
    print("Surface = {}".format(s))
    pm = [0.5*pe[0,0]+0.5*pe[2,0], 0.5*pe[0,1]+0.5*pe[2,1]]
    print("Point inside = ({},{},{})".format(pm[0], pm[1], 10+f(pm[0], pm[1])))

    #A zmax value is required for the top of the finite volume mesh. If the z-coordinate of the top is lower than the terrain coordinate
    #in a particular point, then pMesh will not be capable of creating the mesh.     
    if zmax == 0:
        zmax = (z.max()-z.min())*1.1

    ztop = [f(p1[0], p1[1])+zmax, 
            f(p2[0], p2[1])+zmax, 
            f(p3[0], p3[1])+zmax, 
            f(p4[0], p4[1])+zmax]
    f2 = interpolate.interp2d([p1[0], p2[0], p3[0], p4[0]],[p1[1], p2[1], p3[1], p4[1]], ztop) #f2 is a function that interpolates the top patch

##########################################  Mesh creation  #############################################################
  
    if exactcopy: #The grid of the original raster-file is used   
        if p1in is not None or p2in is not None or p3in is not None or p4in is not None:  #Only a fraction of the raster area is used          
            g = lambda p: [int((p[0]-xllcenter)/cs),int((yllcenter-p[1])/cs+nrows)]
            i1,j1 = g(p1)
            i2,j2 = g(p2)
            i3,j3 = g(p3)
            i4,j4 = g(p4)
            
            imin = min(i1,i2,i3,i4)
            imax = max(i1,i2,i3,i4)
            jmin = min(j1,j2,j3,j4)
            jmax = max(j1,j2,j3,j4)
            
            xres = imax-imin+1
            yres = jmax-jmin+1
            
            xll = x[imin]
            yll = y[jmax]
            
            xx = [float(xll+i*cs) for i in range(xres)]
            yy = [float(yll+(yres-i-1)*cs) for i in range(yres)] 
            yglob, xglob = np.meshgrid(yy, xx) 
        
        else: #The whole raster area is used
            yglob, xglob = np.meshgrid(y, x)
    else: #A different resolution is used with different limits
        if p1in is not None or p2in is not None or p3in is not None or p4in is not None:
            xloc = np.linspace(-1, 1, xres)
            yloc = np.linspace(-1, 1, yres)
            yyloc, xxloc = np.meshgrid(yloc, xloc)
    
            xglob = [] 
            yglob = []
    
            for i in range(len(xxloc)):
                xglob.append([])
                yglob.append([])
                for j in range(len(xxloc[i])):
                    p = getX(np.array([xxloc[i][j], yyloc[i][j]]), pe)
                    xglob[-1].append(p[0])
                    yglob[-1].append(p[1])

        #The mesh will not be a quadrilateral, being its shape defined by the points introduced in maxY, minY, maxX and minX           
        elif (maxY != [] or minY != [] or maxX != [] or minX != []): 
            xloc = np.linspace(0, 1, xres) 
            yloc = np.linspace(0, 1, yres)
                         
            maxY_points = get_boundary_points(maxY, xloc)
            minY_points = get_boundary_points(minY, xloc)
            maxX_points = get_boundary_points(maxX, yloc)
            minX_points = get_boundary_points(minX, yloc)
            
            xglob = [] 
            yglob = []
    
            for i in range(len(xloc)):
                xglob.append([])
                yglob.append([])
                for j in range(len(yloc)):
                    pw = minX_points[j]
                    pe = maxX_points[j]
                    ps = minY_points[i]
                    pn = maxY_points[i]
                    p = pw/(xloc[i]**2+10**-9)+pe/((1-xloc[i])**2+10**-9)+ps/(yloc[j]**2+10**-9)+pn/((1-yloc[j])**2+10**-9)
                    p /= 1/(xloc[i]**2+10**-9)+1/((1-xloc[i])**2+10**-9)+1/(yloc[j]**2+10**-9)+1/((1-yloc[j])**2+10**-9)
                    xglob[-1].append(p.x)
                    yglob[-1].append(p.y)


###########################################  Zones writing  #############################################################

    print("xres = {}".format(xres))
    print("yres = {}".format(yres))

    #File is written in the current directory
    stl = open(outfile+'.obj', 'w')
    stl.write('# Wavefront OBJ file written '+strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')
    stl.write("\n")
    
    nZones = 2
    if alphamaxY is not None:
        nZones = nZones + int(len(alphamaxY)/2)+1
    else:
        nZones = nZones + 1
    if alphaminY is not None:
        nZones = nZones + int(len(alphaminY)/2)+1
    else:
        nZones = nZones + 1
    if alphamaxX is not None:
        nZones = nZones + int(len(alphamaxX)/2)+1
    else:
        nZones = nZones + 1
    if alphaminX is not None:
        nZones = nZones + int(len(alphaminX)/2)+1
    else:
        nZones = nZones + 1
            
    nPoints = xres*yres*2
    nFaces = ((xres-1)*(yres-1)*2+(yres-1)*2+(xres-1)*2)*2
    
    stl.write('# points:  {} \n'.format(nPoints))
    stl.write('# faces:   {} \n'.format(nFaces))
    stl.write('# zones:   {} \n'.format(nZones))

    nfaces_terrain = (xres-1)*(yres-1)*2
    nfaces_x = (yres-1)*2
    nfaces_y = (xres-1)*2

    stl.write('#   0  terrain  (nFaces: {})\n'.format(nfaces_terrain))
    stl.write('#   1  maxZ  (nFaces: {})\n'.format(nfaces_terrain))

    face_number = 1
              
    if alphamaxX is not None:
        nfaces_xf = nfaces_x
        for i in range(int(len(alphamaxX)/2)):
            nfaces_i = int(round((alphamaxX[2*i+1]-alphamaxX[2*i])*nfaces_x,0))
            face_name = 'maxX'+str(i+1)
            face_number += 1
            stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_i))
            nfaces_xf = nfaces_xf - nfaces_i
        face_name = 'maxX'+str(int(len(alphamaxX)/2)+1)
        face_number += 1
        stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_xf))
    else:
        face_number += 1
        stl.write('#   {}  maxX  (nFaces: {})\n'.format(face_number,nfaces_x))

    if alphaminX is not None:
        nfaces_xf = nfaces_x
        for i in range(int(len(alphaminX)/2)):
            nfaces_i = int(round((alphaminX[2*i+1]-alphaminX[2*i])*nfaces_x,0))
            face_name = 'minX'+str(i+1)
            face_number += 1
            stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_i))
            nfaces_xf = nfaces_xf - nfaces_i
        face_name = 'minX'+str(int(len(alphaminX)/2)+1)
        face_number += 1
        stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_xf))
    else:
        face_number += 1
        stl.write('#   {}  minX  (nFaces: {})\n'.format(face_number,nfaces_x))

    if alphamaxY is not None:
        nfaces_yf = nfaces_y
        for i in range(int(len(alphamaxY)/2)):
            nfaces_i = int(round((alphamaxY[2*i+1]-alphamaxY[2*i])*nfaces_y,0))
            face_name = 'maxY'+str(i+1)
            face_number += 1
            stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_i))
            nfaces_yf = nfaces_yf - nfaces_i
        face_name = 'maxY'+str(int(len(alphamaxY)/2)+1)
        face_number += 1
        stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_yf))
    else:
        face_number += 1
        stl.write('#   {}  maxY  (nFaces: {})\n'.format(face_number,nfaces_y))

    if alphaminY is not None:
        nfaces_yf = nfaces_y
        for i in range(int(len(alphaminY)/2)):
            nfaces_i = int(round((alphaminY[2*i+1]-alphaminY[2*i])*nfaces_y,0))
            face_name = 'minY'+str(i+1)
            face_number += 1
            stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_i))
            nfaces_yf = nfaces_yf - nfaces_i
        face_name = 'minY'+str(int(len(alphaminY)/2)+1)
        face_number += 1
        stl.write('#   {}  {}  (nFaces: {})\n'.format(face_number, face_name,nfaces_yf))
    else:
        face_number += 1
        stl.write('#   {}  minY  (nFaces: {})\n'.format(face_number,nfaces_y))

###########################################  Points writing  #############################################################

    stl.write("\n")
    stl.write('# <points count="{}"> \n'.format(nPoints))        
    print("writing vertexes...")
    
    l = 0
    ls = []
    for i in range(xres):
        ls.append([])
        for j in range(yres):
            ls[-1].append([])
            ls[-1][-1].append(l)
            l += 1
            stl.write("v {} {} {} \n".format(round(xglob[i][j],3), round(yglob[i][j],3),round(float(f(xglob[i][j], yglob[i][j])),2)))
            ls[-1][-1].append(l)
            l += 1
            stl.write("v {} {} {} \n".format(round(xglob[i][j],3), round(yglob[i][j],3),round(float(f2(xglob[i][j], yglob[i][j])),2)))
    stl.write("# </points>\n")
    stl.write("\n")

###########################################  Faces writing  #############################################################

    print("writing terrain faces...")
    stl.write('# <faces count="{}">\n'.format(nFaces))
    stl.write('g terrain\n')
    
    for i in range(xres-1):
        for j in range(yres-1):
            if (i+j)%2==0:
                stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i+1][j][0]+1,ls[i+1][j+1][0]+1))
                stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i+1][j+1][0]+1,ls[i][j+1][0]+1))
            else:
                stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i+1][j][0]+1,ls[i][j+1][0]+1))
                stl.write("f {0} {2} {1} \n".format(ls[i][j+1][0]+1,ls[i+1][j][0]+1,ls[i+1][j+1][0]+1))
    
    print("writing g maxZ faces...")
    stl.write('g maxZ\n')
    
    for i in range(xres-1):
        for j in range(yres-1):
            stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j+1][1]+1))
            stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i+1][j+1][1]+1,ls[i][j+1][1]+1))

    if alphamaxX is not None:
        for k in range(int(len(alphamaxX)/2)):
            face_name = 'maxX'+str(k+1)
            print("writing g "+ face_name + " faces...")
            stl.write('g {}\n'.format(face_name))
    
            for j in range(yres-1):
                if j/(yres-1)>=alphamaxX[2*k] and j/(yres-1)<alphamaxX[2*k+1] :
                    i = xres-1
                    stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
                    stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))
    

        face_name = 'maxX'+str(int(len(alphamaxX)/2)+1)
        print("writing g "+ face_name + " faces...")
        stl.write('g {}\n'.format(face_name))
        
        for j in range(yres-1):
            t = 1
            for k in range(int(len(alphamaxX)/2)):
                if j/(yres-1)>=alphamaxX[2*k] and j/(yres-1)<alphamaxX[2*k+1] :
                    t = 0
                
            if t ==1:
                i = xres-1
                stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
                stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))

    else:
        print("writing g maxX faces...")
        stl.write('g maxX\n')
        
        for j in range(yres-1):
            i = xres-1
            stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
            stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))


    if alphaminX is not None:
        for k in range(int(len(alphaminX)/2)):
            face_name = 'minX'+str(k+1)
            print("writing g "+ face_name + " faces...")
            stl.write('g {}\n'.format(face_name))
    
            for j in range(yres-1):
                if j/(yres-1)>=alphaminX[2*k] and j/(yres-1)<alphaminX[2*k+1] :
                    i = 0
                    stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
                    stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))
   

        face_name = 'minX'+str(int(len(alphaminX)/2)+1)
        print("writing g "+ face_name + " faces...")
        stl.write('g {}\n'.format(face_name))
        
        for j in range(yres-1):
            t = 1
            for k in range(int(len(alphaminX)/2)):
                if j/(yres-1)>=alphaminX[2*k] and j/(yres-1)<alphaminX[2*k+1] :
                    t = 0
                
            if t ==1:
                i = 0
                stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
                stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))

    else:
        print("writing g minX faces...")
        stl.write('g minX\n')
        
        for j in range(yres-1):
            i = 0
            stl.write("f {0} {2} {1} \n".format(ls[i][j][0]+1,ls[i][j+1][0]+1,ls[i][j][1]+1))
            stl.write("f {0} {2} {1} \n".format(ls[i][j][1]+1,ls[i][j+1][0]+1,ls[i][j+1][1]+1))

    if alphamaxY is not None:
        for k in range(int(len(alphamaxY)/2)):
            face_name = 'maxY'+str(k+1)
            print("writing g "+ face_name + " faces...")
            stl.write('g {}\n'.format(face_name))
    
            for i in range(xres-1):
                if i/(xres-1)>=alphamaxY[2*k] and i/(xres-1)<alphamaxY[2*k+1] :
                    j = yres-1
                    stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
                    stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))
          

        face_name = 'maxY'+str(int(len(alphamaxY)/2)+1)
        print("writing g "+ face_name + " faces...")
        stl.write('g {}\n'.format(face_name))
        
        for i in range(xres-1):
            t = 1
            for k in range(int(len(alphamaxY)/2)):
                if i/(xres-1)>=alphamaxY[2*k] and i/(xres-1)<alphamaxY[2*k+1] :
                    t = 0
                
            if t ==1:
                j = yres-1
                stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
                stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))

    else:
        print("writing g maxY faces...")
        stl.write('g maxY\n')
        
        for i in range(xres-1):
            j = yres-1
            stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
            stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))
 

    if alphaminY is not None:
        for k in range(int(len(alphaminY)/2)):
            face_name = 'minY'+str(k+1)
            print("writing g "+ face_name + " faces...")
            stl.write('g {}\n'.format(face_name))
    
            for i in range(xres-1):
                if i/(xres-1)>=alphaminY[2*k] and i/(xres-1)<alphaminY[2*k+1] :
                    j = 0
                    stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
                    stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))

        face_name = 'minY'+str(int(len(alphaminY)/2)+1)
        print("writing g "+ face_name + " faces...")
        stl.write('g {}\n'.format(face_name))
        
        for i in range(xres-1):
            t = 1
            for k in range(int(len(alphaminY)/2)):
                if i/(xres-1)>=alphaminY[2*k] and i/(xres-1)<alphaminY[2*k+1] :
                    t = 0
                
            if t ==1:
                j = 0
                stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
                stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))

    else:
        print("writing g minY faces...");
        stl.write('g minY\n')
        
        for i in range(xres-1):
            j = 0
            stl.write("f {} {} {} \n".format(ls[i][j][0]+1,ls[i][j][1]+1,ls[i+1][j][0]+1))
            stl.write("f {} {} {} \n".format(ls[i][j][1]+1,ls[i+1][j][1]+1,ls[i+1][j][0]+1))

    stl.write('# </faces>')  
    stl.close()

#########################################################################################################################

if __name__ == "__main__":
   main(sys.argv[1:])