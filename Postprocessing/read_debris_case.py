#!/usr/bin/env python3

'''
Description
    This code is used to process the results of a debrisfaSavageHutterFoam simulation. If the simulation was performed
    in parallel then it is necessary to reconstruct the fields with the reconstructPar utility or the Par_reconstructPar.py.
Author
    Álvaro González Bilbao alvaro.gonzalez.bilbao@gmail.com

@version: 1.31 (21/06/21)
'''

import numpy as np
import sys  
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import matplotlib.colors as mcolors
import os
import math
from operator import itemgetter
import argparse
import shutil

##########################################  Definitions  ##############################################################

def clean_brackets(list):
    k = 0
    while True:
        repeat = 'False'
        if list[k]=='[':
            list.remove('[')
            repeat = 'True'
        elif list[k]==']':
            list.remove(']')
            repeat = 'True'
        elif list[k][-1]==']':
            list[k] = list[k][:-1]
            repeat = 'True'    
        elif list[k][0]=='[':
            list[k] = list[k][1:]
            repeat = 'True'  
        elif list[k]=='(':
            list.remove('(')
            repeat = 'True'
        elif list[k]==')':
            list.remove(')')
            repeat = 'True'
        elif list[k][-1]==')':
            list[k] = list[k][:-1]
            repeat = 'True'    
        elif list[k][0]=='(':
            list[k] = list[k][1:]
            repeat = 'True' 
        if repeat == 'False':
            k = k+1
        if k == len(list):
            break

def clean_dots(list):
    k = 0
    while True:
        repeat = 'False'
        if list[k]==';':
            list.remove(';')
            repeat = 'True'
        elif list[k][-1]==';':
            list[k] = list[k][:-1]
            repeat = 'True'
        if repeat == 'False':
            k = k+1
        if k == len(list):
            break

def clean_spaces_string(list):
    k = 0
    while True:
        if list[k] == '\n':
            list.remove('\n')
        elif list[k].split()[0] == '//':
            list.remove(list[k])
        else: 
            k += 1
        if k == len(list):
            break

def clean_spaces_list(list):
    k = 0
    while True:
        if len(list[k]) == 0:
            list.remove([])
        else: 
            k += 1
        if k == len(list):
            break

def read_entry(list,dict):
    dict[list[0]] = list[1]
    
def read_dimensionedScalar(value,dict):
    dict['name']=value[0]
    dimensions_list = []
    for j in range(7):
        dimensions_list.append(float(value[1+j]))
    dict['dimensions']=dimensions_list
    dict['value']=float(value[-1])
    
def read_vector(value,list):
    clean_brackets(value)
    if len(value)==3:
        for i in range(len(value)):
            value[i] = float(value[i])
    list.append(value)

def read_constant(list, dict):
    value = list[1:]
    clean_brackets(value)
    subdict = {}
    sublist = []
    if len(value)== 9:
        read_dimensionedScalar(value,subdict)
        dict[list[0]] = subdict
    else:
        read_vector(value,sublist)
        dict[list[0]] = sublist
        
def read_subdictionary(list,output,name,type):
    subdict = {}
    while True:
        pop = list.pop(0)
        if pop[0] == '}':
            break
        clean_brackets(pop)
        if len(pop)==2:
            read_entry(pop,subdict)
        elif len(pop)>2:
            read_constant(pop,subdict)
        else:
            if list[0][0]=='{':
                list.remove(list[0])
                read_subdictionary(list,subdict,pop[0],'dict')
            elif list[0][0]=='(':
                list.remove(list[0])
                read_sublist(list,subdict,pop[0],'dict')
    if type == 'dict':
        output[name] = subdict
    elif type == 'list':
        output.append(name)
        output.append(subdict)

def read_sublist(list,output,name,type):
    sublist = []
    while True:
        pop = list.pop(0)
        if pop[0] == ')':
            break
        clean_brackets(pop)
        if len(pop)==1:
            if list[0][0]=='{':
                list.remove(list[0])
                read_subdictionary(list,sublist,pop[0],'list')
            else:
                sublist.append(pop[0])
        elif len(pop) > 1:
            read_vector(pop,sublist)
    if type == 'dict':
        output[name] = sublist
    elif type == 'list':
        output.append(name)
        output.append(sublist)

def read_dictionary(list,dict):
    while True:
        pop = list.pop(0)
        if len(pop)==2:
            read_entry(pop,dict)
        elif len(pop)==1:
            if list[0][0]=='{':
                list.remove(list[0])
                read_subdictionary(list,dict,pop[0],'dict')
            elif list[0][0]=='(':
                list.remove(list[0])
                read_sublist(list,dict,pop[0],'dict')
        elif len(pop)>2:
            read_constant(pop,dict)          
        if len(list) == 0:
            break
        
def read_face(List):
    subList = []
    if len(List) < 10:
        subList.append(int(List[0][2:]))
    else:
        subList.append(int(List[0][3:]))
    for i in range(len(List)-2):
        subList.append(int(List[i+1]))
    subList.append(int(List[-1][:-1]))
    return subList
       
def get_input(path,list, OF_header = True):
    input = []
    with open(path) as f:
        for line in f:
            input.append(line)
    if OF_header == True:
        input= input[16:(len(input))]
    clean_spaces_string(input)
    for i in range(len(input)):
        list.append(input[i].split())
        clean_dots(list[i])
    clean_spaces_list(list)
    
def get_input_number(path,list):
    input = []
    with open(path) as f:
        i = 0
        for line in f:
            input.append(line)
            i += 1
            if i == 19:
                break
    input= input[16:(len(input))]
    clean_spaces_string(input)
    for i in range(len(input)):
        list.append(input[i].split())
        clean_dots(list[i])
    clean_spaces_list(list)

def get_number_faces(path):  
    list = []
    get_input_number(path+'/0/faFaces',list)  
    return int(list[0][0])

def get_number_edges(path):  
    list = []
    get_input_number(path+'/0/edgeNeighbour',list)  
    return int(list[0][0])

def read_transportProperties(path,dict): 
    list = []
    get_input(path+'/constant/transportProperties',list)    
    read_dictionary(list,dict)
    
def read_releaseArea(path,dict):
    list = []
    get_input(path+'/constant/releaseArea',list)    
    read_dictionary(list,dict)

def read_releaseFlow(path,dict):
    list = []
    get_input(path+'/constant/releaseFlow',list)    
    read_dictionary(list,dict)

def read_controlDict(path,dict): 
    list = []
    get_input(path+'/system/controlDict',list)    
    read_dictionary(list,dict)

def read_meshDict(path,dict): 
    list = []
    get_input(path+'/system/meshDict',list)    
    read_dictionary(list,dict)

def correct_EdgeLabels(list):
    new_list = []
    for i in range(len(list)):
        if (list[i][0] == 'edgeLabels' and len(list[i])>2):
            new_list.append(list[i][:2])
            problem_list = list[i][2:]
            index = problem_list[0].index('(')
            new_list.append([problem_list[0][:index]])
            new_list.append(['('])
            if len(problem_list) > 2:
                new_list.append([problem_list[0][(index+1):]])
                for i in range(len(problem_list)-2):
                    new_list.append([problem_list[i+1]])
                new_list.append([problem_list[-1][:-1]])
            else:
                new_list.append([problem_list[0][(index+1):-1]])
            new_list.append([')'])
        else:
            new_list.append(list[i])
    return new_list
    
def read_faBoundary(path,dict):  
    list = []
    get_input(path+'/constant/faMesh/faBoundary',list)
    list = correct_EdgeLabels(list)
    subdict = {}
    read_dictionary(list,subdict)
    for key in subdict.keys():
        subsubdict = {}
        n = 0
        for j in range(int(len(subdict[key])/2)):
            for sub_key in subdict[key][j*2+1].keys():
                k = str(sub_key)
                if k.isnumeric() == True:
                    n += 1
                    sub_list = subdict[key][j*2+1][k]
                    for i in range(len(sub_list)):
                        sub_list[i] = int(sub_list[i])
                    subsubdict[subdict[key][j*2]] = sub_list
        dict['boundaries'] = subsubdict
        dict['number'] = n
    
def read_faces(path,dict):  
    list = []
    get_input(path+'/0/faFaces',list)  
    faces = len(list)-3
    subList = []
    for i in range(faces):
        subList.append(read_face(list[i+2]))    
    dict['faces'] = subList
    dict['number'] = faces

def read_points(path,dict):  
    list = []
    get_input(path+'/0/faPoints',list)
    subdict = {}
    points = str(len(list)-3)
    read_dictionary(list,subdict)
    dict['points'] = subdict[points]
    dict['number'] = int(points)

def read_edgeOwner(path,dict):  
    list = []
    get_input(path+'/0/edgeOwner',list)
    subdict = {}
    edges = str(len(list)-3)
    read_dictionary(list,subdict)
    dict['edges'] = subdict[edges]
    dict['number'] = int(edges)
    
def read_edgeNeighbour(path,dict):  
    list = []
    get_input(path+'/0/edgeNeighbour',list)    
    subdict = {}
    edges = str(len(list)-3)
    read_dictionary(list,subdict)
    dict['edges'] = subdict[edges]
    dict['number'] = int(edges)

def number_processors(path):
    listdir = os.listdir(path)
    n_processors = 0
    for i in range(len(listdir)):
        if listdir[i][:9] == 'processor':
            if n_processors < int(listdir[i][9:]):
                n_processors = int(listdir[i][9:])
    return (n_processors+1)
    
def create_time(path, time):
    del time [:]
    listdir = os.listdir(path)
    k = 0
    while True:
        if is_number(listdir[k])==True:
            if float(listdir[k]) == int(float(listdir[k])):
                time.append(int(float(listdir[k])))
            else:
                time.append(float(listdir[k]))
            k += 1
        else:
            listdir.remove(listdir[k])
        if k == len(listdir):
            break
    time.sort()

def read_proc_number_faces(path, n_proc, output):
    for p in range(n_proc):
        list = []
        get_input_number(path+'/processor'+str(p)+'/constant/faMesh/faceLabels',list) 
        output['processor'+str(p)] = int(list[0][0])

def read_proc_faBoundary(path, n_proc, output):
    for p in range(n_proc):
        dict = {}
        read_faBoundary(path+'/processor'+str(p), dict)
        output['processor'+str(p)] = dict

def read_proc_faceaddr(path, n_proc, output):
    for p in range(n_proc):
        list = []
        get_input(path+'/processor'+str(p)+'/constant/faMesh/faceProcAddressing',list)
        list = list[2:-1]
        for i in range(len(list)):
            list[i] = int(list[i][0])
        output['processor'+str(p)] = list               

def read_proc_edgeaddr(path, n_proc, output):
    for p in range(n_proc):
        list = []
        get_input(path+'/processor'+str(p)+'/constant/faMesh/edgeProcAddressing',list)
        list = list[2:-1]
        for i in range(len(list)):
            list[i] = int(list[i][0])
        output['processor'+str(p)] = list 

def read_proc_edgeOwners(path, n_proc, dict):                    
    for p in range(n_proc):
        list = []
        get_input(path+'/processor'+str(p)+'/0/edgeOwner',list)
        list = list[2:-1]
        for i in range(len(list)):
            list[i] = int(list[i][0])
        dict['processor'+str(p)] = list

def read_proc_number_edges(path, n_proc, output, proc_fB):
    for p in range(n_proc):
        dict = {}
        list = []
        get_input_number(path+'/processor'+str(p)+'/constant/faMesh/edgeProcAddressing',list) 
        dict['Total'] = int(list[0][0])
        sum = 0
        for key in proc_fB['processor'+str(p)]['boundaries']:
            sum += len(proc_fB['processor'+str(p)]['boundaries'][key])
        dict['internal'] = dict['Total'] - sum
        output['processor'+str(p)] = dict
                       
def read_lp(path,dict):
    input = []
    with open(path) as f:
        for line in f:
            input.append(line)
    clean_spaces_string(input)
    list = []
    for i in range(len(input)):
        list.append(input[i].split())
        clean_dots(list[i])
    clean_spaces_list(list)
    
    n_Channels = int(list[0][1])
    
    for i in range(n_Channels):
        subdict = {}
        name = str(list[1+5*i][0])
        x_list = list[3+5*i]
        y_list = list[4+5*i]
        
        for i in range(len(x_list)):
            x_list[i] = float(x_list[i])
            y_list[i] = float(y_list[i])
        
        subdict['x_list'] = x_list
        subdict['y_list'] = y_list
        dict[name] = subdict

def read_tp(path,dict):
    input = []
    with open(path) as f:
        for line in f:
            input.append(line)
    clean_spaces_string(input)
    list = []
    for i in range(len(input)):
        list.append(input[i].split())
        clean_dots(list[i])
    clean_spaces_list(list)
    
    n_Channels = int(list[0][1])
    
    for i in range(n_Channels):
        name = str(list[1+13*i][0])
        side1 = str(list[3+13*i][0])
        x_list_side1 = list[5+13*i]
        y_list_side1 = list[6+13*i]
        side2 = str(list[8+13*i][0])
        x_list_side2 = list[10+13*i]
        y_list_side2 = list[11+13*i]
        
        list_tp = []
        for i in range(len(x_list_side1)):
            if side1 == 'Left':
                vp1 = vector(float(x_list_side1[i]), float(y_list_side1[i]), 0)
                vp2 = vector(float(x_list_side2[i]), float(y_list_side2[i]), 0)
            elif side2 == 'Left':
                vp1 = vector(float(x_list_side2[i]), float(y_list_side2[i]), 0)
                vp2 = vector(float(x_list_side1[i]), float(y_list_side1[i]), 0)                
            list_tp.append([vp1,vp2])
            
        dict[name] = list_tp

class areaField:
        def __init__(self,type,internalField, boundaryField, dimension):
                self.ty = type
                self.iF= internalField
                self.bF= boundaryField
                self.dim = dimension
                                                
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

def distxy(v1,v2):
        return round(math.sqrt((v1.x-v2.x)**2+(v1.y-v2.y)**2),6)

def dist(v1,v2):
        return round(math.sqrt((v1.x-v2.x)**2+(v1.y-v2.y)**2+(v1.z-v2.z)**2),6)

class dimensions:
        def __init__(self,List):
                self.M= float(List[0][0])
                self.L= float(List[0][1])
                self.T= float(List[0][2])
                self.Te= float(List[0][3])
                self.Mo= float(List[0][4])
                self.C= float(List[0][5])
                self.I= float(List[0][6])
              
        def List(self):
                return [self.M,self.L,self.T,self.Te,self.Mo,self.C,self.I]
    
class internalField:
        def __init__(self,dict,type,number, correct = True):
                self.field = []
                self.u = ''
                self.correct_input(dict,type,number, correct)
              
        def correct_input(self,dict,type,n, correct):
                if correct == True:
                    self.u = dict['internalField'][0][0]
                    if self.u == 'nonuniform':
                        list = dict[str(n)] 
                        if type == 'scalar':
                            for i in range(n):
                                self.field.append(float(list[i]))
                        elif type == 'vector':
                            for i in range(n):
                                self.field.append(vector(list[i][0],list[i][1],list[i][2]))
                    elif self.u == 'uniform':
                        if type == 'scalar':
                            self.field.append(float(dict['internalField'][0][1]))
                        elif type == 'vector':
                            self.field.append(vector(dict['internalField'][0][1:][0],dict['internalField'][0][1:][1],dict['internalField'][0][1:][2]))
                else:
                    self.u = dict['u']
                    self.field = dict['field']

def correct_bFEdgeLabels(dict, type):
        new_dict = {}
        for key in dict.keys():
            if dict[key]['type'] == 'calculated' or dict[key]['type'] == 'fixedValue' or dict[key]['type'] == 'processor':
                if dict[key]['value'][0][1] == 'List<'+type+'>' and len(dict[key]['value'][0])>2:
                    problem_list = dict[key]['value'][0][2:]
                    sub_dict = {}
                    sub_dict['type'] = dict[key]['type']
                    sub_dict['value'] = [dict[key]['value'][0][:2]]
                    if type == 'vector':
                        sub_list = []
                        for i in range(int(len(problem_list)/3)):
                            if i > 0:
                                sub_list.append([float(problem_list[3*i]),float(problem_list[3*i+1]),float(problem_list[3*i+2])])
                            else:
                                index = problem_list[3*i].index('(')+2
                                sub_list.append([float(problem_list[3*i][index:]),float(problem_list[3*i+1]),float(problem_list[3*i+2])])
                        sub_dict[str(int(len(problem_list)/3))] = sub_list
                    elif type == 'scalar':
                        sub_list = []
                        for i in range(len(problem_list)):
                            if i > 0:
                                sub_list.append(float(problem_list[i]))
                            else:
                                index = problem_list[i].index('(')+1
                                sub_list.append(float(problem_list[i][index:]))
                        sub_dict[str(len(problem_list))] = sub_list
                    new_dict[key] = sub_dict
                else:
                    new_dict[key] = dict[key]
            else:
                new_dict[key] = dict[key]
        return new_dict

class boundaryField:
        def __init__(self,dict,type,faBoundary, correct = True):
                self.d= dict
                self.correct_input(type,faBoundary, correct)
                         
        def correct_input(self,type,faBoundary, correct):
                if correct == True:
                    self.d = correct_bFEdgeLabels(self.d, type)
                    pop_list = []
                    for key in self.d.keys():
                        if (key in faBoundary['boundaries'].keys()) == True:
                            if self.d[key]['type'] == 'fixedValue' or self.d[key]['type'] == 'processor' or self.d[key]['type'] == 'calculated':
                                if self.d[key]['value'][0][0] == 'uniform':
                                    value = self.d[key].pop('value')[0]
                                    subdict = {}
                                    subdict['uniform'] = 'yes'
                                    sublist = []
                                    if faBoundary['boundaries'][key] != 0:
                                        if type == 'vector':
                                            v = vector(value[1:][0],value[1:][1],value[1:][2])
                                        elif type == 'scalar':
                                            v= float(value[1:][0])
                                        sublist.append(v)                            
                                    subdict['List'] = sublist
                                    self.d[key]['value'] = subdict
                                elif self.d[key]['value'][0][0] == 'nonuniform':
                                    pop = self.d[key].pop('value')
                                    subdict = {}
                                    subdict['uniform'] = 'no'
                                    sublist = []
                                    if faBoundary['boundaries'][key] != 0:
                                        value = self.d[key].pop(str(len(faBoundary['boundaries'][key])))
                                        for i in range(len(value)):
                                            if type == 'vector':
                                                v = vector(value[i][0],value[i][1],value[i][2])
                                            elif type == 'scalar':
                                                v= float(value[i])
                                            sublist.append(v)
                                    subdict['List'] = sublist
                                    self.d[key]['value'] = subdict
                        else:
                            pop_list.append(key)
                    for key in pop_list:
                        self.d.pop(key)

#Right now, read_areaField and read_edgeField are the same, probably in the past there were some differences. I will conserve both for the moment.  //AG                         
def read_areaField(path,time,name,type,output,number_faces,faBoundary):
    list = []
    get_input(path+'/'+str(time)+'/'+name,list)
    field = {}
    read_dictionary(list,field)
    dimension = dimensions(field['dimensions'])
    boundaryfield = boundaryField(field['boundaryField'],type,faBoundary)
    internalfield = internalField(field,type,number_faces)
    output[str(time)] = areaField(type,internalfield,boundaryfield,dimension)

def read_edgeField(path,time,name,type,output,number_edges,faBoundary):
    list = []
    get_input(path+'/'+str(time)+'/'+name,list)
    field = {}
    read_dictionary(list,field)
    dimension = dimensions(field['dimensions'])
    boundaryfield = boundaryField(field['boundaryField'],type,faBoundary)
    internalfield = internalField(field,type,number_edges)
    output[str(time)] = areaField(type,internalfield,boundaryfield,dimension)
       
def read_h(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'h','scalar',output, number_faces,faBoundary)
                
def read_pb(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'pb','scalar',output, number_faces,faBoundary)
        
def read_Cv(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'Cv','scalar',output, number_faces,faBoundary)
                
def read_deltac0(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'deltac0','scalar',output, number_faces,faBoundary)

def read_deltah0(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'deltah0','scalar',output, number_faces,faBoundary)

def read_Us(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'Us','vector',output, number_faces,faBoundary)

def read_tau(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'tau','vector',output, number_faces,faBoundary)
        
def read_n(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'n','vector',output, number_faces ,faBoundary)

def read_he(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'he','scalar',output, number_faces,faBoundary)

def read_c(path, time, output, number_faces,faBoundary):
    read_areaField(path, time, 'c','vector',output, number_faces ,faBoundary)

def read_ec(path,output, number_edges,faBoundary):
    read_edgeField(path, 0, 'ec','vector',output, number_edges,faBoundary)

def read_A(path,time, output, number_faces,faBoundary):
    read_areaField(path, time, 'A','scalar',output, number_faces ,faBoundary)
        
def read_phi2s(path, time, output, number_edges,faBoundary):
    read_edgeField(path, time, 'phi2s','scalar',output, number_edges ,faBoundary)
        
def read_Q(path, time, output, number_edges, faBoundary):
    read_edgeField(path, time, 'Q','scalar',output, number_edges, faBoundary)

def read_c_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'c','vector',output, n_processors, proc_fB, proc_nF)
        
def read_ec_proc(path,output,n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, 0, 'ec','vector',output, n_processors, proc_fB, proc_nE) 
       
def read_Q_proc(path, time, output, n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, time, 'Q', 'scalar', output, n_processors, proc_fB, proc_nE)

def read_phi2s_proc(path, time, output, n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, time, 'phi2s', 'scalar', output, n_processors, proc_fB, proc_nE)

def read_areaField_proc(path, time, name, type, output, n_processors, proc_fB, proc_nF):
    dict = {}
    for p in range(n_processors):
        list = []
        get_input(path+'/processor'+str(p)+'/'+str(time)+'/'+name,list)
        field = {}
        read_dictionary(list,field)
        dimension = dimensions(field['dimensions'])
        boundaryfield = boundaryField(field['boundaryField'],type,proc_fB['processor'+str(p)])
        internalfield = internalField(field,type,proc_nF['processor'+str(p)])
        dict['processor'+str(p)] = areaField(type,internalfield,boundaryfield,dimension)
    output[str(time)] = dict

def read_edgeField_proc(path, time, name, type, output, n_processors, proc_fB, proc_nE):
    dict = {}
    for p in range(n_processors):
        list = []
        get_input(path+'/processor'+str(p)+'/'+str(time)+'/'+name,list)
        field = {}
        read_dictionary(list,field)
        dimension = dimensions(field['dimensions'])
        boundaryfield = boundaryField(field['boundaryField'],type,proc_fB['processor'+str(p)])
        internalfield = internalField(field,type,proc_nE['processor'+str(p)]['internal'])
        dict['processor'+str(p)] = areaField(type,internalfield,boundaryfield,dimension)
    output[str(time)] = dict

#When running the code in parallel the reconstructed fields are created in ascending order,
#that is why for an edge shared by proc0 and proc1 the value that will remain is the one in proc1 //AG
def correct_edge_proc(internal_field, procBoundary_owner, procBoundary_edgeList):                       
    if internal_field.u == 'nonuniform':
        for patch_key in procBoundary_owner.keys():
            for i in range(len(procBoundary_owner[patch_key])):
                first_owner = procBoundary_owner[patch_key][i]
                if first_owner == 1:
                    internal_field.field[procBoundary_edgeList[patch_key][i]] *= -1


##########################################  runCase class  ##############################################################

#Class runCase is used to read the outputs of the OpenFOAM simulation. Just a few calculations are performed inside this class //AG  
class runCase:
        def __init__(self,path,lp_name='',tp_name='',h='on',pb='off',Cv='on',deltaz0='on',Us='on',tau='off',phi2s='off',Q='on',c='off',n='off',he='off', A='off'):
                self.t= []
                self.h = {}               
                self.pb = {}
                self.Cv = {}
                self.deltah0 = {}
                self.deltac0 = {}
                self.Us = {}
                self.tau = {}
                self.Q = {}
                self.phi2s = {}
                self.n = {}
                self.he = {}
                self.c = {}
                self.A = {}
                
                self.p = path
                self.nF = get_number_faces(self.p)
                self.nE = {}  #number Edges 
                self.np = number_processors(self.p)
                self.ec = {}  #edgeCentres
                self.eO = {}  #edgeOwner
                self.eN = {}  #edgeNeighbor
                self.faces  = {}
                self.points = {}
                self.fB = {}  #faBoundary
                self.lp = {}  #longitudinal profile
                self.tp = {}  #transversal profile
                                
                self.tps = {}  #transportProperties dictionary
                self.rA  = {}  #releaseArea dictionary
                self.rF  = {}  #releaseFlow dictionary
                self.mD  = {}  #meshDict dictionary               
                self.lp_name = lp_name
                self.tp_name = tp_name
                self.app = ''
                self.sT  = ''
                self.eT  = ''
                self.wI  = ''
                self.maxCo = ''
                self.hmin  = -1
                self.rho_w = -1
                self.rho_s = -1
                self.rho_b = -1
                self.eM = ''
                self.eM_coeffs = {}
               
                self.h_flag  = h
                self.pb_flag = pb
                self.Cv_flag = Cv
                self.deltaz0_flag = deltaz0
                self.Us_flag  = Us
                self.tau_flag = tau
                self.Q_flag = Q
                self.phi2s_flag = phi2s
                self.c_flag  = c
                self.A_flag  = A
                self.n_flag  = n
                self.he_flag = he
                
                self.xi = 0
                self.xf = 0
                self.yi = 0
                self.yf = 0
                
                self.proc_owner = {}
                self.proc_edgeList = {}
                self.proc_fB = {}  #processors faBoundary
                self.proc_nE = {}  #processors number of edges
                self.proc_nF = {}
                self.proc_ec = {}
                self.proc_edgeaddr = {}
                self.proc_faceaddr = {}
                self.proc_edgeOwners = {}
                
                self.get_case()

        def get_case(self):
                self.get_faBoundary()
                self.get_nEdges()
                self.get_faces()
                self.get_points()
                self.get_lp()
                self.get_time()
                self.get_ec()
                self.get_limits()
                self.get_edgeOwner()
                self.get_edgeNeighbour()
                self.get_parallel_data()
                self.get_hmin()
                self.get_n()
                self.get_c()
                self.get_A()
            
        def get_faBoundary(self, report = True):
                if report == True:
                    print ('Reading faBoundary ...')
                read_faBoundary(self.p,self.fB)
                
        def get_nEdges(self):
                self.nE['internal'] = get_number_edges(self.p)
                sum = 0
                for key in self.fB['boundaries']:
                    sum += len(self.fB['boundaries'][key])
                self.nE['Total'] = self.nE['internal'] + sum
                          
        def get_faces(self, report = True):
                if report == True:
                    print('Reading faces ...')
                read_faces(self.p,self.faces)
                
        def get_points(self, report = True):
                if report == True:
                    print('Reading points ...')
                read_points(self.p,self.points)
                         
        def get_lp(self, report = True):
                if self.lp_name != '':
                    if report == True:
                        print ('Reading longitudinal profiles ...')
                    read_lp(self.p+'/'+self.lp_name, self.lp)                

        def get_tp(self, report = True):
                if self.lp_name != '':
                    if report == True:
                        print ('Reading transversal profiles ...')
                    read_tp(self.p+'/'+self.tp_name, self.tp) 
         
        def get_time(self, report = True):
                if report == True:
                    print ('Reading times ...')
                create_time(self.p, self.t)

        def get_ec(self, report = True):
                if report == True:
                    print('Reading ec field ...')
                read_ec(self.p,self.ec, self.nE['internal'], self.fB)
                field = [vector(-9999,-9999,-9999)]*self.nE['Total']
                for i in range(self.nE['internal']):
                    field[i] = self.ec['0'].iF.field[i]
                for patch_key in self.ec['0'].bF.d.keys():
                    list = self.ec['0'].bF.d[patch_key]['value']['List']
                    for i in range(len(list)):
                        field[self.fB['boundaries'][patch_key][i]] = list[i]
                self.ec['0'].iF.field = field
                
        def get_limits(self, report = True):
                if report == True:
                    print('Reading region limits ...')
                if self.ec == {}:
                    self.get_ec(False)

                self.xi = math.inf
                self.xf = -math.inf
                self.yi = math.inf
                self.yf = -math.inf
                               
                dict = self.ec['0'].bF.d
                for key in dict.keys():
                    list = dict[key]['value']['List']
                    sublist_x = []
                    sublist_y = []
                    for i in range(len(list)):
                        sublist_x.append(list[i].x)
                        sublist_y.append(list[i].y)
                    min_x = min(sublist_x)
                    max_x = max(sublist_x)
                    min_y = min(sublist_y)
                    max_y = max(sublist_y)

                    if min_x <= self.xi:
                        self.xi = min_x
                    elif max_x >= self.xf:
                        self.xf = max_x
                    if min_y <= self.yi:
                        self.yi = min_y
                    elif max_y >= self.yf:
                        self.yf = max_y

        def get_hmin(self, report = True):
                if report == True:
                    print('Reading hmin ...')
                dict = {}
                read_transportProperties(self.p,dict)
                self.hmin = dict['hmin']['value']

        def get_edgeOwner(self, report = True):
                if report == True:
                    print('Reading edgeOwners ...')
                read_edgeOwner(self.p,self.eO)
            
        def get_edgeNeighbour(self, report = True):
                if report == True:
                    print('Reading edgeNeighbours ...')
                read_edgeNeighbour(self.p,self.eN) 

        def get_proc_faBoundary(self):          
                read_proc_faBoundary(self.p, self.np, self.proc_fB)

        def get_proc_nF(self):
                read_proc_number_faces(self.p, self.np, self.proc_nF)

        def get_proc_nE(self):
                read_proc_number_edges(self.p, self.np, self.proc_nE, self.proc_fB)
            
        def get_proc_ec(self):
                read_ec_proc(self.p,self.proc_ec,self.np, self.proc_fB, self.proc_nE)

        def get_proc_edgeaddr(self):          
                read_proc_edgeaddr(self.p, self.np, self.proc_edgeaddr)

        def get_proc_faceaddr(self):          
                read_proc_faceaddr(self.p, self.np, self.proc_faceaddr)

        def get_proc_edgeOwners(self):          
                read_proc_edgeOwners(self.p, self.np, self.proc_edgeOwners)

        def get_parallel_data(self, report = True):
                if self.np > 1 and ((self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True) or (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True)):
                    if report == True:
                        print('Reading parallel information for the '+ str(self.np) +' processors ...')
                    self.get_proc_faBoundary()
                    self.get_proc_nF()
                    self.get_proc_nE()
                    self.get_proc_ec()
                    self.get_proc_edgeaddr()
                    self.get_proc_faceaddr()
                    self.get_proc_edgeOwners()
                    self.order_proc_faBoundary()
                    self.order_proc_edgeOwners()

                    procBoundaries = {}
                    for proc_key in self.proc_fB.keys():
                        for patch_key in self.proc_fB[proc_key]['boundaries'].keys():
                            if patch_key[:12] == 'procBoundary':
                                proc1 = patch_key[12: patch_key.find('t')]
                                proc2 = patch_key[patch_key.find('t')+2:]
                                if proc1 < proc2:
                                    edges = []
                                    owner = []
                                    dict = {}
                                    indexes = self.proc_fB[proc_key]['boundaries'][patch_key]
                                    for i in indexes:
                                        edges.append(self.proc_edgeaddr[proc_key][i])
                                        owner.append(self.proc_faceaddr[proc_key][self.proc_edgeOwners[proc_key][i]])
                                    dict['list'] = edges
                                    dict['owner'] = owner
                                    procBoundaries[patch_key] = dict
 
                    if self.eO == {}:
                        self.get_edgeOwner(False)
                    if self.eN == {}:
                        self.get_edgeNeighbour(False)

                    n_problems = 0                       
                    for key in procBoundaries.keys():
                        list = []                        
                        for i in range(len(procBoundaries[key]['list'])):
                            if procBoundaries[key]['owner'][i] == int(self.eO['edges'][procBoundaries[key]['list'][i]]):
                                list.append(1)
                            elif procBoundaries[key]['owner'][i] == int(self.eN['edges'][procBoundaries[key]['list'][i]]):
                                list.append(2)
##I am almost sure that this problem will never happen again. Anyway I think it is convenient to keep this just in case //AG
                            else:
                                list.append(0)
                                n_problems += 1
                        self.proc_owner[key] = list
                        
                    if n_problems != 0:
                        print(str(n_problems)+' problematic edges has been detected when reading the parallel data ...')
                                        
                    for key in procBoundaries.keys():
                        self.proc_edgeList[key] = procBoundaries[key]['list']

                self.clean_proc_dicts()                                                                                    

        def order_proc_faBoundary(self):
                indexes_internal = []
                for proc_key in self.proc_fB.keys():  
                    for patch_key in self.proc_fB[proc_key]['boundaries'].keys():
                        for i in range(len(self.proc_fB[proc_key]['boundaries'][patch_key])):
                            index = self.proc_fB[proc_key]['boundaries'][patch_key][i]            
                            if self.proc_edgeaddr[proc_key][index] < self.nE['internal']:
                                if not(self.proc_edgeaddr[proc_key][index] in indexes_internal):
                                    indexes_internal.append(self.proc_edgeaddr[proc_key][index])                                    
                indexes_internal.sort()
                
                ec_dict = {}
                ec_list1 = []
                ec_list2 = []
                for i in indexes_internal:
                    ec_list1.append(i)
                    ec_list2.append(self.ec['0'].iF.field[i])
                ec_dict['internal'] = [ec_list1, ec_list2]
                ec_patch_dict = {}
                for key in self.ec['0'].bF.d.keys():
                    ec_patch_list1 = []
                    ec_patch_list2 = []
                    for i in range(len(self.ec['0'].bF.d[key]['value']['List'])):
                        ec_patch_list1.append(self.fB['boundaries'][key][i])
                        ec_patch_list2.append(self.ec['0'].bF.d[key]['value']['List'][i])
                    ec_patch_dict[key] = [ec_patch_list1, ec_patch_list2]
                ec_dict['boundary'] = ec_patch_dict

                ec_proc_dict = {}
                for proc_key in self.proc_ec['0'].keys():
                    patch_dict = {}
                    for patch_key in self.proc_ec['0'][proc_key].bF.d.keys():
                        proc_patch_list1 = []
                        proc_patch_list2 = []
                        for i in range(len(self.proc_ec['0'][proc_key].bF.d[patch_key]['value']['List'])):
                            index = self.proc_fB[proc_key]['boundaries'][patch_key][i]
                            proc_patch_list1.append(index)
                            proc_patch_list2.append(self.proc_ec['0'][proc_key].bF.d[patch_key]['value']['List'][i])
                        patch_dict[patch_key] = [proc_patch_list1, proc_patch_list2]
                    ec_proc_dict[proc_key] = patch_dict
                                            
                for key_proc in ec_proc_dict.keys():
                    for key_patch in ec_proc_dict[key_proc].keys():
                        if key_patch in ec_dict['boundary'].keys():
                            for i in range(len(ec_proc_dict[key_proc][key_patch][0])):
                                index = ec_dict['boundary'][key_patch][1].index(ec_proc_dict[key_proc][key_patch][1][i])
                                self.proc_edgeaddr[key_proc][ec_proc_dict[key_proc][key_patch][0][i]] = ec_dict['boundary'][key_patch][0][index]
                        else:
                            for i in range(len(ec_proc_dict[key_proc][key_patch][0])):
                                index = ec_dict['internal'][1].index(ec_proc_dict[key_proc][key_patch][1][i])
                                self.proc_edgeaddr[key_proc][ec_proc_dict[key_proc][key_patch][0][i]] = ec_dict['internal'][0][index]

#Here only the edgeOwners related to the patches of type 'processor' are corrected. The reason is that correcting the others is complicated
#The truth is that the edgeOwners of all the boundaries are wrong, but here only a few are used //AG                               
        def order_proc_edgeOwners(self):
                for proc_key in self.proc_fB.keys():
                    for patch_key in self.proc_fB[proc_key]['boundaries'].keys():
                        if patch_key[:12] == 'procBoundary':
                            for i in range(len(self.proc_fB[proc_key]['boundaries'][patch_key])):
                                index = self.proc_edgeaddr[proc_key][self.proc_fB[proc_key]['boundaries'][patch_key][i]]
                                owner = int(self.eO['edges'][index])
                                try:
                                    proc_owner = self.proc_faceaddr[proc_key].index(owner)
                                except:
                                    ngb = int(self.eN['edges'][index])
                                    proc_owner = self.proc_faceaddr[proc_key].index(ngb)
                                self.proc_edgeOwners[proc_key][self.proc_fB[proc_key]['boundaries'][patch_key][i]] = proc_owner
                    
        def get_h(self, report = True, time = -1):                
                if (self.h_flag == 'on' or self.h_flag == 'yes' or self.h_flag == True):
                    if report == True:
                        print('Reading h field ...')
                    if time == -1:
                        for t in self.t:
                            read_h(self.p, t, self.h,self.nF, self.fB)
                    else:
                        read_h(self.p, time, self.h,self.nF, self.fB)
                                            
        def get_pb(self, report = True, time = -1):                
                if (self.pb_flag == 'on' or self.pb_flag == 'yes' or self.pb_flag == True):
                    if report == True:
                        print('Reading pb field ...')
                    if time == -1:
                        for t in self.t:
                            read_pb(self.p, t, self.pb,self.nF, self.fB)
                    else:
                        read_pb(self.p, time, self.pb,self.nF, self.fB)

        def get_Cv(self, report = True, time = -1):               
                if (self.Cv_flag == 'on' or self.Cv_flag == 'yes' or self.Cv_flag == True):
                    if report == True:
                        print('Reading Cv field ...')
                    if time == -1:
                        for t in self.t:
                            read_Cv(self.p, t, self.Cv,self.nF, self.fB)
                    else:
                        read_Cv(self.p, time, self.Cv,self.nF, self.fB)

        def get_deltaz0(self, report = True, time = -1):
                if (self.deltaz0_flag == 'on' or self.deltaz0_flag == 'yes' or self.deltaz0_flag == True):
                    if self.tps == {}:
                        self.get_transportProperties(False)
                    if self.tps['terrainModification'] == 'on' or self.tps['terrainModification'] == True:
                        self.get_deltac0(True, time)
                    elif self.tps['entrainmentModel'] != 'entrainmentOff' and self.tps['depositionModel'] != 'depositionOff':
                        self.get_deltah0(True, time)
                        
                    self.clean_tps()

        def get_deltah0(self, report = True, time = -1):
                if (self.deltaz0_flag == 'on' or self.deltaz0_flag == 'yes' or self.deltaz0_flag == True):
                    if report == True:
                        print('Reading deltah0 field ...')
                    if time == -1:
                        for t in self.t:
                            read_deltah0(self.p, t, self.deltah0,self.nF, self.fB)
                    else:
                        read_deltah0(self.p, time, self.deltah0,self.nF, self.fB)

        def get_deltac0(self, report = True, time = -1):                
                if (self.deltaz0_flag == 'on' or self.deltaz0_flag == 'yes' or self.deltaz0_flag == True):
                    if report == True:
                        print('Reading deltac0 field ...')
                    if time == -1:
                        for t in self.t:
                            read_deltac0(self.p, t, self.deltac0,self.nF, self.fB)
                    else:
                        read_deltac0(self.p, time, self.deltac0,self.nF, self.fB)

        def get_Us(self, report = True, time = -1):                
                if (self.Us_flag == 'on' or self.Us_flag == 'yes' or self.Us_flag == True):
                    if report == True:
                        print('Reading Us field ...')
                    if time == -1:
                        for t in self.t:
                            read_Us(self.p, t, self.Us,self.nF, self.fB)
                        self.check_h_U_values(False)
                    else:
                        read_Us(self.p, time, self.Us,self.nF, self.fB)
                        self.check_h_U_values(False, time)
                                
        def get_tau(self, report = True, time = -1):                
                if (self.tau_flag == 'on' or self.tau_flag == 'yes' or self.tau_flag == True):
                    if report == True:
                        print('Reading tau field ...')
                    if time == -1:
                        for t in self.t:
                            read_tau(self.p, t, self.tau,self.nF, self.fB)
                    else:
                        read_tau(self.p, time, self.tau,self.nF, self.fB)

        def get_Q(self, report = True, time = -1):               
                if (self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True):
                    if report == True:
                        print('Reading Q field ...')
                    if time == -1:
                        for t in self.t:
                            read_Q(self.p, t, self.Q, self.nE['internal'], self.fB)
                    else:
                        read_Q(self.p, time, self.Q, self.nE['internal'], self.fB)
                    if self.np>1:
                        self.correct_parallel_Q(False, time)

        def get_phi2s(self, report = True, time = -1):                
                if (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True):
                    if report == True:
                        print('Reading phi2s field ...')
                    if time == -1:
                        for t in self.t:
                            read_phi2s(self.p, t, self.phi2s, self.nE['internal'], self.fB)
                    else:
                        read_phi2s(self.p, time, self.phi2s, self.nE['internal'], self.fB)
                    if self.np>1:
                        self.correct_parallel_phi2s(False, time)
                  
        def get_c(self, report = True, time = -1):
                if report == True:
                    print('Reading c field ...')
                if (self.c_flag == 'on' or self.c_flag == 'yes' or self.c_flag == True):
                    if time == -1:
                        for t in self.t:
                            read_c(self.p, t, self.c,self.nF, self.fB)
                    else:
                        read_c(self.p, time, self.c,self.nF, self.fB)
                else:    
                    read_c(self.p,'0',self.c,self.nF, self.fB)

        def get_n(self, report = True, time = -1):
                if report == True:
                    print('Reading n field ...')
                if (self.n_flag == 'on' or self.n_flag == 'yes' or self.n_flag == True):
                    if time == -1:
                        for t in self.t:
                            read_n(self.p, t, self.n,self.nF, self.fB)
                    else:
                        read_n(self.p, time, self.n,self.nF, self.fB)
                else:    
                    read_n(self.p,'0',self.n,self.nF, self.fB)

        def get_he(self, report = True, time = -1):
                if report == True:
                    print('Reading he field ...')             
                if (self.he_flag == 'on' or self.he_flag == 'yes' or self.he_flag == True):
                    if time == -1:
                        for t in self.t:
                            read_he(self.p, t, self.he,self.nF, self.fB)
                    else:
                        read_he(self.p, time, self.he,self.nF, self.fB)
                else:
                    read_he(self.p,'0',self.he,self.nF, self.fB)
 
        def get_A(self, report = True, time = -1):
                if report == True:
                    print('Reading A field ...')             
                if (self.A_flag == 'on' or self.A_flag == 'yes' or self.A_flag == True):
                    if time == -1:
                        for t in self.t:
                            read_A(self.p, t, self.A,self.nF, self.fB)
                    else:
                        read_A(self.p, time, self.A,self.nF, self.fB)
                else:
                    read_A(self.p,'0',self.A,self.nF, self.fB)
                                              
        def check_h_U_values(self, report = True, time = -1):
                if (self.h_flag == 'on' or self.h_flag == 'yes' or self.h_flag == True) and (self.Us_flag == 'on' or self.Us_flag == 'yes' or self.Us_flag == True):
                    if report == True:
                        print('Checking h and U values...')
                    if self.hmin == -1:
                        self.get_hmin(False)
                    if time == -1:
                        for t in self.t:
                            self.correct_h_U_values(t)
                    else:
                        self.correct_h_U_values(time)

        def correct_h_U_values(self, time):
                if (str(time) in self.h.keys()) == False:
                    self.get_h(False, time)
                h = self.h[str(time)].iF.field
                for i in range(len(h)):
                    if h[i]<= self.hmin*10:
                        self.Us[str(time)].iF.field[i] = vector(0,0,0)            
                         
        def correct_parallel_Q(self, report = True, time = -1):
                if (self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True):
                    if report == True:
                        print('Correcting Q field in parallel simulation ...')
                    if time == -1:
                        for t in self.t:
                            correct_edge_proc(self.Q[str(t)].iF, self.proc_owner, self.proc_edgeList)
                    else:
                        correct_edge_proc(self.Q[str(time)].iF, self.proc_owner, self.proc_edgeList)
                        
        def correct_parallel_phi2s(self, report = True, time = -1):
                if (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True):
                    if report == True:
                        print('Correcting phi2s field in parallel simulation ...')
                    if time == -1:               
                        for t in self.t:
                            correct_edge_proc(self.phi2s[str(t)].iF, self.proc_owner, self.proc_edgeList) 
                    else:
                        correct_edge_proc(self.phi2s[str(time)].iF, self.proc_owner, self.proc_edgeList)

        def get_transportProperties(self, report = True):
                if report == True:
                    print ('Reading transportProperties ...')
                read_transportProperties(self.p,self.tps)

        def get_densities(self, report = True):
                if report == True:
                    print('Reading rho_w y rho_s ...')
                dict = {}
                read_transportProperties(self.p,dict)
                self.rho_w = dict['rho_w']['value']
                self.rho_s = dict['rho_s']['value']
                self.rho_b = dict['rho_b']['value']

        def get_entrainment_coeffs(self, report = True):
                if report == True:
                    print('Reading entrainmentModel coefficients ...')
                dict = {}
                read_transportProperties(self.p,dict)
                self.eM = dict['entrainmentModel']
                
                if self.eM != 'entrainmentOff':
                    self.eM_coeffs = dict[self.eM+'Coeffs']
                               
        def get_controlDict(self, report = True):
                if report == True:
                    print ('Reading controlDict ...')
                dict = {}
                read_controlDict(self.p,dict)
                self.app = dict['application']
                self.sT  = dict['startTime']
                self.eT  = dict['endTime']
                self.wI  = dict['writeInterval']
                self.maxCo = dict['maxCo']
                                        
        def get_releaseArea(self, report = True):
                if report == True:
                    print ('Reading releaseArea ...')
                read_releaseArea(self.p,self.rA)

        def get_releaseFlow(self, report = True):
                if report == True:
                    print ('Reading releaseFlow ...')
                read_releaseFlow(self.p,self.rF)

        def get_meshDict(self, report = True):
                if report == True:
                    print ('Reading meshDict ...')
                read_meshDict(self.p,self.mD)
                        
        def report_Case(self, out_path):
                fi = open(out_path, 'a')

                if self.mD == {}:
                    self.get_meshDict(False)
                if self.app == '':
                    self.get_controlDict(False)
                if self.tps == {}:
                    self.get_transportProperties(False)
                if self.rA == {}:
                    self.get_releaseArea(False)                   
                if self.rF == {}:
                    self.get_releaseFlow(False)
                if self.A == {}:
                    self.get_A(False,0)
                    
                min_d = round(math.sqrt(4*min(self.A['0'].iF.field)/math.pi),3)
                max_d = round(math.sqrt(4*max(self.A['0'].iF.field)/math.pi),3)
                avg_d = round(math.sqrt(4*sum(self.A['0'].iF.field)/len(self.A['0'].iF.field)/math.pi),3)
                    
                fi.write("\n")              
                fi.write("Summary of the inputs to the simulation"+"\n")
                fi.write("\n")
                fi.write("SurfaceFile = " + self.mD['surfaceFile'][21:-1]+"\n")
                fi.write("MaxCellSize = " + self.mD['maxCellSize']+' and local refinement = '+ self.mD['localRefinement']['terrain']['additionalRefinementLevels']+"\n")
                fi.write("Maximum diameter = "+ str(max_d) + ", Minimun diameter = "+ str(min_d) + ", Average diameter = "+ str(avg_d) +"\n")
                fi.write("Longitudinal profile = "+ self.lp_name +"\n")
                fi.write("Transversal profile = "+ self.tp_name +"\n") 
                fi.write("Application = " + self.app+"\n")
                fi.write("startTime = " + self.sT+ ", endTime = " + self.eT+", writeInterval = " + self.wI+" and maximum Courant number = " + self.maxCo+"\n")
                fi.write("Number of processors = " + str(self.np)+ " and Number of elements = " + str(self.nF)+"\n")
                
                write_dictionary(fi, self.rF, 'release Flow')
                write_dictionary(fi, self.rA, 'release Area')
                write_dictionary(fi, self.tps, 'transport Properties')

                self.mD = {}
                self.rA = {}
                self.rF = {}
                
                self.app = ''
                self.sT  = ''
                self.eT  = ''
                self.wI  = ''
                
                self.clean_tps()
                
                fi.close() 

        def clean_fields(self, report = True, time = -1):
                if time == -1:
                    if report == True:
                        print('Cleaning fields ...')
                    for t in self.t:
                        self.clean_areaFields(t)
                        self.clean_edgeFields(t)
                else:
                    if report == True:
                        print('Cleaning fields for time = '+ str(time)+' ...')
                    self.clean_areaFields(time)
                    self.clean_edgeFields(time)
                    
        def clean_areaFields(self, t):
                self.h.pop(str(t), None)
                self.Cv.pop(str(t), None)
                self.pb.pop(str(t), None)
                self.deltac0.pop(str(t), None)
                self.deltah0.pop(str(t), None)
                self.Us.pop(str(t), None)
                self.tau.pop(str(t), None)
                
        def clean_edgeFields(self, t):
                self.Q.pop(str(t), None)
                self.phi2s.pop(str(t), None)
                
        def clean_n(self, t = 0):
                self.n.pop(str(t), None)
                
        def clean_c(self, t = 0):
                self.c.pop(str(t), None)

        def clean_A(self, t = 0):
                self.A.pop(str(t), None)

        def clean_he(self, t = 0):
                self.he.pop(str(t), None)
                
        def clean_edgeON(self):
                self.eO = {}
                self.eN = {}
                
        def clean_faces(self):
                self.faces = {}
                
        def clean_points(self):
                self.points = {}
                
        def clean_ec(self, t= 0):        
                self.ec.pop(str(t), None)

        def clean_lp(self):        
                self.lp = {} 

        def clean_tp(self):        
                self.tp = {} 
                
        def clean_tps(self):        
                self.tps = {} 
                
        def clean_eM(self):
                self.eM = ''
                self.eM_coeffs = {}               

        def clean_proc_dicts(self):
                self.proc_fB = {}
                self.proc_nF = {}
                self.proc_nE = {}
                self.proc_ec = {}
                self.proc_faceaddr = {}
                self.proc_edgeaddr = {}
                self.proc_edgeOwners = {}


##########################################  write functions  ##############################################################

def write_dictionary(fi, d, name):
        fi.write("\n")                                
        fi.write("Inputs in the " + name + " dictionary" +"\n")
        fi.write("\n")
        
        for key in d.keys():
            lvl = 1
            if type(d[key]) is dict:
                write_subdict(d[key], fi, key, lvl)
            elif type(d[key]) is str:
                write_str(d[key], fi, key, lvl)
            elif type(d[key]) is list:
                write_list(d[key], fi, key, lvl)
            elif is_number(d[key]) == True:
                write_float(d[key], fi, key, lvl)

def write_subdict(d, fi, name, lvl):
        fi.write("   "*lvl+str(name) + "\n")
        fi.write("   "*lvl+"{" + "\n")
                    
        for key in d.keys():
            if type(d[key]) is dict:
                write_subdict(d[key], fi, key, lvl+1)
            elif type(d[key]) is str:
                write_str(d[key], fi, key, lvl+1)
            elif key == 'dimensions':
                write_dims(d[key], fi, key, lvl+1)
            elif type(d[key]) is list:
                write_list(d[key], fi, key, lvl+1)
            elif is_number(d[key]) == True:
                write_float(d[key], fi, key, lvl+1)
                
        fi.write("   "*lvl+"}"+"\n")
          
def write_str(string, fi, key, lvl):
        fi.write("   "*lvl+str(key) + "  "+string + "\n")

def write_float(flt, fi, key, lvl):
        fi.write("   "*lvl+str(key) + "  "+str(flt) + "\n")

def write_list(l, fi, key, lvl):
        fi.write("   "*lvl+str(key)+ "\n")
        fi.write("   "*lvl+"["+ "\n")
        for i in range(len(l)):
            if type(l[i]) is dict:
                write_subdict(l[i], fi, '', lvl+1)
            elif type(l[i]) is str:
                write_str(l[i], fi, '', lvl+1)
            elif key == 'dimensions':
                write_dims(l[i], fi, '', lvl+1)
            elif type(l[i]) is list:
                if len(l[i]) == 3:
                    write_vector(l[i], fi, '', lvl+1)
                else:
                    write_list(l[i], fi, '', lvl+1)
            elif is_number(l[i]) == True:
                write_float(l[i], fi, '', lvl+1)

        fi.write("   "*lvl+"]"+ "\n")
                
def write_dims(dims, fi, key, lvl):
        dim = "  [" +str(dims[0])
        for i in range(len(dims)-1):
            dim += ", "+str(dims[i+1])
        dim += "]"
        fi.write("   "*lvl+str(key)+ dim + "\n")                

def write_vector(l, fi, key, lvl):
        v = "  [" +str(l[0])
        for i in range(len(l)-1):
            v += ", "+str(l[i+1])
        v += "]"
        fi.write("   "*lvl+v + "\n") 


##########################################  outPut class  ##############################################################

           
#This class needs as input a runCase object, which has all the functions to read inputs. Class output performs all the calculations and write the results. //AG             
class outPut:
        def __init__(self,runCase,d_x,d_y,r_x,r_y,alpha=0,n_iterations=0,alpha_field=0,n_iterations_field=0,dist =0,n_tp=100,rank=0,Sm='off',rho='off',rcg='off',M='off',Vsed ='off',V='off'):
                self.r  = runCase
                self.x  = np.arange(runCase.xi,runCase.xf+d_x,d_x)
                self.y  = np.arange(runCase.yi,runCase.yf+d_y,d_y)
                self.z  = []
                self.he = [] 
                self.nz = [] #Probably is not useful //AG
                self.n0 = [] #Probably is not useful //AG
                self.t  = runCase.t
                self.p  = runCase.p #path
                
                self.face_points = []  #for every face we have the points that form it
                self.edge_faces  = []  #for every face we have the edges that form it
                self.face_neighbour = []  #for every face we have its neighbors
                self.edges_sense = {}  #defined to deal with the sign problem with Q and phi2s for parallel simulations in the patches type processor.
                self.edges_list  = {}  #defined to deal with the sign problem with Q and phi2s for parallel simulations in the patches type processor.
                
                self.h  = []
                self.pb = []
                self.Cv = []
                self.deltah0 = []
                self.deltac0 = []
                self.Us  = []
                self.tau = []
                self.Q = {}
                self.phi2s = {}  
                
                self.z_tp  = {}
                self.he_tp = {}
                self.h_tp  = {}
                self.Us_m_tp = {}
                self.deltac0_tp = {}
                self.deltah0_tp = {}

                self.index_matrix = []
                self.alpha_matrix = []
                self.lp = {}
                self.tp = {}
                self.tp_index_matrix = {}
                self.tp_alpha_matrix = {}
                
                self.Sm  = []
                self.rho = []
                self.rcg = []
                self.M = []
                self.V = []
                self.Vsed = []
                               
                self.rF_V = []  #releaseFlow volume input
                self.rF_M = []  #releaseFlow mass input
                
                self.alpha = alpha
                self.niter = n_iterations
                self.alphafield = alpha_field
                self.niterfield = n_iterations_field
                self.n_tp = n_tp

                self.rank = rank
                self.t_size = []
                
                self.Sm_flag  = Sm
                self.rho_flag = rho
                self.rcg_flag = rcg
                self.M_flag = M
                self.V_flag = V
                self.Vsed_flag = Vsed
                                
                self.get_case(r_x,r_y, dist)

        def get_case(self,r_x,r_y,dist):
                self.get_index_matrix(r_x,r_y)
                self.get_alpha_matrix()
                self.get_longitudinal_profiles()
                self.get_transversal_profiles(dist)
                self.get_tp_index_matrix()
                self.get_tp_alpha_matrix()                              
                self.get_face_points()
                self.get_edge_faces()
                self.get_edges()
                self.get_face_neighbour()
                self.get_scalar_flux_edges()

##########################################  main functions  ##############################################################

        def order_c(self,c_x,c_y):
                list = []
                for c in self.r.c['0'].iF.field:
                    list.append(c.List())
                for i in range(len(list)):
                    list[i].append(i)
                list_x= sorted(list,key=itemgetter(0))
                list_y= sorted(list,key=itemgetter(1))
                for i in range(len(list)):
                    c_x.append([list_x[i][0],list_x[i][3]])
                    c_y.append([list_y[i][1],list_y[i][3]])
                    
        def correct_xy(self, r_x, r_y, c_x, c_y):
                break_x = False
                break_y = False
                while True:                    
                    if (self.x[-1]-r_x) >= c_x[-1][0]:
                        self.x = np.delete(self.x, -1)
                    else:
                        break_x = True
                    if (self.y[-1]-r_y) >= c_y[-1][0]:
                        self.y = np.delete(self.y, -1)                   
                    else:
                        break_y = True
                    if break_x == True and break_y == True:
                        break
                   
        def get_index_matrix(self, r_x, r_y): #the result is a matrix with the labels of the nearer faces to every (x,y) coordinate
                if self.rank == 0:
                    print('Creating index matrix ...')
                if self.r.c == {}:
                    self.r.get_c(False)                                
                c_x= []; c_y = []
                self.order_c(c_x,c_y)
                list_x = []; list_y = []
                self.correct_xy(r_x, r_y, c_x, c_y)
                get_list(list_x, self.x, r_x, c_x, len(self.r.c['0'].iF.field), 200) #list_x now has a list of labels of the faces whose c.x value is close to every x coordinate
                get_list(list_y, self.y, r_y, c_y, len(self.r.c['0'].iF.field), 200) #list_y now has a list of labels of the faces whose c.y value is close to every y coordinate
                
                for i in range(len(self.x)):
                    sublist = []
                    sublist_x = list_x[i]
                    for j in range(len(self.y)):
                        list = []
                        sublist_y = list_y[j]
                        for k in range(len(sublist_x)):
                            for value in sublist_y[k]:
                                if value in sublist_x[k]:
                                    list.append(value)
                        sublist.append(list)                               
                    self.index_matrix.append(sublist)
                                                    
        def get_alpha_matrix(self): #Generates a matrix which elements are list object, each one with an alpha<1 value.
                if self.rank == 0:
                    print('Creating alpha matrix ...')                
                for i in range(len(self.x)):
                    list = []
                    for j in range(len(self.y)):
                        sublist = []
                        l = 0
                        v1 = vector(self.x[i],self.y[j],0)
                        for index in self.index_matrix[i][j]:
                            l += (1/(distxy(v1,self.r.c['0'].iF.field[index])+10**-6))**2
                        for index in self.index_matrix[i][j]:
                            alpha = ((1/(distxy(v1,self.r.c['0'].iF.field[index])+10**-6))**2)/l
                            sublist.append(alpha)
                        list.append(sublist)  
                    self.alpha_matrix.append(list)
                                    
        def get_z_interpolation(self): #calculate the z-matrix
                if self.rank == 0:
                    print('Creating z interpolation ...')
                max_z = -math.inf
                z = np.zeros((len(self.x),len(self.y)))
                for i in range(len(self.x)):
                    for j in range(len(self.y)):                        
                        for k in range(len(self.index_matrix[i][j])):
                            index = self.index_matrix[i][j][k]
                            z[i,j] += self.r.c['0'].iF.field[index].z*self.alpha_matrix[i][j][k]
                        if z[i,j] != 0:
                            max_z = max(max_z, z[i,j])
                for i in range(len(self.x)):
                    for j in range(len(self.y)):                        
                        if z[i,j] == 0:
                            z[i,j] = max_z
                            
                self.z = np.zeros((len(self.x),len(self.y)))
                surface_smoother(z,self.z,self.alpha,self.niter) #smooth the z-matrix using the values in the neighbors N/S/E/W

        def get_he_interpolation(self, report = True):
                if report == True:
                    print('Creating he interpolation ...')
                if self.r.he == {}:
                    self.r.get_he(False)                    
                self.he = np.zeros((len(self.x),len(self.y)))
                self.get_scalarinput_interpolation(self.r.he, self.r.get_he, self.he, self.alphafield, self.niterfield, 0)
                            
        def get_nz_interpolation(self): #calculate the nz-matrix
                if self.rank == 0:
                    print('Creating nz interpolation ...')
                if self.r.n == {}:
                    self.r.get_n(False)
                self.n0 = self.r.n['0'].iF.field
                nz = np.zeros((len(self.x),len(self.y)))
                for i in range(len(self.x)):
                    for j in range(len(self.y)):                        
                        for k in range(len(self.index_matrix[i][j])):
                            index = self.index_matrix[i][j][k]
                            nz[i,j] += -self.n0[index].List()[2]*self.alpha_matrix[i][j][k]
                
                self.nz = np.zeros((len(self.x),len(self.y)))
                surface_smoother(nz,self.nz,self.alpha,self.niter) #smooth the nz-matrix using the values in the neighbors N/S/E/W

                self.r.clean_n()
                self.n0 = []

        def get_longitudinal_profiles(self):
                if self.r.lp_name != '':
                    if self.rank == 0:
                        print('Creating longitudinal profiles ...')
                    if self.r.lp == {}:
                        self.r.get_lp(False)
                    for key in self.r.lp.keys():
                        list_coords = []
                        list_dists = [0]
                        dict = {}
                        for i in range(len(self.r.lp[key]['x_list'])):
                            x = self.r.lp[key]['x_list'][i]
                            y = self.r.lp[key]['y_list'][i]
                            list_coords.append([x,y])
                            v = vector(x,y,0)
                                
                            if i > 0:
                                vi = vector(self.r.lp[key]['x_list'][i-1], self.r.lp[key]['y_list'][i-1], 0)
                                list_dists.append(distxy(v,vi)+list_dists[-1])
                                
                        dict['list_coords'] = list_coords
                        dict['list_dists'] = list_dists
                        self.lp[key] = dict
                                                       
        def get_transversal_profiles(self, dist):
                if self.r.lp_name != '':
                    if self.rank == 0:
                        print('Creating transversal profiles ...')
                    if self.r.lp == {}:
                        self.r.get_lp(False)
                    if self.r.tp_name == '':
                        for key in self.r.lp.keys():
                            x_lp = self.r.lp[key]['x_list']
                            y_lp = self.r.lp[key]['y_list']
                            list_tp = []
        
                            for i in range(len(x_lp)-2):
                                v1 = vector(x_lp[i], y_lp[i], 0)
                                v2 = vector(x_lp[i+1], y_lp[i+1], 0)
                                v3 = vector(x_lp[i+2], y_lp[i+2], 0)
                                d1 = v1-v2
                                d2 = v3-v2
                                theta = d2.angle()-d1.angle()
                                if theta < 0: theta += 2*math.pi
                                elif theta >= 2*math.pi: theta -= 2*math.pi
                                theta = theta/2+d2.angle()
        
                                vp1 = vector(math.cos(theta), math.sin(theta),0)*dist*0.5+v2
                                vp2 = vector(math.cos(theta+math.pi), math.sin(theta+math.pi),0)*dist*0.5+v2  
    
                                if vp1.x < self.x[0] or vp1.x > self.x[-1] or  vp1.y < self.y[0] or vp1.y > self.y[-1]:                           
                                    if vp1.x < self.x[0]:
                                        alpha = (self.x[0]-v2.x)/(vp1.x-v2.x)
                                        vp1 = v2+(vp1-v2)*alpha*0.9 #0.9 is included because if the final point of the tp is in the border there will be problems further in the code //AG
                                    elif vp1.x > self.x[-1]:
                                        alpha = (self.x[-1]-v2.x)/(vp1.x-v2.x)
                                        vp1 = v2+(vp1-v2)*alpha*0.9
                                    if vp1.y < self.y[0]:
                                        alpha = (self.y[0]-v2.y)/(vp1.y-v2.y)
                                        vp1 = v2+(vp1-v2)*alpha*0.9
                                    elif vp1.y > self.y[-1]:                                    
                                        alpha = (self.y[-1]-v2.y)/(vp1.y-v2.y)
                                        vp1 = v2+(vp1-v2)*alpha*0.9
                                                                        
                                    vp2 = v2*2-vp1
    
                                if vp2.x < self.x[0] or vp2.x > self.x[-1] or  vp2.y < self.y[0] or vp2.y > self.y[-1]: 
                                    if vp2.x < self.x[0]:
                                        alpha = (self.x[0]-v2.x)/(vp2.x-v2.x)
                                        vp2 = v2+(vp2-v2)*alpha*0.9
                                    elif vp2.x > self.x[-1]:
                                        alpha = (self.x[-1]-v2.x)/(vp2.x-v2.x)
                                        vp2 = v2+(vp2-v2)*alpha*0.9
                                    if vp2.y < self.y[0]:
                                        alpha = (self.y[0]-v2.y)/(vp2.y-v2.y)
                                        vp2 = v2+(vp2-v2)*alpha*0.9
                                    elif vp2.y > self.y[-1]:                                    
                                        alpha = (self.y[-1]-v2.y)/(vp2.y-v2.y)
                                        vp2 = v2+(vp2-v2)*alpha*0.9
                                        
                                    vp1 = v2*2-vp2
                                    
                                sublist_tp = []
                                for j in range(self.n_tp+1):
                                    sublist_tp.append(vp1+(vp2-vp1)*j/self.n_tp)
                                list_tp.append(sublist_tp)
                               
                            self.tp[key] = list_tp
                    else:
                        if self.r.tp == {}:
                            self.r.get_tp(False)
                        for key in self.r.lp.keys():
                            list_tp = []
                            for i in range(len(self.r.tp[key])):
                                vp1 = self.r.tp[key][i][0]
                                vp2 = self.r.tp[key][i][1]                      
                                sublist_tp = []
                                
                                for j in range(self.n_tp+1):
                                    sublist_tp.append(vp1+(vp2-vp1)*j/self.n_tp)
                                list_tp.append(sublist_tp)
                               
                            self.tp[key] = list_tp                    
                        
                        self.r.clean_tp()
                    self.r.clean_lp()

        def get_tp_index_matrix(self):
                if self.r.lp_name != '':
                    if self.rank == 0:
                        print('Creating transversal profiles index matrix ...')
                    for key in self.tp.keys():
                        tp = self.tp[key]
                        tp_index_matrix = []
                        for i in range(len(tp)):
                            sub_index = []
                            for j in range(self.n_tp+1):
                                 v = tp[i][j]
                                 [i_x, i_y] = self.get_coords(v)
                                 list = []
                                 for k in range(len(self.index_matrix[i_x][i_y])):
                                     list.append(self.index_matrix[i_x][i_y][k])
                                 for k in range(len(self.index_matrix[i_x+1][i_y])):
                                     list.append(self.index_matrix[i_x+1][i_y][k])
                                 for k in range(len(self.index_matrix[i_x+1][i_y+1])):
                                     list.append(self.index_matrix[i_x+1][i_y+1][k])
                                 for k in range(len(self.index_matrix[i_x][i_y+1])):
                                     list.append(self.index_matrix[i_x][i_y+1][k])
                                 
                                 list.sort()
                                 new_list = [list[0]]
                                 for k in range(len(list)-1):
                                     if list[k+1] != new_list[-1]:
                                         new_list.append(list[k+1])
                                 sub_index.append(new_list)
                            tp_index_matrix.append(sub_index)
                        
                        self.tp_index_matrix[key] = tp_index_matrix

        def get_tp_alpha_matrix(self): #Generates a matrix which elements are list objects, each one with a alpha<1 value.
                if self.r.lp_name != '':
                    if self.rank == 0:
                        print('Creating transversal profiles alpha matrix ...') 
                    for key in self.tp.keys():
                        tp = self.tp[key]
                        tp_alpha_matrix = []
                        for i in range(len(tp)):
                            sub_index = []
                            for j in range(self.n_tp+1):
                                list = []
                                v = tp[i][j]
                                l = 0
                                for index in self.tp_index_matrix[key][i][j]:
                                    l += (1/(distxy(v,self.r.c['0'].iF.field[index])+10**-6))**2
                                k = 0
                                while True:
                                    index = self.tp_index_matrix[key][i][j][k]
                                    alpha = ((1/(distxy(v,self.r.c['0'].iF.field[index])+10**-6))**2)/l
                                    if alpha < 0.01:
                                        pop = self.tp_index_matrix[key][i][j].pop(k)
                                    else:
                                        k += 1
                                    if k == len(self.tp_index_matrix[key][i][j]):
                                        break
                                l = 0
                                for index in self.tp_index_matrix[key][i][j]:
                                    l += (1/(distxy(v,self.r.c['0'].iF.field[index])+10**-6))**2                                
                                for index in self.tp_index_matrix[key][i][j]:
                                    alpha = ((1/(distxy(v,self.r.c['0'].iF.field[index])+10**-6))**2)/l                                    
                                    list.append(round(alpha,6))                             
                                sub_index.append(list)
                            tp_alpha_matrix.append(sub_index)
                        
                        self.tp_alpha_matrix[key] = tp_alpha_matrix   

        def get_face_points(self): #for very face we have a list with the points that form it
                if self.rank == 0:
                    print('Creating face-points list ...')
                if self.r.faces == {}:
                    self.r.get_faces(False)
                if self.r.points == {}:
                    self.r.get_points(False)
                for i in range(self.r.nF):
                    List = self.r.faces['faces'][i]
                    subList = []
                    for j in range(len(List)):
                        p = self.r.points['points'][List[j]]
                        subList.append(vector(p[0],p[1],p[2]))
                    self.face_points.append(subList)
                    
                self.r.clean_faces()
                self.r.clean_points()
                
        def get_edge_faces(self): #for very face we have a list with the edges that form it
                if self.rank == 0:
                    print('Creating edge-faces list ...')
                if self.r.eO == {}:
                    self.r.get_edgeOwner(False)
                if self.r.eN == {}:
                    self.r.get_edgeNeighbour(False)                
                for i in range(self.r.nF):
                    self.edge_faces.append([])
                list_eO = self.r.eO['edges']
                list_eN = self.r.eN['edges']
                for k in range(len(list_eO)):
                    index = int(list_eO[k])
                    self.edge_faces[index].append(k)
                for k in range(len(list_eN)):
                    index = int(list_eN[k])
                    self.edge_faces[index].append(k)
                                        
        def get_edges(self): #creates self.ec, which is a list of vectors. For every edge one vector, including the boundary edges.
                if self.rank == 0:    
                    print('Creating ordered list of edges ...')
                if self.r.ec == {}:
                    self.r.get_ec(False)
                self.order_edges(self.r.ec['0'].iF.field)
                self.r.clean_ec()
                
        def order_edges(self, ec): #order self.edge_faces, now the labels follow the order of face_points
                for i in range(len(self.edge_faces)):
                    list = []
                    for j in range(len(self.edge_faces[i])-1):
                        point = (self.face_points[i][j]+self.face_points[i][j+1])/2
                        dist = 10**6
                        for k in range(len(self.edge_faces[i])):
                            if distxy(point,ec[self.edge_faces[i][k]]) <= dist:
                                dist = distxy(point,ec[self.edge_faces[i][k]])
                                index = self.edge_faces[i][k]
                        list.append(index)
                        
                    point = (self.face_points[i][-1]+self.face_points[i][0])/2
                    dist = 10**6
                    for k in range(len(self.edge_faces[i])):
                        if distxy(point,ec[self.edge_faces[i][k]]) < dist:
                            dist = distxy(point,ec[self.edge_faces[i][k]])
                            index = self.edge_faces[i][k]
                    list.append(index)
                    self.edge_faces[i] = list
                    
        def get_face_neighbour(self): #creates a list of lists, where every list have the labels of the neighbour faces
                if self.rank == 0:    
                    print('Creating face neighbours list ...')
                if self.r.eO == {}:
                    self.r.get_edgeOwner(False)
                if self.r.eN == {}:
                    self.r.get_edgeNeighbour(False)                
                for i in range(self.r.nF):
                    self.face_neighbour.append([])                
                list_eO = self.r.eO['edges']
                list_eN = self.r.eN['edges']
                for k in range(len(list_eN)):
                    self.face_neighbour[int(list_eO[k])].append(int(list_eN[k]))
                    self.face_neighbour[int(list_eN[k])].append(int(list_eO[k]))

        def get_scalar_flux_edges(self):
                if self.r.lp_name != '':
                    if (self.r.Q_flag == 'on' or self.r.Q_flag == 'yes' or self.r.Q_flag == True) or (self.r.phi2s_flag == 'on' or self.r.phi2s_flag == 'yes' or self.r.phi2s_flag == True):
                        if self.rank == 0:
                            print('Calculating edges in flux-transversal profiles ...')            
                        points = self.get_points_tps()
                        dx = math.sqrt(4/math.pi*(self.x[-1]-self.x[0])*(self.y[-1]-self.y[0])/self.r.nF)
        
                        if self.r.eO == {}:
                            self.r.get_edgeOwner(False)
                        if self.r.eN == {}:
                            self.r.get_edgeNeighbour(False) 
                        
                        for key in points.keys():
                            edges_sense_list = []
                            edges_list = []
                            for p in range(len(points[key])):
                                list_faces = self.get_faces_tp(points[key][p][0],points[key][p][1],dx)                        
                                [list_points,list_faces] = self.get_points_tp(list_faces,points[key][p][0],points[key][p][1])
                                [list_edges,list_points] = self.get_edges_tp(list_faces,list_points)
                                edges_sense = self.get_edges_sense(list_edges,list_points, dx)
                                
                                edges_list.append(list_edges)
                                edges_sense_list.append(edges_sense)
                                
                            self.edges_sense[key] = edges_sense_list
                            self.edges_list[key] = edges_list
                            
                        self.r.clean_edgeON()                    


######################################  lp and tp plot functions  ##############################################################

        def plot_longitudinal_profiles(self):
                figure = plt.figure(figsize = (5,10))
                plt.grid(True)
                plt.title('Longitudinal profiles',fontsize = 17); plt.xlabel('y (m)'); plt.ylabel('x (m)')
                for key in self.lp.keys():               
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    plt.plot(coords[:,0],coords[:,1])

        def plot_zylongitudinal_profiles(self, entrainmentZones = False):
                figure = plt.figure(figsize = (5,10))
                plt.grid(True)
                z = self.z.transpose()
                plt.contourf(self.x,self.y,z,cmap = plt.get_cmap('hot'))
                plt.colorbar()
                plt.title('Longitudinal profiles',fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                for key in self.lp.keys():               
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    line1, = plt.plot(coords[:,0],coords[:,1], lw = 5, label = 'longitudinal profile')
                    
                if entrainmentZones == True:
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)
                    eZlist = []
                    for element in self.r.eM_coeffs['entrainmentZones']:
                        if type(element) == dict:
                            eZlist.append(element)
                
                    for i in range(len(eZlist)):
                        off = eZlist[i]['offset'][0][:-1]
                        polygon = eZlist[i]['vertices']
                        rectangle = np.zeros((len(polygon)+1,2))
                        
                        for j in range(len(polygon)):
                            rectangle[j,0] = polygon[j][0]-off[0]
                            rectangle[j,1] = polygon[j][1]-off[1]
                            
                        rectangle[-1,0] = rectangle[0,0]
                        rectangle[-1,1] = rectangle[0,1] 
                            
                        line2, = plt.plot(rectangle[:,0], rectangle[:,1], lw = 5, color='black', label = 'entrainment zone '+ str(i+1))
                        
                    self.r.clean_eM()
                                                
                plt.legend(handler_map= {}) 

        def plot_hylongitudinal_profiles(self, t, entrainmentZones = False):
                figure = plt.figure(figsize = (5,10))
                plt.grid(True)
                j = self.t.index(t)
                field_t = self.h[j,:,:].transpose()
                plt.contourf(self.x,self.y,field_t,cmap = plt.get_cmap('hot'))
                plt.colorbar()
                plt.title('Longitudinal profiles',fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                for key in self.lp.keys():               
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    line1, = plt.plot(coords[:,0],coords[:,1], lw = 5, label = 'longitudinal profile')
                    
                if entrainmentZones == True:
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)
                    eZlist = []
                    for element in self.r.eM_coeffs['entrainmentZones']:
                        if type(element) == dict:
                            eZlist.append(element)
                
                    for i in range(len(eZlist)):
                        off = eZlist[i]['offset'][0][:-1]
                        polygon = eZlist[i]['vertices']
                        rectangle = np.zeros((len(polygon)+1,2))
                        
                        for j in range(len(polygon)):
                            rectangle[j,0] = polygon[j][0]-off[0]
                            rectangle[j,1] = polygon[j][1]-off[1]
                            
                        rectangle[-1,0] = rectangle[0,0]
                        rectangle[-1,1] = rectangle[0,1] 
                            
                        line2, = plt.plot(rectangle[:,0], rectangle[:,1], lw = 5, color='black', label = 'entrainment zone '+ str(i+1))
                        
                    self.r.clean_eM()
                                                
                plt.legend(handler_map= {}) 

        def plot_hytransversal_profiles(self, t, entrainmentZones = False):
                figure = plt.figure(figsize = (5,10))
                plt.grid(True)
                j = self.t.index(t)
                field_t = self.h[j,:,:].transpose()
                plt.contourf(self.x,self.y,field_t,cmap = plt.get_cmap('hot'))
                plt.colorbar()
                plt.title('Transversal Profiles',fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                for key in self.lp.keys():               
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    line1, = plt.plot(coords[:,0],coords[:,1], lw = 5, label = 'longitudinal profile')

                    for i in range(len(self.tp[key])):
                        vertexes = self.tp[key][i]
                        plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])
                    
                if entrainmentZones == True:
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)
                    eZlist = []
                    for element in self.r.eM_coeffs['entrainmentZones']:
                        if type(element) == dict:
                            eZlist.append(element)
                
                    for i in range(len(eZlist)):
                        off = eZlist[i]['offset'][0][:-1]
                        polygon = eZlist[i]['vertices']
                        rectangle = np.zeros((len(polygon)+1,2))
                        
                        for j in range(len(polygon)):
                            rectangle[j,0] = polygon[j][0]-off[0]
                            rectangle[j,1] = polygon[j][1]-off[1]
                            
                        rectangle[-1,0] = rectangle[0,0]
                        rectangle[-1,1] = rectangle[0,1] 
                            
                        line2, = plt.plot(rectangle[:,0], rectangle[:,1], lw = 5, color='black', label = 'entrainment zone '+ str(i+1))
                        
                    self.r.clean_eM()
                                                
                plt.legend(handler_map= {}) 
                                       
        def plot_transversal_profiles(self):
                figure = plt.figure(figsize = (8,20))
                plt.grid(True)
                plt.title('Transversal  Profiles',fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                for key in self.lp.keys():
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    plt.plot(coords[:,0],coords[:,1])
                    
                    for i in range(len(self.tp[key])):
                        vertexes = self.tp[key][i]
                        plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])

        def plot_zytransversal_profiles(self):
                figure = plt.figure(figsize = (8,20))
                plt.grid(True)
                z = self.z.transpose()
                plt.contourf(self.x,self.y,z,cmap = plt.get_cmap('hot'))
                plt.colorbar()
                plt.title('Transversal  Profiles',fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                for key in self.lp.keys():
                    coords = np.zeros((len(self.lp[key]['list_coords']),2))
                    for i in range(len(coords)):
                        coords[i,0] = self.lp[key]['list_coords'][i][0]
                        coords[i,1] = self.lp[key]['list_coords'][i][1]
                    plt.plot(coords[:,0],coords[:,1], lw = 5)
                    
                    for i in range(len(self.tp[key])):
                        vertexes = self.tp[key][i]
                        plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])

        def plot_transversal_profile(self, key, list, limit = True, entrainmentZones = False):
                figure = plt.figure(figsize = (4,10))
                plt.grid(True)

                coords = np.zeros((len(self.lp[key]['list_coords']),2))
                for i in range(len(coords)):
                    coords[i,0] = self.lp[key]['list_coords'][i][0]
                    coords[i,1] = self.lp[key]['list_coords'][i][1]
                line1, = plt.plot(coords[:,0],coords[:,1], lw = 5, label = 'longitudinal profile')
                if limit == True:
                    plt.xlim([self.x[0], self.x[-1]])
                    plt.ylim([self.y[0], self.y[-1]])

                names = ''                
                for i in range(len(list)):
                    names += str(list[i])+' '
                    vertexes = self.tp[key][list[i]]
                    plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y])

                plt.title('Transversal  Profile n = '+ names,fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')

                if entrainmentZones == True:
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)
                    eZlist = []
                    for element in self.r.eM_coeffs['entrainmentZones']:
                        if type(element) == dict:
                            eZlist.append(element)
                
                    for i in range(len(eZlist)):
                        off = eZlist[i]['offset'][0][:-1]
                        polygon = eZlist[i]['vertices']
                        rectangle = np.zeros((len(polygon)+1,2))
                        
                        for j in range(len(polygon)):
                            rectangle[j,0] = polygon[j][0]-off[0]
                            rectangle[j,1] = polygon[j][1]-off[1]
                            
                        rectangle[-1,0] = rectangle[0,0]
                        rectangle[-1,1] = rectangle[0,1] 
                            
                        line2, = plt.plot(rectangle[:,0], rectangle[:,1], lw = 5, color='black', label = 'entrainment zone '+ str(i+1))
                        
                    self.r.clean_eM()
                                                
                plt.legend(handler_map= {})

        def plot_zytransversal_profile(self, key, list, entrainmentZones = False, save = False):
                figure = plt.figure(figsize = (6,12))
                plt.grid(True)
                z = self.z.transpose()
                plt.contourf(self.x,self.y,z,cmap = plt.get_cmap('hot'))
                plt.colorbar()

                coords = np.zeros((len(self.lp[key]['list_coords']),2))
                for i in range(len(coords)):
                    coords[i,0] = self.lp[key]['list_coords'][i][0]
                    coords[i,1] = self.lp[key]['list_coords'][i][1]
                line1, = plt.plot(coords[:,0],coords[:,1], lw = 5, label = 'longitudinal profile')

                for i in range(len(list)):
                    vertexes = self.tp[key][list[i]]
                    plt.plot([vertexes[0].x, vertexes[-1].x] ,[vertexes[0].y, vertexes[-1].y], lw = 3)

                if entrainmentZones == True:
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)
                    eZlist = []
                    for element in self.r.eM_coeffs['entrainmentZones']:
                        if type(element) == dict:
                            eZlist.append(element)
                
                    for i in range(len(eZlist)):
                        off = eZlist[i]['offset'][0][:-1]
                        polygon = eZlist[i]['vertices']
                        rectangle = np.zeros((len(polygon)+1,2))
                        
                        for j in range(len(polygon)):
                            rectangle[j,0] = polygon[j][0]-off[0]
                            rectangle[j,1] = polygon[j][1]-off[1]
                            
                        rectangle[-1,0] = rectangle[0,0]
                        rectangle[-1,1] = rectangle[0,1] 
                            
                        line2, = plt.plot(rectangle[:,0], rectangle[:,1], lw = 5, color='black', label = 'entrainment zone '+str(i+1))
                        
                    self.r.clean_eM()
                                                
                plt.legend(handler_map= {})
                if save == True:
                    path = self.p+'/images'
                    figure.savefig(path+'/'+'zplot')
                    plt.close(figure)


##########################################  tp interpolation functions  ##############################################################                          
                                                                  
        def get_scalarinput_tpinterpolation(self, field, output_field, key, time, ix = -1):
                if time == -1:
                    for t in range(len(self.t)):
                        self.get_scalar_field_tpinterpolation(field[str(self.t[t])],output_field[t,:,:], key, ix)
                else:
                    self.get_scalar_field_tpinterpolation(field[str(time)],output_field, key, ix)

        def get_vectorinput_tpinterpolation(self, field, output_field, key, time):
                if time == -1:
                    for t in range(len(self.t)):
                        vector_field = np.zeros((len(self.tp[key]),self.n_tp+1,3))
                        self.get_scalar_field_tpinterpolation(field[str(self.t[t])],vector_field[:,:,0], key, 0)
                        self.get_scalar_field_tpinterpolation(field[str(self.t[t])],vector_field[:,:,1], key, 1)
                        self.get_scalar_field_tpinterpolation(field[str(self.t[t])],vector_field[:,:,2], key, 2)
                        
                        for i in range(len(self.tp[key])):
                            for j in range(self.n_tp+1):
                                output_field[t,i,j] = round(vector(vector_field[i,j,0],vector_field[i,j,1],vector_field[i,j,2]).mag(),3)
                else:
                    vector_field = np.zeros((len(self.tp[key]),self.n_tp+1,3))
                    self.get_scalar_field_tpinterpolation(field[str(time)],vector_field[:,:,0], key, 0)
                    self.get_scalar_field_tpinterpolation(field[str(time)],vector_field[:,:,1], key, 1)
                    self.get_scalar_field_tpinterpolation(field[str(time)],vector_field[:,:,2], key, 2)                    

                    for i in range(len(self.tp[key])):
                        for j in range(self.n_tp+1):
                            output_field[i,j] = round(vector(vector_field[i,j,0],vector_field[i,j,1],vector_field[i,j,2]).mag(),3)
                    
        def get_scalar_field_tpinterpolation(self, field, output_field, key, ix = -1):
                for i in range(output_field.shape[0]):
                    for j in range(output_field.shape[1]):
                        value = 0
                        if ix == -1:
                            if field.iF.u == 'nonuniform':
                                for k in range(len(self.tp_index_matrix[key][i][j])):
                                    index = self.tp_index_matrix[key][i][j][k]
                                    alpha = self.tp_alpha_matrix[key][i][j][k]                                
                                    value += field.iF.field[index]*alpha
                            else:
                                value = field.iF.field[0]
                        else:
                            if field.iF.u == 'nonuniform':
                                for k in range(len(self.tp_index_matrix[key][i][j])):
                                    index = self.tp_index_matrix[key][i][j][k]
                                    alpha = self.tp_alpha_matrix[key][i][j][k]
                                    value += field.iF.field[index].List()[ix]*alpha
                            else:
                                value = field.iF.field[0].List()[ix]
                        output_field[i,j] = round(value,6)

        def get_z_tpinterpolation(self, report = True):
                if self.r.lp_name != '':
                    if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                        if report == True:
                            print('Interpolating z field in transversal profiles ...')
                        for key in self.tp.keys():
                            z_tp = np.zeros((len(self.tp[key]),self.n_tp+1)) 
                            self.get_scalarinput_tpinterpolation(self.r.c, z_tp, key, 0, 2)
                            self.z_tp[key] = curve_smoother(z_tp,self.alpha,self.niter)
                        
        def get_he_tpinterpolation(self, report = True):
                if self.r.lp_name != '':
                    if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                        if report == True:
                            print('Interpolating he field in transversal profiles ...')
                        if self.r.he == {}:
                            self.r.get_he(False)                          
                        for key in self.tp.keys():
                            he_tp = np.zeros((len(self.tp[key]),self.n_tp+1)) 
                            self.get_scalarinput_tpinterpolation(self.r.he, he_tp, key, 0)
                            self.he_tp[key] = he_tp
                            
                        self.r.clean_he()

        def get_h_tpinterpolation(self, report = True, time = -1):
                if self.r.lp_name != '':
                    if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                        if report == True:
                            print('Interpolating h field in transversal profiles ...')
                        for key in self.tp.keys():
                            if time == -1:
                                h_tp = np.zeros((len(self.t),len(self.tp[key]),self.n_tp+1))
                            else:
                                h_tp = np.zeros((len(self.tp[key]),self.n_tp+1))
                            self.get_scalarinput_tpinterpolation(self.r.h, h_tp, key, time)
                            self.h_tp[key] = h_tp
    
        def get_deltaz0_tpinterpolation(self, report = True, time = -1):
                if self.r.lp_name != '':
                    if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                        if self.r.tps == {}:
                            self.r.get_transportProperties(False)
                        if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                            self.get_deltac0_tpinterpolation(True, time)
                        elif self.r.tps['entrainmentModel'] != 'entrainmentOff' and self.r.tps['depositionModel'] != 'depositionOff':
                            self.get_deltah0_tpinterpolation(True, time)
                            
                        self.r.clean_tps()

        def get_deltah0_tpinterpolation(self, report = True, time = -1):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if report == True:
                        print('Interpolating deltah0 field in transversal profiles ...')
                    for key in self.tp.keys():
                        if time == -1:
                            deltah0_tp = np.zeros((len(self.t),len(self.tp[key]),self.n_tp+1))
                        else:
                            deltah0_tp = np.zeros((len(self.tp[key]),self.n_tp+1))
                        self.get_scalarinput_tpinterpolation(self.r.deltah0, deltah0_tp, key, time)
                        self.deltah0_tp[key] = deltah0_tp

        def get_deltac0_tpinterpolation(self, report = True, time = -1):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if report == True:
                        print('Interpolating deltac0 field in transversal profiles...')
                    for key in self.tp.keys():
                        if time == -1:
                            deltac0_tp = np.zeros((len(self.t),len(self.tp[key]),self.n_tp+1))
                        else:
                            deltac0_tp = np.zeros((len(self.tp[key]),self.n_tp+1))
                        self.get_scalarinput_tpinterpolation(self.r.deltac0, deltac0_tp, key, time)
                        self.deltac0_tp[key] = deltac0_tp

        def get_Us_m_tpinterpolation(self, report = True, time = -1):
                if self.r.lp_name != '':
                    if (self.r.Us_flag == 'on' or self.r.Us_flag == 'yes' or self.r.Us_flag == True):
                        if report == True:
                            print('Interpolating Us field in transversal profiles ...')
                        for key in self.tp.keys():
                            if time == -1:
                                Us_m_tp = np.zeros((len(self.t),len(self.tp[key]),self.n_tp+1))
                            else:
                                Us_m_tp = np.zeros((len(self.tp[key]),self.n_tp+1))
                            self.get_vectorinput_tpinterpolation(self.r.Us, Us_m_tp, key, time)
                            self.Us_m_tp[key] = Us_m_tp


##########################################  tp plot functions  ##############################################################                          

        def plot_h_tp(self, key, t, index, xmin = -1, xmax = -1):
                j = self.t.index(t)
                h_tp = self.h_tp[key][j,index,:]
                he_tp = self.he_tp[key][index,:]
                z_tp = self.z_tp[key][index,:]
                
                dists = [0]*len(self.tp[key][index])
                dist_max = distxy(self.tp[key][index][0], self.tp[key][index][-1])
            
                for i in range(len(dists)):
                    dists[i] = distxy(self.tp[key][index][0], self.tp[key][index][i])-dist_max/2

                if xmin != -1:
                    imin = int((xmin-dists[0])/dist_max*len(dists))
                else:
                    imin = 0
                if xmax != -1:
                    imax = int((xmax-dists[0])/dist_max*len(dists))+1
                else:
                    imax = len(dists)
                    
                figure = plt.figure(figsize = (12,5))
                plt.grid(True)
                tit = 't = '+ str(t) +'s'
                plt.title(tit,fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('z (m)')
                
                line1, = plt.plot(dists[imin:imax], z_tp[imin:imax],'orange', label = 'initial bottom')
                line2, = plt.plot(dists[imin:imax], (z_tp-he_tp)[imin:imax], 'y', label = 'initial erodible layer')

                if self.deltac0_tp != {}:
                    deltac0_tp = self.deltac0_tp[key][j,index,:]
                    line3, = plt.plot(dists[imin:imax], (h_tp+z_tp+deltac0_tp)[imin:imax], 'b', label = 'flow height')     
                    line4, = plt.plot(dists[imin:imax], (deltac0_tp+z_tp)[imin:imax], 'r', label = 'actual bottom')
                elif self.deltah0_tp != {}:
                    deltah0_tp = self.deltah0_tp[key][j,index,:]
                    line3, = plt.plot(dists[imin:imax], (h_tp+z_tp)[imin:imax], 'b', label = 'flow height')     
                    line4, = plt.plot(dists[imin:imax], (deltah0_tp+z_tp)[imin:imax], 'r', label = 'actual erodible layer')
                else:
                    line3, = plt.plot(dists[imin:imax], (h_tp+z_tp)[imin:imax], 'b', label = 'flow height') 
                                                                
                plt.legend(handler_map= {})                  

        def plot_h_tp_video(self, key, index, clean = True, xmin = -1, xmax = -1):
                j1='imagen00'
                j2='imagen0'
                j3='imagen'
                
                make_dir(self.p, 'images')
                path = self.p+'/images'
                dir = 'index'+str(index)
                clean_dir(path, dir)
                make_dir(path, dir)
                z_tp = self.z_tp[key][index,:]
                he_tp = self.he_tp[key][index,:]

                dists = [0]*len(self.tp[key][index])
                dist_max = distxy(self.tp[key][index][0], self.tp[key][index][-1])
            
                for i in range(len(dists)):
                    dists[i] = distxy(self.tp[key][index][0], self.tp[key][index][i])-dist_max/2

                if xmin != -1:
                    imin = int((xmin-dists[0])/dist_max*len(dists))
                else:
                    imin = 0
                if xmax != -1:
                    imax = int((xmax-dists[0])/dist_max*len(dists))+1
                else:
                    imax = len(dists)

                if self.deltac0_tp != {}:
                    deltac0_tp = self.deltac0_tp[key][-1,index,:]
                else:
                    deltac0_tp = z_tp*0
                   
                max_z = max((z_tp+deltac0_tp)[imin:imax])+0.2
                min_z = min(min((z_tp+deltac0_tp)[imin:imax]),min((z_tp-he_tp)[imin:imax]))-0.2
                                   
                for t in range(len(self.t)):
                    figure = plt.figure(figsize = (12,5))
                    plt.grid(True)
                    tit = 't = '+ str(self.t[t]) +'s'
                    plt.title(tit,fontsize = 17)
                    plt.xlabel('x (m)'); plt.ylabel('z (m)')
                    plt.ylim([min_z, max_z])
                    
                    h_tp = self.h_tp[key][t,index,:]
                    
                    line1, = plt.plot(dists[imin:imax], z_tp[imin:imax], 'orange', label = 'initial bottom') 
                    line2, = plt.plot(dists[imin:imax], (z_tp-he_tp)[imin:imax], 'y', label = 'initial erodible layer')

                    if self.deltac0_tp != {}:
                        deltac0_tp = self.deltac0_tp[key][t,index,:]
                        line3, = plt.plot(dists[imin:imax], (h_tp+z_tp+deltac0_tp)[imin:imax], 'b', label = 'flow height')     
                        line4, = plt.plot(dists[imin:imax], (deltac0_tp+z_tp)[imin:imax], 'r', label = 'actual bottom')
                    elif self.deltah0_tp != {}:
                        deltah0_tp = self.deltah0_tp[key][t,index,:]
                        line3, = plt.plot(dists[imin:imax], (h_tp+z_tp)[imin:imax], 'b', label = 'flow height')     
                        line4, = plt.plot(dists[imin:imax], (deltah0_tp+z_tp)[imin:imax], 'r', label = 'actual erodible layer')
                    else:
                        line3, = plt.plot(dists[imin:imax], (h_tp+z_tp)[imin:imax], 'b', label = 'flow height') 

                    plt.legend(handler_map= {})
                    
                    if t>=10 and t<100:
                        d=j2+str(t)
                    elif t<10:
                        d=j1+str(t)
                    else:
                        d=j3+str(t)
                        
                    figure.savefig(path+'/'+dir+'/'+d)
                    plt.close(figure)

                os.chdir(self.p+'/images/'+dir)
                os.system("ffmpeg -framerate 3 -i imagen%03d.png {}.mp4".format('h_tp_'+str(index)))

                if clean == True:
                    shutil.copy(self.p+'/images/'+dir+'/'+'h_tp_'+str(index)+'.mp4', self.p+'/images/'+'h_tp_'+str(index)+'.mp4')
                    clean_dir(path, dir)

        def plot_h_tp_time(self, key, index, save = False):
                self.plot_fieldtp_time(self.h_tp, key, index, 'h', 'm', 1, save)

        def plot_Us_m_tp_time(self, key, index, save = False):
                self.plot_fieldtp_time(self.Us_m_tp, key, index, 'Us_m', 'm/s', 1, save)
                
        def plot_deltac0_tp_time(self, key, index, sign, save = False):
                self.plot_fieldtp_time(self.deltac0_tp, key, index, 'deltac0', 'm', sign, save)

        def plot_deltah0_tp_time(self, key, index, sign, save = False):
                self.plot_fieldtp_time(self.deltah0_tp, key, index, 'deltah0', 'm', sign, save)

        def plot_fieldtp_time(self, field, key, index, name, units, sign, save):
                field_tp_time = np.zeros((len(self.t),1))
                for t in range(len(self.t)):
                    field_tp_time[t,0] = max(field[key][t,index,:]*sign)*sign
                    
                fig = plt.figure(figsize = (8,8))
                plt.plot(self.t,field_tp_time[:,0])
                plt.grid(True)
                tit = 'index = '+ str(index)
                plt.title(tit,fontsize = 17); plt.xlabel('t (s)'); plt.ylabel(name +' ('+ units+')')                

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+name+'_tp_time_index_'+str(index))
                    plt.close(fig)


##########################################  field interpolation functions  ##########################################################                          

        def get_h_interpolation(self, report = True, time = -1):
                if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                    if report == True:
                        print('Interpolating h field ...')
                    if time == -1:
                        self.h = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.h = np.zeros((len(self.x),len(self.y)))
                    self.get_scalarinput_interpolation(self.r.h, self.r.get_h, self.h, self.alphafield, self.niterfield, time)

        def get_Cv_interpolation(self, report = True, time = -1):
                if (self.r.Cv_flag == 'on' or self.r.Cv_flag == 'yes' or self.r.Cv_flag == True):
                    if report == True:
                        print('Interpolating Cv field ...')
                    if time == -1:
                        self.Cv = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.Cv = np.zeros((len(self.x),len(self.y)))
                    self.get_scalarinput_interpolation(self.r.Cv, self.r.get_Cv, self.Cv, self.alphafield, self.niterfield, time)

        def get_deltaz0_interpolation(self, report = True, time = -1):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if self.r.tps == {}:
                        self.r.get_transportProperties(False)
                    if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                        self.get_deltac0_interpolation(True, time)
                    elif self.r.tps['entrainmentModel'] != 'entrainmentOff' and self.r.tps['depositionModel'] != 'depositionOff':
                        self.get_deltah0_interpolation(True, time)
                        
                    self.r.clean_tps()
                               
        def get_deltah0_interpolation(self, report = True, time = -1):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if report == True:
                        print('Interpolating deltah0 field ...')
                    if time == -1:
                        self.deltah0 = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.deltah0 = np.zeros((len(self.x),len(self.y)))
                    self.get_scalarinput_interpolation(self.r.deltah0, self.r.get_deltah0, self.deltah0, self.alphafield, self.niterfield, time)

        def get_deltac0_interpolation(self, report = True, time = -1):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if report == True:
                        print('Interpolating deltac0 field ...')
                    if time == -1:
                        self.deltac0 = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.deltac0 = np.zeros((len(self.x),len(self.y)))
                    self.get_scalarinput_interpolation(self.r.deltac0, self.r.get_deltac0, self.deltac0, self.alphafield, self.niterfield, time)

        def get_pb_interpolation(self, report = True, time = -1):
                if (self.r.pb_flag == 'on' or self.r.pb_flag == 'yes' or self.r.pb_flag == True):
                    if report == True:
                        print('Interpolating pb field ...')
                    if time == -1:
                        self.pb = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.pb = np.zeros((len(self.x),len(self.y)))
                    self.get_scalarinput_interpolation(self.r.pb, self.r.get_pb, self.pb, self.alphafield, self.niterfield, time)
                                                                       
        def get_Us_interpolation(self, report = True, time = -1):
                if (self.r.Us_flag == 'on' or self.r.Us_flag == 'yes' or self.r.Us_flag == True):
                    if report == True:
                        print('Interpolating Us field ...')
                    if time == -1:
                        self.Us = np.zeros((len(self.t),len(self.x),len(self.y),3))
                        self.Us_m = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.Us = np.zeros((len(self.x),len(self.y),3))
                        self.Us_m = np.zeros((len(self.x),len(self.y)))
                    self.get_vectorinput_interpolation(self.r.Us, self.r.get_Us, self.Us, self.Us_m, self.alphafield, self.niterfield, time)   
                                                                                                                
        def get_tau_interpolation(self, report = True, time = -1):
                if (self.r.tau_flag == 'on' or self.r.tau_flag == 'yes' or self.r.tau_flag == True):
                    if report == True:
                        print('Interpolating tau field...')
                    if time == -1:
                        self.tau = np.zeros((len(self.t),len(self.x),len(self.y),3))
                        self.tau_m = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.tau = np.zeros((len(self.x),len(self.y),3))
                        self.tau_m = np.zeros((len(self.x),len(self.y)))
                    self.get_vectorinput_interpolation(self.r.tau, self.r.get_tau, self.tau, self.tau_m, self.alphafield, self.niterfield, time) 

        def get_scalarinput_interpolation(self, runCase_field, runCase_func, outPut_field, alpha, n_iterations, time):
                if time == -1:
                    if runCase_field == {}:
                        runCase_func(False) 
                    for t in range(len(self.t)):
                        self.get_scalar_field_interpolation(runCase_field,outPut_field[t,:,:],t,alpha,n_iterations)
                else:
                    if (str(time) in runCase_field.keys()) == False:
                        runCase_func(False, time)
                    t = self.t.index(time)
                    self.get_scalar_field_interpolation(runCase_field,outPut_field,t,alpha,n_iterations)

        def get_vectorinput_interpolation(self, runCase_field, runCase_func, outPut_field, outPut_field_m, alpha, n_iterations, time):
                if time == -1:
                    if runCase_field == {}:
                        runCase_func(False) 
                    for t in range(len(self.t)):
                        self.get_vector_field_interpolation(runCase_field, outPut_field[t,:,:,:],t,alpha,n_iterations)

                    for t in range(len(self.t)):
                        for i in range(len(self.x)):
                            for j in range(len(self.y)):
                                outPut_field_m[t,i,j] = round(math.sqrt(outPut_field[t,i,j,0]**2+outPut_field[t,i,j,1]**2+outPut_field[t,i,j,2]**2),6) 
                else:
                    if (str(time) in runCase_field.keys()) == False:
                        runCase_func(False, time)
                    t = self.t.index(time)
                    self.get_vector_field_interpolation(runCase_field,outPut_field,t,alpha,n_iterations)
                    
                    for i in range(len(self.x)):
                        for j in range(len(self.y)):
                            outPut_field_m[i,j] = round(math.sqrt(outPut_field[i,j,0]**2+outPut_field[i,j,1]**2+outPut_field[i,j,2]**2),6) 
                                                                                                      
        def get_scalar_field_interpolation(self,runCase_field,output_field,t,alpha,n_iterations):
                field = runCase_field[str(self.t[t])].iF.field
                field_u = runCase_field[str(self.t[t])].iF.u
                scalar_field = np.zeros((len(self.x),len(self.y)))
                if field_u == 'uniform':
                    scalar_field = np.ones((len(self.x),len(self.y)))*field[0]
                else:
                    for i in range(len(self.x)):
                        for j in range(len(self.y)):                        
                            for k in range(len(self.index_matrix[i][j])):
                                scalar_field[i,j] += field[self.index_matrix[i][j][k]]*self.alpha_matrix[i][j][k]

                surface_smoother(scalar_field,output_field,alpha,n_iterations)
      
        def get_vector_field_interpolation(self,runCase_field,output_field,t,alpha,n_iterations):
                field = runCase_field[str(self.t[t])].iF.field
                field_u = runCase_field[str(self.t[t])].iF.u
                vector_field = np.zeros((len(self.x),len(self.y),3))
                if field_u == 'uniform':
                    vector_field[:,:,0] = np.ones((len(self.x),len(self.y)))*field[0].x
                    vector_field[:,:,1] = np.ones((len(self.x),len(self.y)))*field[0].y
                    vector_field[:,:,2] = np.ones((len(self.x),len(self.y)))*field[0].z
                else:
                    for i in range(len(self.x)):
                        for j in range(len(self.y)):                        
                            for k in range(len(self.index_matrix[i][j])):
                                vector_field[i,j,0] += field[self.index_matrix[i][j][k]].x*self.alpha_matrix[i][j][k]
                                vector_field[i,j,1] += field[self.index_matrix[i][j][k]].y*self.alpha_matrix[i][j][k]
                                vector_field[i,j,2] += field[self.index_matrix[i][j][k]].z*self.alpha_matrix[i][j][k]

                surface_smoother(vector_field[:,:,0],output_field[:,:,0],alpha,n_iterations)
                surface_smoother(vector_field[:,:,1],output_field[:,:,1],alpha,n_iterations)
                surface_smoother(vector_field[:,:,2],output_field[:,:,2],alpha,n_iterations)


##########################################  field plot functions  ##########################################################

        def plot_z(self):
                self.plot_field(-1, self.z)

        def plot_he(self):
                self.plot_field(-1, self.he)
               
        def plot_nz(self):
                self.plot_field(-1, self.nz)

        def plot_h_image(self,t,im_name,off_x = 0,off_y = 0,n_intervals = 8, save = False):
                self.plot_field_image(self.h,t,im_name,'h', 'm',off_x,off_y,n_intervals,save)

        def plot_Us_image(self,t,im_name,off_x = 0,off_y = 0,n_intervals = 8, save = False):
                self.plot_field_image(self.Us_m,t,im_name,'Us', 'm/s',off_x,off_y,n_intervals,save)

        def plot_Cv_image(self,t,im_name,off_x = 0,off_y = 0,n_intervals = 8, save = False):
                self.plot_field_image(self.Cv,t,im_name,'Cv', '-',off_x,off_y,n_intervals,save)

        def plot_deltac0_image(self,t,im_name,off_x = 0,off_y = 0,n_intervals = 8, save = False):
                self.plot_field_image(self.deltac0,t,im_name,'deltac0', 'm',off_x,off_y,n_intervals,save)

        def plot_field_image(self,field,t,im_name, name, units, off_x, off_y, n_intervals, save):
                im = plt.imread(self.p+'/'+im_name)
                fig, ax = plt.subplots(figsize=(14, 6))
                extent = [self.x[0]+off_x,self.x[-1]+off_x,self.y[0]+off_y,self.y[-1]+off_y]
                ax.imshow(im,extent = extent)
                plt.grid()
                                
                x, y, arrow_length = 0.9, 0.2, 0.1
                ax.annotate('N', xy=(x, y), xytext=(x, y-arrow_length),
                    arrowprops=dict(facecolor='black', width=5, headwidth=15),
                    ha='center', va='center', fontsize=20,
                    xycoords=ax.transAxes)

                colors = []
                for c in np.linspace(0,1,n_intervals+1):
                    if len(colors) > 1:
                        colors.append((c,abs(0.5-c),1-c,1))
                    elif len(colors) == 1:
                        colors.append((c,abs(0.5-c),1-c,0.5))
                    else:
                        colors.append((c,abs(0.5-c),1-c,0))
                mycmap = mcolors.ListedColormap(colors)
                
                j = self.t.index(t)
                field_t = field[j,:,:].transpose()
                field_inv = field_t*0
                nrows = field_t.shape[0]
                for i in range(nrows):
                    field_inv[i,:] = field_t[nrows-1-i,:]

                max_field = round(np.amax(field_inv),1)
                a = [0,0.001]
                ticks = [0]
                for i in range(n_intervals):
                    a.append(round(max_field*(i+1)/(n_intervals),2))
                    ticks.append(round(max_field*(i+1)/(n_intervals),2))
                    
                a = np.array([a])
                norm = mcolors.BoundaryNorm(a[0], len(a[0])-1)
                img= plt.imshow(field_inv,cmap=mycmap,extent=extent,norm = norm)
                cbar = fig.colorbar(img, ticks = ticks, spacing = 'proportional',label= name+' ('+units+')')
                ticks_names = []
                for tick in ticks: ticks_names.append(str(tick))
                cbar.ax.set_yticklabels(ticks)
                plt.xlabel('x (m)'); plt.ylabel('y (m)')

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+name+'_t_'+str(t))
                    plt.close(fig)
                
        def plot_field(self,t,field, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                figure = plt.figure(figsize = (5,10))
                plt.grid(True)
                if t != -1:
                    j = self.t.index(t)
                    field_t = field[j,:,:].transpose()
                    if clip == False:
                        x = self.x
                        y = self.y
                    elif clip == True:
                        imin = int((xmin-self.x[0])/(self.x[1]-self.x[0]))
                        imax = int((xmax-self.x[0])/(self.x[1]-self.x[0]))+1
                        jmin = int((ymin-self.y[0])/(self.y[1]-self.y[0]))
                        jmax = int((ymax-self.y[0])/(self.y[1]-self.y[0]))+1
                        x = self.x[imin:imax]
                        y = self.y[jmin:jmax]
                        field_t = field_t[jmin:jmax, imin:imax]
                    if lim == True:
                        for i in range(len(x)):
                            for j in range(len(y)):
                                field_t[j,i] = max(min(field_t[j,i], lim_max), lim_min)
                    plt.contourf(x,y,field_t,cmap = plt.get_cmap('hot'))
                    if lim == True: ##IMPROVE
                        rg = np.arange(lim_min, lim_max+0.1, 0.1)
                        plt.colorbar(ticks=rg, label='digit value')
                        plt.clim(vmin = lim_min, vmax = lim_max)
                    else:
                        plt.colorbar(label='digit value')
                    tit = 't = '+ str(t) +'s'
                    plt.title(tit,fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                else:
                    field_t = field.transpose()
                    plt.contourf(self.x,self.y,field_t,cmap = plt.get_cmap('hot'))
                    plt.xlabel('x (m)'); plt.ylabel('y (m)')
                    plt.colorbar()
                                        
        def plot_field_video(self, field, name, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                j1='imagen00'
                j2='imagen0'
                j3='imagen'
                
                make_dir(self.p, 'images')
                path = self.p+'/images'
                clean_dir(path, name)
                make_dir(path, name)
                
                for t in range(len(self.t)):
                    figure = plt.figure(figsize = (5,10))
                    plt.grid(True)
                    field_t = field[t,:,:].transpose()
                    if clip == False:
                        x = self.x
                        y = self.y
                    elif clip == True:
                        imin = int((xmin-self.x[0])/(self.x[1]-self.x[0]))
                        imax = int((xmax-self.x[0])/(self.x[1]-self.x[0]))+1
                        jmin = int((ymin-self.y[0])/(self.y[1]-self.y[0]))
                        jmax = int((ymax-self.y[0])/(self.y[1]-self.y[0]))+1
                        x = self.x[imin:imax]
                        y = self.y[jmin:jmax]
                        field_t = field_t[jmin:jmax, imin:imax]
                    if lim == True:
                        for i in range(len(x)):
                            for j in range(len(y)):
                                field_t[j,i] = max(min(field_t[j,i], lim_max), lim_min)
                    plt.contourf(x,y,field_t,cmap = plt.get_cmap('hot'))
                    if lim == True: ##This is not working //AG
                        rg = np.arange(lim_min, lim_max+0.1, 0.1)
                        plt.colorbar(ticks=rg, label='digit value')
                        plt.clim(vmin = lim_min, vmax = lim_max)
                    else:
                        plt.colorbar(label='digit value')
                    tit = 't = '+ str(self.t[t]) +'s'
                    plt.title(tit,fontsize = 17); plt.xlabel('x (m)'); plt.ylabel('y (m)')
                   
                    if t>=10 and t<100:
                        d=j2+str(t)
                    elif t<10:
                        d=j1+str(t)
                    else:
                        d=j3+str(t)
                        
                    figure.savefig(self.p+'/images/'+name+'/'+d)
                    plt.close(figure)

                os.chdir(self.p+'/images/'+name)
                os.system("ffmpeg -framerate 3 -i imagen%03d.png {}.mp4".format(name))
                
                if clean == True:
                    shutil.copy(self.p+'/images/'+name+'/'+name+'.mp4', self.p+'/images/'+name+'.mp4')
                    clean_dir(path, name)
 
        def plot_h_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.h, 'h', clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_pb_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.pb, 'pb',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_Cv_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.Cv, 'Cv',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_deltah0_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.deltah0,'deltah0',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_deltac0_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.deltac0,'deltac0',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_Us_m_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.Us_m, 'Us_m',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_tau_m_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.tau_m, 'tau_m',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_Sm_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.Sm, 'Sm', lim,  clean, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_rho_video(self, clean = True, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field_video(self.rho, 'rho',  clean, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_h(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.h, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_Us_x(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Us[:,:,:,0], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_Us_y(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Us[:,:,:,1], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_Us_z(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Us[:,:,:,2], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_Us_m(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Us_m, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                                
        def plot_pb(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.pb, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
        
        def plot_tau_x(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.tau[:,:,:,0], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_tau_y(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.tau[:,:,:,1], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_tau_z(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.tau[:,:,:,2], lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_tau_m(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.tau_m, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)
                
        def plot_Cv(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Cv, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_deltah0(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.deltah0, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_deltac0(self, t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.deltac0, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_Sm(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.Sm, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)

        def plot_rho(self,t, lim = False, lim_min = 0, lim_max = 1, clip = False, xmin = 0, xmax = 0, ymin = 0, ymax = 0):
                self.plot_field(t,self.rho, lim, lim_min, lim_max, clip, xmin, xmax, ymin, ymax)


##########################################  fluxes functions  ##########################################################

##The following three functions where useful in the past.
#        def plot_face(self,index,plot = True):
#                face_points = self.face_points[index]
#                x_list = []
#                y_list = []
#                for p in face_points:
#                    x_list.append(p.x)
#                    y_list.append(p.y)
#                    
#                x_list.append(face_points[0].x)
#                y_list.append(face_points[0].y)
#                
#                if plot == True:
#                    fig = plt.figure(figsize = (10,5))
#                    plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                    plt.grid(True)
#                    plt.plot(y_list,x_list)
#                else:
#                    return [x_list,y_list]
#
#        def plot_faces(self, list_faces, zoom = True):
#                fig = plt.figure(figsize = (10,4))
#                plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                plt.grid(True)
#                if zoom != True:
#                    plt.ylim([self.x[0], self.x[-1]])
#                    plt.xlim([self.y[0], self.y[-1]])
#
#                for index in list_faces:
#                    [x_list,y_list] = self.plot_face(index,False)
#                    plt.plot(y_list,x_list)
#
#        def plot_faces_points(self, list_faces, list_points, zoom = True):
#                fig = plt.figure(figsize = (10,4))
#                plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                plt.grid(True)
#                if zoom != True:
#                    plt.ylim([self.x[0], self.x[-1]])
#                    plt.xlim([self.y[0], self.y[-1]])
#
#                for index in list_faces:
#                    [x_list,y_list] = self.plot_face(index,False)
#                    plt.plot(y_list,x_list)
#                    
#                for p in list_points:
#                    plt.plot(p.y,p.x, 'o', color = 'black')

        def get_phi2s(self, report = True, time = -1):
                if self.r.lp_name != '':
                    if (self.r.phi2s_flag == 'on' or self.r.phi2s_flag == 'yes' or self.r.phi2s_flag == True):
                        if report == True:
                            print('Creating phi2s ...')
                        if time == -1:
                            if self.r.phi2s == {}:
                                self.r.get_phi2s(False)
                        else:
                            if (str(time) in self.r.phi2s.keys()) == False:
                                self.r.get_phi2s(False, time)
                        self.get_scalar_flux(self.r.phi2s, self.phi2s, time)
                                    
        def get_Q(self, report = True, time = -1):
                if self.r.lp_name != '':
                    if (self.r.Q_flag == 'on' or self.r.Q_flag == 'yes' or self.r.Q_flag == True):
                        if report == True:
                            print('Creating Q ...')
                        if time == -1:
                            if self.r.Q == {}:
                                self.r.get_Q(False)
                        else:
                            if (str(time) in self.r.Q.keys()) == False:
                                self.r.get_Q(False, time)
                        self.get_scalar_flux(self.r.Q, self.Q, time)
                                   
        def get_scalar_flux(self, flux, dict, time):                
                for key in self.edges_list.keys():
                    if time == -1:
                        flux_lp = np.zeros((len(self.t),len(self.edges_list[key])))
                        for edges in range(len(self.edges_list[key])):                        
                            sub_flux = np.zeros((len(self.t),len(self.edges_list[key][edges])))
                            for t in range(len(self.t)):                    
                                sub_flux[t,:] = get_flux_tP(flux,self.edges_list[key][edges],self.t[t],self.edges_sense[key][edges])                        
                                sum_flux = 0
                                for k in range(sub_flux.shape[1]):
                                    sum_flux += sub_flux[t,k]
                                flux_lp[t,edges] = sum_flux
                    else:
                        flux_lp = np.zeros((1,len(self.edges_list[key])))
                        for edges in range(len(self.edges_list[key])):
                            sub_flux = get_flux_tP(flux,self.edges_list[key][edges],time,self.edges_sense[key][edges])
                            sum_flux = 0
                            for k in range(sub_flux.shape[1]):
                                sum_flux += sub_flux[0,k]
                            flux_lp[0,edges] = sum_flux                        
                                                          
                    dict[key] = flux_lp                                              

        def get_points_tps(self): #creates point tuples
                points = {}
                for key in self.tp.keys():
                    sub_points = []
                    for i in range(len(self.tp[key])):
                        sub_points.append([self.tp[key][i][0], self.tp[key][i][-1]])
                    points[key] = sub_points
                return points

        def get_faces_tp(self,point1,point2,delta): #returns a list of faces with the chain of faces that connects point1 and point2
                if point2.x != point1.x:
                    m = (point2.y-point1.y)/(point2.x-point1.x)
                    if point2.x > point1.x:
                        vpi = point1
                        vpf = point2
                    else:
                        vpi = point2
                        vpf = point1
                    list_points =[vpi]
                    list_faces = [self.get_face(vpi)]
                    while True:
                        vp = vector(list_points[-1].x+delta,m*delta+list_points[-1].y,0)
                        if vp.x > vpf.x:
                            if list_faces[-1] == self.get_face(vpf):
                                break
                            else:
                                delta = delta/2
                        else:
                            face_vp = self.get_face(vp)
                            if face_vp != list_faces[-1]:
                                if self.is_face_neighbour(face_vp,list_faces[-1]) == True:
                                    list_points.append(vp)
                                    list_faces.append(face_vp)
                                else:
                                    delta = delta/2
                            else:
                                delta = 1.5*delta
                else:
                    if point2.y > point1.y:
                        vpi = point1
                        vpf = point2
                    else:
                        vpi = point2
                        vpf = point1
                    list_points =[vpi]
                    list_faces = [self.get_face(vpi)]
                    while True:
                        vp = vector(list_points[-1].x, delta+list_points[-1].y,0)
                        if vp.y > vpf.y:
                            if list_faces[-1] == self.get_face(vpf):
                                break
                            else:
                                delta = delta/2
                        else:
                            face_vp = self.get_face(vp)
                            if face_vp != list_faces[-1]:
                                if self.is_face_neighbour(face_vp,list_faces[-1]) == True:
                                    list_points.append(vp)
                                    list_faces.append(face_vp)
                                else:
                                    delta = delta/2
                            else:
                                delta = 1.5*delta                    
                        
                return list_faces

        def get_points_tp(self,list_faces,point1,point2):#returns list_points, which is a list of the chain of points followed between point1 and point2
                #and list_faces with no duplicate indices
                if point2.x != point1.x:
                    m = (point2.y-point1.y)/(point2.x-point1.x)
                    if point2.x > point1.x:
                        vpi = point1
                        vpf = point2
                    else:
                        vpi = point2
                        vpf = point1
                else:
                    m = math.inf
                    if point2.y > point1.y:
                        vpi = point1
                        vpf = point2
                    else:
                        vpi = point2
                        vpf = point1  
                    
                list_points = []
                list_points.append(self.face_points[list_faces[0]][self.get_point_in_face(list_faces[0],vpi)])
                
                for i in range(len(list_faces)-1):
                    edge_index = self.get_common_edge(list_faces[i],list_faces[i+1])
                    point_1 = self.face_points[list_faces[i]][edge_index] #first point of the common edge
                    if edge_index != len(self.face_points[list_faces[i]])-1:
                        point_2 = self.face_points[list_faces[i]][edge_index+1] #second point of the common edge
                    else:
                        point_2 = self.face_points[list_faces[i]][0]

                    if get_closer_point(m,vpi,point_1,point_2) == 1:
                        point_f = point_1
                    elif get_closer_point(m,vpi,point_1,point_2) == 2:
                        point_f = point_2
                        
                    list_points.append(point_f)
                
                list_points.append(self.face_points[list_faces[-1]][self.get_point_in_face(list_faces[-1],vpf)])
                
                i = 0
                while True:
                    if list_points[i] == list_points[i+1]:
                        point = list_points.pop(i)
                        face = list_faces.pop(i)
                    else:
                        i += 1
                    if i == len(list_points)-1:
                        break
                                          
                return [list_points,list_faces]
      
        def get_edges_tp(self,list_faces,list_points): #we start with a non-connected list of points, the result is a connected list of points and edges
                list_edges = []
                new_list_points = []
                new_list_points.append(list_points[0])
                for i in range(len(list_faces)):
                    index_point_i = self.get_point_in_face(list_faces[i],list_points[i])
                    index_point_f = self.get_point_in_face(list_faces[i],list_points[i+1])
                    
                    if index_point_i != index_point_f:
                        n = len(self.face_points[list_faces[i]])
                        sense = choose_sense(n,index_point_i,index_point_f)
                        if sense == 1:
                            while True:
                                list_edges.append(self.edge_faces[list_faces[i]][index_point_i])
                                index_point_i +=1
                                if index_point_i == n:
                                    index_point_i = 0
                                new_list_points.append(self.face_points[list_faces[i]][index_point_i])
                                if index_point_i == index_point_f:
                                    break
                        else:
                            while True:
                                index_point_i -=1
                                if index_point_i ==-1:
                                    index_point_i = n-1
                                list_edges.append(self.edge_faces[list_faces[i]][index_point_i])
                                new_list_points.append(self.face_points[list_faces[i]][index_point_i])
                                if index_point_i == index_point_f:
                                    break
                                
                return [list_edges, new_list_points] 

        def get_edges_sense(self,list_edges,list_points, dx): #returns the sense of each edge in list_edges
                edges_sense = []
                for i in range(len(list_points)-1):
                    face_owner = int(self.r.eO['edges'][list_edges[i]])
                    face_neighbour = int(self.r.eN['edges'][list_edges[i]])
                    dx_i = dx
                    while True:
                        point_t = get_point_t(list_points[i], list_points[i+1], dx_i)
                        face = self.get_face(point_t)
                        if face  == -1:
                            dx_i *= 1.5
                        else:
                            if face == face_owner:
                                sense = 1
                                break
                            elif face == face_neighbour:
                                sense = -1
                                break
                            else:
                                dx_i /= 1.2
                    edges_sense.append(sense)
                return edges_sense

        def get_face(self,v): #return the face to which a point v belongs to
                face_index = -1
                [i,j] = self.get_coords(v)
                list_index = []
                for k in range(len(self.index_matrix[i][j])):
                    if not(self.index_matrix[i][j][k] in list_index):
                        list_index.append(self.index_matrix[i][j][k])
                i +=1
                for k in range(len(self.index_matrix[i][j])):
                    if not(self.index_matrix[i][j][k] in list_index):
                        list_index.append(self.index_matrix[i][j][k])
                j +=1
                for k in range(len(self.index_matrix[i][j])):
                    if not(self.index_matrix[i][j][k] in list_index):
                        list_index.append(self.index_matrix[i][j][k])
                i -=1
                for k in range(len(self.index_matrix[i][j])):
                    if not(self.index_matrix[i][j][k] in list_index):
                        list_index.append(self.index_matrix[i][j][k])
                poligons = []
                for k in range(len(list_index)):
                    poligons.append(self.face_points[list_index[k]])
                for k in range(len(poligons)):
                    if face_in(poligons[k],v) == True:
                        face_index = list_index[k]
                        return face_index
                if face_index == -1:
                    return face_index
                      
        def get_coords(self,vector): #return the (i,j) indices of a vector using the previously defined grid
                i = int((vector.x-self.x[0])/(self.x[1]-self.x[0]))
                j = int((vector.y-self.y[0])/(self.y[1]-self.y[0]))
                return [i,j]

        def get_point_in_face(self,index,v): #knowing that v belogs to index, which is a face, we determine which is the closest point in the face (index) to v
                list = self.face_points[index]
                dist = 10**6
                for i in range(len(list)):
                    if distxy(list[i],v)<=dist:
                        dist = distxy(list[i],v)
                        point_index = i
                return point_index
            
        def is_face_neighbour(self,index1,index2): #returns true if the face index1 is a neighbour of the face index2
                if index2 in self.face_neighbour[index1]:
                    return True
                else:
                    return False

        def get_common_edge(self,index1,index2): #returns the index of the common edge between two faces
                for i in range(len(self.edge_faces[index1])):
                    if self.edge_faces[index1][i] in self.edge_faces[index2]:
                        return i


########################################  fluxes plot functions  #########################################################

        def plot_Q_tp(self,key_list,index,deltax_axis = 20,deltay_axis = 20,save = False):
                self.plot_lpfield_tp(key_list,index,self.Q,'Q','m³/s',deltax_axis,deltay_axis,save)

        def plot_phi2s_tp(self,key_list,index,deltax_axis = 20,deltay_axis = 20,save = False):
                self.plot_lpfield_tp(key_list,index,self.phi2s,'\u03C1'+ 'Q', 'kg/s',deltax_axis,deltay_axis,save)

        def plot_lpfield_tp(self, key_list, index, field_dict, name, units, deltax_axis, deltay_axis,save):
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.yaxis.set_major_locator(MultipleLocator(deltay_axis))
                ax.xaxis.set_major_locator(MultipleLocator(deltax_axis))
                field = np.zeros((len(self.t),1))
                for key in key_list:
                    if key != key_list[0]:
                        if self.lp[key_list[0]]['list_coords'][1:-1][index] != self.lp[key]['list_coords'][1:-1][index]:
                            field[:,0] += field_dict[key][:,index]
                    else:
                        field[:,0] += field_dict[key][:,index]                  
                plt.plot(self.t ,abs(field[:,0]))
                plt.xlim([0, self.t[-1]])
                plt.ylim([0,  1.05*max(abs(field[:,0]))])
                plt.grid(True)
                plt.title('index = '+str(index),fontsize = 17); plt.xlabel('t (s)'); plt.ylabel(name +' ('+ units+')')

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+name+'_tp_time_index_'+str(index))
                    plt.close(fig)  

        def plot_phi2s(self,t,key_list,deltax_axis = 10,deltay_axis = 20):
                self.plot_lpfield(t, key_list, self.phi2s, '\u03C1'+ 'Q', 'kg/s', deltax_axis, deltay_axis)

        def plot_Q(self,t,key_list,deltax_axis = 10,deltay_axis = 20):
                self.plot_lpfield(t, key_list, self.Q, 'Q', 'm³/s', deltax_axis, deltay_axis)

        def plot_lpfield(self,t, key_list, field_dict, name, units, deltax_axis, deltay_axis):
                fig, ax = plt.subplots(figsize=(10, 5))
                j = self.t.index(t)
                ax.yaxis.set_major_locator(MultipleLocator(deltay_axis))
                ax.xaxis.set_major_locator(MultipleLocator(deltax_axis))
                field = np.zeros((len(self.lp[key_list[0]]['list_dists'][1:-1]),1))
                for key in key_list:
                    for i in range(len(field)):
                        if key != key_list[0]:
                            if self.lp[key_list[0]]['list_coords'][1:-1][i] != self.lp[key]['list_coords'][1:-1][i]:
                                field[i,0] += field_dict[key][j,i]
                        else:
                            field[i,0] += field_dict[key][j,i]                            
                plt.plot(self.lp[key_list[0]]['list_dists'][1:-1] ,abs(field[:,0]))
                plt.xlim([0, self.lp[key_list[0]]['list_dists'][-1]])
                plt.ylim([0,  1.05*max(abs(field[:,0]))])
                plt.grid(True)
                plt.title('t = '+str(t)+ 's',fontsize = 17); plt.xlabel('l (m)'); plt.ylabel(name +' ('+ units+')')

        def plot_phi2s_video(self,key_list,deltax_axis = 10,deltay_axis = 20, clean = True):
                self.plot_lpfield_video(key_list, self.phi2s, '\u03C1'+ 'Q', 'kg/s', deltax_axis, deltay_axis, clean)

        def plot_Q_video(self,key_list,deltax_axis = 10,deltay_axis = 20, clean = True):
                self.plot_lpfield_video(key_list, self.Q, 'Q', 'm³/s', deltax_axis, deltay_axis, clean)
 
        def plot_lpfield_video(self, key_list, field_dict, name, units, deltax_axis, deltay_axis, clean):
                j1='imagen00'
                j2='imagen0'
                j3='imagen'
                
                make_dir(self.p, 'images')
                path = self.p+'/images'
                clean_dir(path, name)
                make_dir(path, name)
                
                field = np.zeros((len(self.t),len(self.lp[key_list[0]]['list_dists'][1:-1])))

                for key in key_list:
                    for i in range(len(self.lp[key_list[0]]['list_dists'][1:-1])):
                        if key != key_list[0]:
                            if self.lp[key_list[0]]['list_coords'][1:-1][i] != self.lp[key]['list_coords'][1:-1][i]:
                                field[:,i] += field_dict[key][:,i]
                        else:
                            field[:,i] += field_dict[key][:,i] 

                field_max = np.amax(abs(field))
                
                for t in range(len(self.t)):
                    fig, ax = plt.subplots(figsize=(10, 5))
                    plt.grid(which='major')
                    ax.yaxis.set_major_locator(MultipleLocator(deltay_axis))
                    ax.xaxis.set_major_locator(MultipleLocator(deltax_axis))
                    
                    plt.plot(self.lp[key_list[0]]['list_dists'][1:-1] ,abs(field[t,:]))
                    plt.title('t = '+str(self.t[t])+ 's',fontsize = 17); plt.xlabel('l (m)'); plt.ylabel(name +' ('+ units+')')
                    plt.ylim([0, field_max*1.05])
                    plt.xlim([0, self.lp[key_list[0]]['list_dists'][-1]])
                    
                    if t>=10 and t<100:
                        d=j2+str(t)
                    elif t<10:
                        d=j1+str(t)
                    else:
                        d=j3+str(t)
                        
                    fig.savefig(self.p+'/images/'+name+'/'+d)
                    plt.close(fig)

                os.chdir(self.p+'/images/'+name)
                os.system("ffmpeg -framerate 4 -i imagen%03d.png {}.mp4".format(name))

                if clean == True:
                    shutil.copy(self.p+'/images/'+name+'/'+name+'.mp4', self.p+'/images/'+name+'.mp4')
                    clean_dir(path, name)

##These functions were used in the past to check the correct behaviour of the code
#        def plot_edges_tp(self, key):
#                figure = plt.figure(figsize = (20,4))
#                plt.grid(True)
#                plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                for i in range(len(self.Q[key]['list_points_tp'])):
#                    sub_x = []
#                    sub_y = []
#                    for j in range(len(self.Q[key]['list_points_tp'][i])):
#                        sub_x.append(self.Q[key]['list_points_tp'][i][j].x)
#                        sub_y.append(self.Q[key]['list_points_tp'][i][j].y)
#                    plt.plot(sub_y,sub_x)
#                    
#        def plot_edge_tp(self, key, list):
#                figure = plt.figure(figsize = (8,8))
#                plt.grid(True)
#                plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                for i in list:
#                    sub_x = []
#                    sub_y = []
#                    for j in range(len(self.Q[key]['list_points_tp'][i])):
#                        sub_x.append(self.Q[key]['list_points_tp'][i][j].x)
#                        sub_y.append(self.Q[key]['list_points_tp'][i][j].y)
#                    plt.plot(sub_y,sub_x)
#
#        def plot_edges_tp_faces(self, key, list_faces, list_edges, zoom = True):
#                fig = plt.figure(figsize = (10,4))
#                plt.xlabel('y (m)'); plt.ylabel('x (m)')
#                plt.grid(True)
#                if zoom != True:
#                    plt.ylim([self.x[0], self.x[-1]])
#                    plt.xlim([self.y[0], self.y[-1]])
#
#                for index in list_faces:
#                    [x_list,y_list] = self.plot_face(index,False)
#                    plt.plot(y_list,x_list)
#
#                for i in list_edges:
#                    sub_x = []
#                    sub_y = []
#                    for j in range(len(self.Q[key]['list_points_tp'][i])):
#                        sub_x.append(self.Q[key]['list_points_tp'][i][j].x)
#                        sub_y.append(self.Q[key]['list_points_tp'][i][j].y)
#                    plt.plot(sub_y,sub_x)


########################################  some extra functions  #########################################################

        def get_Sm(self, report = True, time = -1): # Entrainment rate as a scalar field  
                if (self.Sm_flag == 'on' or self.Sm_flag == 'yes' or self.Sm_flag == True):
                    if report == True:
                        print('Calculating Sm field ...')
                        
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)    
                        
                    if self.r.eM != 'entrainmentOff':                        
                        if time == -1:
                            self.Sm = np.zeros((len(self.t),len(self.x),len(self.y)))
                        else:
                            self.Sm = np.zeros((len(self.x),len(self.y)))
                                            
                        if len(self.Us) == 0:
                            self.get_Us_interpolation(False, time)
                        if len(self.h) == 0:
                            self.get_h_interpolation(False, time)

## It does not work for a case with uniform entrainmentZone and it have not been tested for a case with two or more zones //AG                        
                        if self.r.eM == 'Exponential':
                            E = self.r.eM_coeffs['E']['value']

                            if self.r.eM_coeffs == {}:
                                self.r.get_entrainment_coeffs(False)
                            eZlist = []
                            for element in self.r.eM_coeffs['entrainmentZones']:
                                if type(element) == dict:
                                    eZlist.append(element)
                        
                            poligons = []
                            for i in range(len(eZlist)):
                                off = eZlist[i]['offset'][0][:-1]
                                polygon = eZlist[i]['vertices']
                                rectangle = []
                                
                                for j in range(len(polygon)):
                                    rectangle.append(vector(polygon[j][0]-off[0], polygon[j][1]-off[1], 0))
                                poligons.append(rectangle)    
                       
                            if time == -1:
                                for t in range(len(self.t)):
                                    for i in range(len(self.x)):
                                        for j in range(len(self.y)):
                                            ez = 0
                                            v = vector(self.x[i],self.y[j],0)
                                            for k in range(len(poligons)):
                                                if face_in(poligons[k],v) == True:
                                                    ez = 1
                                            self.Sm[t,i,j] = E*self.h[t,i,j]*self.Us_m[t,i,j]*ez
                            else:
                                for i in range(len(self.x)):
                                        for j in range(len(self.y)):
                                            ez = 0
                                            v = vector(self.x[i],self.y[j],0)
                                            for k in range(len(poligons)):
                                                if face_in(poligons[k],v) == True:
                                                    ez = 1
                                            self.Sm[i,j] = E*self.h[i,j]*self.Us_m[i,j]*ez
                                    
                    self.r.clean_eM()

        def get_rho(self, report = True, time = -1):  # density as a scalar field 
                if (self.rho_flag == 'on' or self.rho_flag == 'yes' or self.rho_flag == True):
                    if report == True:
                        print('Calculating rho field ...')
                    if time == -1:
                        self.rho = np.zeros((len(self.t),len(self.x),len(self.y)))
                    else:
                        self.rho = np.zeros((len(self.x),len(self.y)))
                    
                    if self.r.rho_w == -1 or self.r.rho_s == -1:
                        self.r.get_densities(False)                    
                    
                    if len(self.Cv) == 0:
                        Cv_flag = self.r.Cv_flag
                        self.r.Cv_flag = 'on'
                        self.get_Cv_interpolation(False, time)
                        self.r.Cv_flag = Cv_flag
                                        
                    if time == -1:
                        for t in range(len(self.t)):
                            for i in range(len(self.x)):
                                for j in range(len(self.y)):
                                    self.rho[t,i,j] = self.r.rho_w+self.Cv[t,i,j]*(self.r.rho_s-self.r.rho_w)
                    else:
                        for i in range(len(self.x)):
                                for j in range(len(self.y)):
                                    self.rho[i,j] = self.r.rho_w+self.Cv[i,j]*(self.r.rho_s-self.r.rho_w)

        def get_field_t(self, time, field, func):
                if (str(time) in field.keys()) == False:
                    func(False, time)

                if field[str(time)].iF.u == 'uniform':
                    f = []
                    f = field[str(time)].iF.field*self.r.nF
                else:
                    f = field[str(time)].iF.field
                    
                return f

        def get_rcg(self, report = True, time = -1):  #center of gravity of the whole fluid     
                if (self.rcg_flag == 'on' or self.rcg_flag == 'yes' or self.rcg_flag == True):
                    if report == True:
                        print('Calculating center of gravity of the flux in the whole domain ...')
                    if time == -1:
                        self.rcg = np.zeros((len(self.t),3))
                    else:
                        self.rcg = np.zeros((1,3))
                    
                    if self.r.rho_w == -1 or self.r.rho_s == -1:
                        self.r.get_densities(False)                    
                    if self.r.n == {}:
                        self.r.get_n(False,0)
                    if self.r.A == {}:
                        self.r.get_A(False,0)
                    if self.r.c == {}:
                        self.r.get_c(False,0)
 
                    n = self.r.n['0'].iF.field
                    c = self.r.c['0'].iF.field
                    A = self.r.A['0'].iF.field
                                                           
                    if time == -1:
                        for t in range(len(self.t)):
                            h = self.get_field_t(self.t[t],self.r.h, self.r.get_h)
                            Cv = self.get_field_t(self.t[t],self.r.Cv, self.r.get_Cv)
                                                                                    
                            sum_num = vector(0,0,0)
                            sum_den = 0
                            for i in range(self.r.nF):
                                sum_den += (self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*h[i]
                                sum_num += (c[i]-n[i]*1/2*h[i])*(self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*h[i]
                                    
                            self.rcg[t,0] = sum_num.x/sum_den
                            self.rcg[t,1] = sum_num.y/sum_den
                            self.rcg[t,2] = sum_num.z/sum_den
                    else:                        
                        h = self.get_field_t(time,self.r.h, self.r.get_h)
                        Cv = self.get_field_t(time,self.r.Cv, self.r.get_Cv)
                                                
                        sum_num = vector(0,0,0)
                        sum_den = 0
                        for i in range(self.r.nF):
                            sum_den += (self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*h[i]
                            sum_num += (c[i]-n[i]*1/2*h[i])*(self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*h[i]
                                    
                        self.rcg[0,0] = sum_num.x/sum_den
                        self.rcg[0,1] = sum_num.y/sum_den
                        self.rcg[0,2] = sum_num.z/sum_den

        def get_M(self, report = True, time = -1): #Mass of the whole fluid   
                if (self.M_flag == 'on' or self.M_flag == 'yes' or self.M_flag == True):
                    if report == True:
                        print('Calculating Mass of flux in the whole domain ...')
                    if time == -1:
                        self.M = np.zeros((len(self.t),1))
                    else:
                        self.M = np.zeros((1,1))

                    if self.r.rho_w == -1 or self.r.rho_s == -1:
                        self.r.get_densities(False)                    
                    if self.r.A == {}:
                        self.r.get_A(False,0)
                    if self.r.hmin == -1:
                        self.r.get_hmin(False)
                        
                    A = self.r.A['0'].iF.field
                    
                    if time == -1:
                        for t in range(len(self.t)):
                            h = self.get_field_t(self.t[t],self.r.h, self.r.get_h)
                            Cv = self.get_field_t(self.t[t],self.r.Cv, self.r.get_Cv)
                            
                            sum = 0
                            for i in range(self.r.nF):
                                sum += (self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*(h[i]-self.r.hmin)
                                    
                            self.M[t,0] = sum
                    else:
                        h = self.get_field_t(time,self.r.h, self.r.get_h)
                        Cv = self.get_field_t(time,self.r.Cv, self.r.get_Cv)
                       
                        sum = 0
                        for i in range(self.r.nF):
                            sum += (self.r.rho_w+Cv[i]*(self.r.rho_s-self.r.rho_w))*A[i]*(h[i]-self.r.hmin)
                                    
                        self.M[0,0] = sum

        def get_V(self, report = True, time = -1):        
                if (self.V_flag == 'on' or self.V_flag == 'yes' or self.V_flag == True):
                    if report == True:
                        print('Calculating Volume of flux in the whole domain ...')
                    if time == -1:
                        self.V = np.zeros((len(self.t),1))
                    else:
                        self.V = np.zeros((1,1))

                    if self.r.A == {}:
                        self.r.get_A(False,0)
                    if self.r.hmin == -1:
                        self.r.get_hmin(False)
                        
                    A = self.r.A['0'].iF.field

                    if time == -1:
                        for t in range(len(self.t)):
                            h = self.get_field_t(self.t[t],self.r.h, self.r.get_h)
                                                        
                            sum = 0
                            for i in range(self.r.nF):
                                sum += A[i]*(h[i]-self.r.hmin)
                                    
                            self.V[t,0] = sum
                    else:
                        h = self.get_field_t(time,self.r.h, self.r.get_h)
                                                    
                        sum = 0
                        for i in range(self.r.nF):
                            sum += A[i]*(h[i]-self.r.hmin)
                                    
                        self.V[0,0] = sum

        def get_Vsed(self, report = True, time = -1):        
                if (self.Vsed_flag == 'on' or self.Vsed_flag == 'yes' or self.Vsed_flag == True):
                    if report == True:
                        print('Calculating Volume of sediment eroded in the whole domain ...')
                    if time == -1:
                        self.Vsed = np.zeros((len(self.t),1))
                    else:
                        self.Vsed = np.zeros((1,1))

                    if self.r.A == {}:
                        self.r.get_A(False,0)
                    if self.r.tps == {}:
                        self.r.get_transportProperties(False)
                        
                    A = self.r.A['0'].iF.field
                                                
                    if time == -1:
                        for t in range(len(self.t)):
                            if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                                deltaz0 = self.get_field_t(self.t[t],self.r.deltac0, self.r.get_deltac0)

                                sum = 0
                                for i in range(self.r.nF):
                                    sum += A[i]*deltaz0[i]
                                        
                                self.Vsed[t,0] = sum

                            elif self.r.tps['entrainmentModel'] != 'entrainmentOff' and self.r.tps['depositionModel'] != 'depositionOff':
                                deltaz0 = self.get_field_t(self.t[t],self.r.deltah0, self.r.get_deltah0)
                                                        
                                sum = 0
                                for i in range(self.r.nF):
                                    sum += A[i]*deltaz0[i]
                                        
                                self.Vsed[t,0] = sum
                    else:
                        if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                            deltaz0 = self.get_field_t(time,self.r.deltac0, self.r.get_deltac0)

                            sum = 0
                            for i in range(self.r.nF):
                                sum += A[i]*deltaz0[i]
                                        
                            self.Vsed[0,0] = sum

                        elif self.r.tps['entrainmentModel'] != 'entrainmentOff' and self.r.tps['depositionModel'] != 'depositionOff':
                            deltaz0 = self.get_field_t(time,self.r.deltah0, self.r.get_deltah0)
                                                    
                            sum = 0
                            for i in range(self.r.nF):
                                sum += A[i]*deltaz0[i]
                                        
                            self.Vsed[0,0] = sum
                        
                    self.r.clean_tps()

        def get_rF_VyM(self, report = True):
                if report == True:
                        print('Calculating Volume and Mass entering to the geometric domain ...')
                self.rF_V = np.zeros((len(self.t),1))
                self.rF_M = np.zeros((len(self.t),1))
                if self.r.rF == {}:
                    self.r.get_releaseFlow(False)
                self.r.get_densities(False)

                for element in self.r.rF['boundaries']:
                    if type(element) is dict:
                        hidrogram = np.zeros((len(element['Q']),6))
                        for i in range(len(hidrogram)):
                            hidrogram[i,0] = element['Q'][i][0]
                            hidrogram[i,1] = element['Q'][i][1]
                            hidrogram[i,2] = element['Q'][i][2]
                            
                        for i in range(len(hidrogram)-1):
                            hidrogram[i+1,3] = hidrogram[i,3]+1/2*(hidrogram[i+1,0]-hidrogram[i,0])*(hidrogram[i+1,1]+hidrogram[i,1])
                            hidrogram[i+1,4] = hidrogram[i,4]+1/3*(hidrogram[i+1,2]-hidrogram[i,2])*(hidrogram[i+1,1]-hidrogram[i,1])*(hidrogram[i+1,0]-hidrogram[i,0])+1/2*(hidrogram[i+1,0]-hidrogram[i,0])*(hidrogram[i+1,1]*hidrogram[i,2]+hidrogram[i,1]*hidrogram[i+1,2])
                            hidrogram[i+1,5] = self.r.rho_w*hidrogram[i+1,3]+(self.r.rho_s-self.r.rho_w)*hidrogram[i+1,4]
                                            
                        for t in range(len(self.t)-1):
                            tt = self.t[t+1]
                            i = 0
                            j = len(hidrogram)-1
                            
                            if tt <= hidrogram[j,0]:
                                while (j-i>1):
                                    n = int(i+(j-i)/2)
                                    if  hidrogram[n,0]==tt:
                                        i = n
                                        j = n
                                    elif hidrogram[n,0]>tt:
                                        j = n
                                    else:
                                        i = n
                                Qt = hidrogram[i,1]+(hidrogram[j,1]-hidrogram[i,1])/(hidrogram[j,0]-hidrogram[i,0]+10**-6)*(tt-hidrogram[i,0])
                                Cvt = hidrogram[i,2]+(hidrogram[j,2]-hidrogram[i,2])/(hidrogram[j,0]-hidrogram[i,0]+10**-6)*(tt-hidrogram[i,0])
                            else:
                                Qt = hidrogram[j,1]
                                Cvt = hidrogram[j,2]          
            
                            self.rF_V[t+1,0] = hidrogram[i,3]+1/2*(tt-hidrogram[i,0])*(Qt+hidrogram[i,1])
                            self.rF_M[t+1,0] = hidrogram[i,5]+(self.r.rho_s-self.r.rho_w)*(tt-hidrogram[i,0])*(1/3*(Qt-hidrogram[i,1])*(Cvt-hidrogram[i,2])+1/2*(Qt*hidrogram[i,2]+Cvt*hidrogram[i,1]))+self.r.rho_w*(1/2*(tt-hidrogram[i,0])*(Qt+hidrogram[i,1]))

                self.r.rF = {}


########################################  some extra plot functions  #########################################################

        def plot_M(self,firstvalue= 0,save = False):
                self.plot_timefield(self.M, 'M', 'kg', firstvalue, save)

        def plot_V(self,firstvalue= 0,save = False):
                self.plot_timefield(self.V, 'V', 'm³', firstvalue, save)
                
        def plot_Vsed(self,firstvalue= 0,save = False):
                self.plot_timefield(self.Vsed, 'Vsed', 'm³', firstvalue, save)

        def plot_rcg_x(self,firstvalue= 0,save = False):
                self.plot_timefield(self.rcg[:,0], 'rcg_x', 'm', firstvalue, save, True, self.x[0], self.x[-1])

        def plot_rcg_y(self,firstvalue= 0,save = False):
                self.plot_timefield(self.rcg[:,1], 'rcg_y', 'm', firstvalue, save, True, self.y[0], self.y[-1])

        def plot_rcg_z(self,firstvalue= 0,save = False):
                self.plot_timefield(self.rcg[:,2], 'rcg_z', 'm', firstvalue, save)                
                
        def plot_timefield(self, field, name, units, firstvalue, save, clip = False, ymin = 0, ymax = 0):
                fig, ax = plt.subplots(figsize=(5, 10))
                plt.plot(self.t[firstvalue:] ,field[firstvalue:])
                plt.xlim([0, self.t[-1]])
                if clip == True:
                    plt.ylim([ymin,  ymax])
                plt.grid(True)
#                plt.title('',fontsize = 17); 
                plt.xlabel('t (s)'); plt.ylabel(name +' ('+ units+')')

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+name)
                    plt.close(fig)

        def plot_M_net(self, save = False, clip = False, ymin = 0, ymax = 0):
                if self.r.rho_b == -1:
                    self.r.get_densities(False)
                fig, ax = plt.subplots(figsize=(5, 10))
                line1, = plt.plot(self.t ,self.rF_M[:,0], label = 'input mass')
                line2, = plt.plot(self.t ,self.M[:,0], label = 'total mass')
                line3, = plt.plot(self.t ,-self.Vsed[:,0]*self.r.rho_b, label = 'estimated incorporated mass')
                line4, = plt.plot(self.t ,self.M[:,0]-self.rF_M[:,0], label = 'incorporated mass')

                plt.xlim([0, self.t[-1]])
                if clip == True:
                    plt.ylim([ymin,  ymax])
                plt.grid(True)
#                plt.title('',fontsize = 17); 
                plt.xlabel('t (s)'); plt.ylabel('m (kg)')
                
                plt.legend(handler_map= {})

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+'M_net')
                    plt.close(fig)            

        def plot_V_net(self, save = False, clip = False, ymin = 0, ymax = 0):
                fig, ax = plt.subplots(figsize=(5, 10))
                line1, = plt.plot(self.t ,self.rF_V[:,0], label = 'input volume')
                line2, = plt.plot(self.t ,self.V[:,0], label = 'total volume')
                line3, = plt.plot(self.t ,-self.Vsed[:,0], label = 'estimated incorporated volume')
                line4, = plt.plot(self.t ,self.V[:,0]-self.rF_V[:,0], label = 'incorporated volume')

                plt.xlim([0, self.t[-1]])
                if clip == True:
                    plt.ylim([ymin,  ymax])
                plt.grid(True)
#                plt.title('',fontsize = 17); 
                plt.xlabel('t (s)'); plt.ylabel('V (m³)')
                
                plt.legend(handler_map= {})

                if save == True:
                    path = self.p+'/images'
                    fig.savefig(path+'/'+'V_net')
                    plt.close(fig) 


########################################  writing functions  #########################################################

        def create_t_size(self, size):
                for i in range(size):
                    sub_list = []
                    for j in range(int(len(self.t)/size)):
                        sub_list.append(self.t[j*size+i])
                    self.t_size.append(sub_list)
                for i in range(len(self.t)-int(len(self.t)/size)*size):
                    self.t_size[i+1].append(self.t[int(len(self.t)/size)*size+i])            

        def write_output(self, path, size, time = -1):
                if self.rank == 0:
                    print('Writing output files ...')
                    self.write_t(path)
                    self.write_x(path)
                    self.write_y(path)
                    self.write_z(path)
                    self.write_lp(path)
                    self.write_tp(path)
                    self.write_he(path)
                    self.write_rF_VyM(path)
                    self.write_z_tp(path)
                    self.write_he_tp(path)
                self.create_dir(path, size, time)
                self.write_fields(path, size, time)              

        def create_dir(self, path, size, time):
                if size == 1:
                    if time == -1:
                        for t in self.t:
                            os.mkdir(path+'/'+str(t))
                    else:
                        os.mkdir(path+'/'+str(time))
                else:
                    if len(self.t_size) == 0:
                        self.create_t_size(size)
                        
                    if time == -1:
                        for t in self.t_size[self.rank]:
                            os.mkdir(path+'/'+str(t))
                        
        def write_x(self, path):
                self.write_list(path, 'x', self.x)

        def write_y(self, path):
                self.write_list(path, 'y', self.y)

        def write_z(self, path):
                if len(self.z) == 0:
                    self.get_z_interpolation()
                fi = open(path+'/'+'z', 'w')
                self.write_field(fi, self.z)
                fi.close()

        def write_he(self, path):
                if len(self.he) == 0:
                    self.get_he_interpolation()
                fi = open(path+'/'+'he', 'w')
                self.write_field(fi, self.he)
                fi.close()

        def write_t(self, path):
                self.write_list(path, 't', self.t)

        def write_rF_VyM(self, path):
                if len(self.rF_M) == 0 or len(self.rF_V) == 0:
                    self.get_rF_VyM()
                self.write_list(path, 'rF_M', self.rF_M[:,0])
                self.write_list(path, 'rF_V', self.rF_V[:,0])

        def write_lp(self, path):
                if self.lp != {}:
                    fi = open(path+'/'+'lp', 'w')
                    for key in self.lp.keys():
                        self.write_lpfield(fi, self.lp[key]['list_dists'][1:-1], key)
                    fi.close()

                    lp_c = {}
                    for key in self.lp.keys():
                        sub_list_x = []
                        sub_list_y = []
                        for i in range(len(self.lp[key]['list_coords'])):
                            sub_list_x.append(self.lp[key]['list_coords'][i][0])
                            sub_list_y.append(self.lp[key]['list_coords'][i][1])
                        lp_c[key] = [sub_list_x, sub_list_y]
                    
                    fi = open(path+'/'+'lp_cx', 'w')
                    for key in self.lp.keys():
                        self.write_lpfield(fi, lp_c[key][0], key)
                    fi.close()  
                    
                    fi = open(path+'/'+'lp_cy', 'w')
                    for key in self.lp.keys():
                        self.write_lpfield(fi, lp_c[key][1], key)
                    fi.close()  

        def write_tp(self, path, nround = 3):
                if self.tp != {}:
                    fi = open(path+'/'+'tp_x', 'w')
                    for key in self.tp.keys():
                        fi.write(str(key)+"\n")
                        fi.write("{"+"\n")
                        for i in range(len(self.tp[key])):
                            for j in range(len(self.tp[key][i])):
                                 fi.write(str(round(self.tp[key][i][j].x, nround))+" ")
                            fi.write("\n") 
                        fi.write("}"+"\n")
                        fi.write("\n")
                    fi.close()
                    fi = open(path+'/'+'tp_y', 'w')
                    for key in self.tp.keys():
                        fi.write(str(key)+"\n")
                        fi.write("{"+"\n")
                        for i in range(len(self.tp[key])):
                            for j in range(len(self.tp[key][i])):
                                 fi.write(str(round(self.tp[key][i][j].y, nround))+" ")
                            fi.write("\n") 
                        fi.write("}"+"\n")
                        fi.write("\n")
                    fi.close()

        def write_fields(self, path, size, time):
                import time as timex
                if time == -1:
                    if size == 1:
                        for t in self.t:
                            print('   Writing output fields for time = '+ str(t)+ ' at time = ' + date(timex.gmtime()))
                            self.write_areaFields(path,t)
                            self.write_edgeFields(path,t)
                            self.write_tpFields(path,t)
                            self.write_timeFields(path,t)
                            self.r.clean_fields(False, t)
                    else:
                        for t in self.t_size[self.rank]:
                            print('   Writing output fields for time = '+ str(t)+ ' at time = ' + date(timex.gmtime())+ ' with rank = ' + str(self.rank))
                            self.write_areaFields(path,t)
                            self.write_edgeFields(path,t)
                            self.write_tpFields(path,t)
                            self.write_timeFields(path,t)
                            self.r.clean_fields(False, t)
                else:
                    print('   Writing output fields for time = '+ str(time)+ ' at time = ' + date(timex.gmtime()))
                    self.write_areaFields(path,time)
                    self.write_edgeFields(path,time)
                    self.write_tpFields(path,time)
                    self.write_timeFields(path,time)
                    self.r.clean_fields(False, time)

        def write_areaFields(self, path, t):
                    self.write_h(path,t)
                    self.write_pb(path,t)
                    self.write_Cv(path,t)
                    self.write_deltaz0(path,t)
                    self.write_Us(path,t)
                    self.write_tau(path,t)
                    self.write_Sm(path,t)
                    self.write_rho(path,t)

        def write_edgeFields(self, path, t):
                    self.write_Q(path,t)
                    self.write_phi2s(path,t)

        def write_tpFields(self, path, t):
                    self.write_h_tp(path,t)
                    self.write_deltaz0_tp(path,t)
                    self.write_Us_m_tp(path,t)

        def write_timeFields(self, path, t):
                    self.write_rcg(path,t)
                    self.write_M(path,t)
                    self.write_V(path,t)
                    self.write_Vsed(path,t)

        def write_rcg(self, path, time):
                if (self.rcg_flag == 'on' or self.rcg_flag == 'yes' or self.rcg_flag == True):
                    self.get_rcg(False, time)
                    self.write_vectortimefield(path, 'rcg', self.rcg, time)

        def write_M(self, path, time):
                if (self.M_flag == 'on' or self.M_flag == 'yes' or self.M_flag == True):
                    self.get_M(False, time)
                    self.write_scalartimefield(path, 'M', self.M, time)

        def write_V(self, path, time):
                if (self.V_flag == 'on' or self.V_flag == 'yes' or self.V_flag == True):
                    self.get_V(False, time)
                    self.write_scalartimefield(path, 'V', self.V, time)

        def write_Vsed(self, path, time):
                if (self.Vsed_flag == 'on' or self.Vsed_flag == 'yes' or self.Vsed_flag == True):
                    if self.r.tps == {}:
                        self.r.get_transportProperties(False)
                    if (self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True) or (self.r.tps['entrainmentModel'] != 'entrainmentOff' and self.r.tps['depositionModel'] != 'depositionOff'):
                        self.get_Vsed(False, time)
                        self.write_scalartimefield(path, 'Vsed', self.Vsed, time)
                    self.r.clean_tps()

        def write_scalartimefield(self, path, name, field, time, nround= 0):
                if time == -1:
                    for t in range(len(self.t)):
                        self.write_scalarvalue(path, name, field[t,0], self.t[t], nround)
                else:
                    self.write_scalarvalue(path, name, field[0,0], time, nround)
                            
        def write_scalarvalue(self, path, name, value, time, nround= 0):
                fi = open(path+'/'+str(time)+'/'+name, 'w')
                fi.write(str(round(value, nround))+"\n")                          
                fi.close()
        
        def write_vectortimefield(self, path, name, field, time, nround = 0):
                if time == -1:
                    for t in range(len(self.t)):
                        self.write_scalarvalue(path, name+'_x', field[t,0], self.t[t], nround)
                        self.write_scalarvalue(path, name+'_y', field[t,1], self.t[t], nround)
                        self.write_scalarvalue(path, name+'_z', field[t,2], self.t[t], nround)
                else:
                    self.write_scalarvalue(path, name+'_x', field[0,0], time, nround)
                    self.write_scalarvalue(path, name+'_y', field[0,1], time, nround)
                    self.write_scalarvalue(path, name+'_z', field[0,2], time, nround) 

        def write_h_tp(self, path, time):
                if self.r.lp_name != '':
                    if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                        self.get_h_tpinterpolation(False, time)
                        self.write_scalartpfield(path, 'h_tp', self.h_tp, time)

        def write_Us_m_tp(self, path, time):
                if self.r.lp_name != '':
                    if (self.r.Us_flag == 'on' or self.r.Us_flag == 'yes' or self.r.Us_flag == True):
                        self.get_Us_m_tpinterpolation(False, time)
                        self.write_scalartpfield(path, 'Us_m_tp', self.Us_m_tp, time)            

        def write_deltaz0_tp(self, path, time):
                if self.r.lp_name != '':
                    if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                        if self.r.tps == {}:
                            self.r.get_transportProperties(False)
                        if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                            self.write_deltac0_tp(path, time)
                        elif self.r.tps['entrainmentModel'] != 'entrainmentOff' or self.r.tps['depositionModel'] != 'depositionOff':
                            self.write_deltah0_tp(path, time)
                            
                        self.r.clean_tps()

        def write_deltah0_tp(self, path, time):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    self.get_deltah0_tpinterpolation(False, time)
                    self.write_scalartpfield(path, 'deltah0_tp', self.deltah0_tp, time)

        def write_deltac0_tp(self, path, time):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    self.get_deltac0_tpinterpolation(False, time)
                    self.write_scalartpfield(path, 'deltac0_tp', self.deltac0_tp, time)

        def write_z_tp(self, path):
                if self.r.lp_name != '':
                    if len(self.z) != 0:
                        self.get_z_tpinterpolation(False)
                        if self.rank == 0:
                            fi = open(path+'/'+'z_tp', 'w')
                            for key in self.z_tp.keys():
                                self.write_tpfield(fi, self.z_tp[key], key, 3)
                            fi.close()

        def write_he_tp(self, path):
                if self.r.lp_name != '':
                    if len(self.he) != 0:
                        self.get_he_tpinterpolation(False)
                        if self.rank == 0:
                            fi = open(path+'/'+'he_tp', 'w')
                            for key in self.he_tp.keys():
                                self.write_tpfield(fi, self.he_tp[key], key, 3)
                            fi.close()

        def write_h(self, path, time):
                if (self.r.h_flag == 'on' or self.r.h_flag == 'yes' or self.r.h_flag == True):
                    self.get_h_interpolation(False, time)
                    self.write_scalarfield(path, 'h', self.h, time)

        def write_pb(self, path, time):
                if (self.r.pb_flag == 'on' or self.r.pb_flag == 'yes' or self.r.pb_flag == True):
                    self.get_pb_interpolation(False, time)
                    self.write_scalarfield(path, 'pb', self.pb, time)

        def write_Cv(self, path, time):
                if (self.r.Cv_flag == 'on' or self.r.Cv_flag == 'yes' or self.r.Cv_flag == True):
                    self.get_Cv_interpolation(False, time)
                    self.write_scalarfield(path, 'Cv', self.Cv, time)

        def write_deltaz0(self, path, time):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    if self.r.tps == {}:
                        self.r.get_transportProperties(False)
                    if self.r.tps['terrainModification'] == 'on' or self.r.tps['terrainModification'] == True:
                        self.write_deltac0(path, time)
                    elif self.r.tps['entrainmentModel'] != 'entrainmentOff' or self.r.tps['depositionModel'] != 'depositionOff':
                        self.write_deltah0(path, time)
                        
                    self.r.clean_tps()

        def write_deltah0(self, path, time):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    self.get_deltah0_interpolation(False, time)
                    self.write_scalarfield(path, 'deltah0', self.deltah0, time, 6) 

        def write_deltac0(self, path, time):
                if (self.r.deltaz0_flag == 'on' or self.r.deltaz0_flag == 'yes' or self.r.deltaz0_flag == True):
                    self.get_deltac0_interpolation(False, time)
                    self.write_scalarfield(path, 'deltac0', self.deltac0, time, 6) 

        def write_Us(self, path, time):
                if (self.r.Us_flag == 'on' or self.r.Us_flag == 'yes' or self.r.Us_flag == True):
                    self.get_Us_interpolation(False, time)
                    self.write_vectorfield(path, 'Us', self.Us, time)
                    self.write_scalarfield(path, 'Us_m', self.Us_m, time)  
                    
        def write_tau(self, path, time):
                if (self.r.tau_flag == 'on' or self.r.tau_flag == 'yes' or self.r.tau_flag == True):
                    self.get_tau_interpolation(False, time)
                    self.write_vectorfield(path, 'tau', self.tau, time)
                    self.write_scalarfield(path, 'tau_m', self.tau_m, time) 

        def write_Sm(self, path, time):
                if (self.Sm_flag == 'on' or self.Sm_flag == 'yes' or self.Sm_flag == True):
                    if self.r.eM_coeffs == {}:
                        self.r.get_entrainment_coeffs(False)    
                    if self.r.eM != 'entrainmentOff':
                        self.get_Sm(False, time)
                        self.write_scalarfield(path, 'Sm', self.Sm, time, 6) 

        def write_rho(self, path, time):
                if (self.rho_flag == 'on' or self.rho_flag == 'yes' or self.rho_flag == True):
                    self.get_rho(False, time)
                    self.write_scalarfield(path, 'rho', self.rho, time) 

        def write_Q(self, path, time):
                if self.r.lp_name != '':
                    if (self.r.Q_flag == 'on' or self.r.Q_flag == 'yes' or self.r.Q_flag == True):
                        self.get_Q(False, time)
                        self.write_scalarlpfield(path, 'Q', self.Q, time, False, 6)
    
        def write_phi2s(self, path, time):
                if self.r.lp_name != '':
                    if (self.r.phi2s_flag == 'on' or self.r.phi2s_flag == 'yes' or self.r.phi2s_flag == True):
                        self.get_phi2s(False, time)
                        self.write_scalarlpfield(path, 'phi2s', self.phi2s, time, False, 6)

        def write_scalarlpfield(self, path, name, field, time, edge = False, nround=3):
                if time == -1:
                    for t in range(len(self.t)):
                        fi = open(path+'/'+str(self.t[t])+'/'+name, 'w')
                        for key in field.keys():
                            if edge == True:
                                self.write_lpfield(fi, field[key]['flux_lp'][t,:], key, nround) 
                            else:
                                self.write_lpfield(fi, field[key][t,:], key, nround)
                        fi.close()
                else:
                    fi = open(path+'/'+str(time)+'/'+name, 'w')
                    for key in field.keys():
                        if edge == True:
                            self.write_lpfield(fi, field[key]['flux_lp'][0,:], key, nround) 
                        else:
                            self.write_lpfield(fi, field[key][0,:], key, nround)
                    fi.close()
                    
        def write_lpfield(self, fi, field, key, nround = 3):
                fi.write(str(key)+"\n")
                fi.write("{"+"\n")
                for i in range(len(field)):
                    fi.write(str(round(field[i], nround))+" ")
                fi.write("\n") 
                fi.write("}"+"\n")
                fi.write("\n")  
                     
        def write_scalartpfield(self, path, name, field, time, nround=3):
                if time == -1:
                    for t in range(len(self.t)):
                        fi = open(path+'/'+str(self.t[t])+'/'+name, 'w')
                        for key in field.keys():
                            self.write_tpfield(fi, field[key][t,:,:], key, nround)                        
                        fi.close()
                else:
                    fi = open(path+'/'+str(time)+'/'+name, 'w')
                    for key in field.keys():
                        self.write_tpfield(fi, field[key], key, nround)                        
                    fi.close()                    

        def write_tpfield(self, fi, field, key, nround):
                fi.write(str(key)+"\n")
                fi.write("{"+"\n")
                for i in range(field.shape[0]):
                    for j in range(field.shape[1]):
                         fi.write(str(round(field[i,j],nround))+" ")
                    fi.write("\n") 
                fi.write("}"+"\n")
                fi.write("\n")  
                
        def write_list(self, path, name, list, nround = 3):
                fi = open(path+'/'+name, 'w')
                for i in range(len(list)):
                    fi.write(str(round(list[i], nround))+" ")
                fi.close()                  

        def write_scalarfield(self, path, name, field, time, nround= 3):
                if time == -1:
                    for t in range(len(self.t)):
                        fi = open(path+'/'+str(self.t[t])+'/'+name, 'w')
                        self.write_field(fi, field[t,:,:], nround)                           
                        fi.close()
                else:
                    fi = open(path+'/'+str(time)+'/'+name, 'w')
                    self.write_field(fi, field, nround)                           
                    fi.close()                    
            
        def write_field(self, fi, field, nround = 3):
                for j in range(len(self.y)):
                    for i in range(len(self.x)):
                         fi.write(str(round(field[i,j], nround))+" ")
                    fi.write("\n")                            

        def write_vectorfield(self, path, name, field, time, nround = 3):
                if time == -1:
                    self.write_scalarfield(path, name+'_x', field[:,:,:,0], time, nround)
                    self.write_scalarfield(path, name+'_y', field[:,:,:,1], time, nround)
                    self.write_scalarfield(path, name+'_z', field[:,:,:,2], time, nround)
                else:
                    self.write_scalarfield(path, name+'_x', field[:,:,0], time, nround)
                    self.write_scalarfield(path, name+'_y', field[:,:,1], time, nround)
                    self.write_scalarfield(path, name+'_z', field[:,:,2], time, nround)                    


########################################  extra functions  #########################################################
                           
def get_list(list, xy, r_xy, c_xy, nTotal, nInter):
    i_inf = 1; i_sup = 1
    ordered_list = []
    for i in range(len(xy)):
        if i_inf == 0:
            i_inf = 1
        if i_sup == len(c_xy)-1:
            i_sup = len(c_xy)-2
        sublist = []
        break_inf = False
        break_sup = False
        while True:
            if break_inf == False:
                if (c_xy[i_inf][0]+r_xy >= xy[i] and
                   c_xy[i_inf-1][0]+r_xy < xy[i]):
                    break_inf = True
                elif c_xy[i_inf-1][0]+r_xy >= xy[i]:
                    i_inf = i_inf-1
                elif c_xy[i_inf][0]+r_xy < xy[i]:
                    i_inf = i_inf+1
                if i_inf == 0:
                    break_inf = True
            
            if break_sup == False:
                if (c_xy[i_sup+1][0]-r_xy > xy[i] and
                   c_xy[i_sup][0]-r_xy <= xy[i]):
                    break_sup = True
                elif c_xy[i_sup][0]-r_xy > xy[i]:
                    i_sup = i_sup-1
                elif c_xy[i_sup+1][0]-r_xy <= xy[i]:
                    i_sup = i_sup+1
                if i_sup == len(c_xy)-1:
                    break_sup = True
                                        
            if break_inf == True and break_sup == True:
                break
        for k in range(i_sup-i_inf+1):
            sublist.append(c_xy[i_inf+k][1])
        sublist.sort()
        ordered_list.append(sublist)
                        
    sub_k = []
    for k in range(nInter):
        sub_k.append(int(nTotal*(k+1)/nInter))
                   
    for i in range(len(ordered_list)):
        k = 0
        j = 0
        sub_list = []
        sub_sub_list = []
        while True:
            if ordered_list[i][j] <= sub_k[k]:
                sub_sub_list.append(ordered_list[i][j])
                j += 1
            else:
                k += 1
                sub_list.append(sub_sub_list)
                sub_sub_list = []
            if j == len(ordered_list[i]):
                sub_list.append(sub_sub_list)
                while True:
                    if len(sub_list) < len(sub_k):
                        sub_list.append([])
                    else:
                        break
                break
                                            
        list.append(sub_list)
       
def surface_smoother(input, output, alpha, n_iterations):
    if n_iterations > 0:
        while True:
            n_iterations -= 1
            for i in range(input.shape[0]-2): 
                for j in range(input.shape[1]-2):
                    output[i+1,j+1] = input[i+1,j+1]*(1-alpha)+alpha*(input[i+2,j+1]+input[i,j+1]+input[i+1,j+2]+input[i+1,j])/4                    
            for i in range(input.shape[0]-2):
                output[i+1,0] = input[i+1,0]*(1-alpha)+alpha*(input[i,0]+input[i+2,0]+2*input[i+1,1])/4      
                output[i+1,input.shape[1]-1] = input[i+1,input.shape[1]-1]*(1-alpha)+alpha*(input[i,input.shape[1]-1]+input[i+2,input.shape[1]-1]+2*input[i+1,input.shape[1]-2])/4      
            for j in range(input.shape[1]-2):
                output[0,j+1] = input[0,j+1]*(1-alpha)+alpha*(input[0,j+2]+input[0,j]+2*input[1,j+1])/4
                output[input.shape[0]-1,j+1] = input[input.shape[0]-1,j+1]*(1-alpha)+alpha*(input[input.shape[0]-1,j+2]+input[input.shape[0]-1,j]+2*input[input.shape[0]-2,j+1])/4
            output[0,0] = input[0,0]
            output[0,input.shape[1]-1] = input[0,input.shape[1]-1]
            output[input.shape[0]-1,0] = input[input.shape[0]-1,0]
            output[input.shape[0]-1,input.shape[1]-1] = input[input.shape[0]-1,input.shape[1]-1]
            if n_iterations == 0:
                break
            else:
                input = output
    else:
        for i in range(input.shape[0]): 
            for j in range(input.shape[1]):
                output[i,j] = round(input[i,j],6)

def curve_smoother(input,  alpha, n_iterations):
    output = np.zeros((input.shape[0],input.shape[1]))
    if n_iterations > 0:
        while True:
            n_iterations -= 1
            for i in range(output.shape[0]):
                for j in range(output.shape[1]-4):
                        output[i,j+2] = input[i,j+2]*(1-alpha)+alpha*(input[i,j]+input[i,j+1]+input[i,j+3]+input[i,j+4])/4
                output[i,0] = input[i,0]
                output[i,1] = input[i,1]*(1-alpha)+alpha*(input[i,0]+input[i,2])/2
                output[i,-2] = input[i,-2]*(1-alpha)+alpha*(input[i,-3]+input[i,-1])/2
                output[i,-1] = input[i,-1]
                
            if n_iterations == 0:
                break
            else:
                input = output
    else:
        output = input
    
    return output

def face_in(List,v):
    is_in = True
    for i in range(len(List)-1):
        if check_in(List[i],List[i+1],v) != True:
            is_in = False
            break
    if check_in(List[-1],List[0],v) != True:
        is_in = False                
    return is_in

def get_closer_point(m,P0,P1,P2): #determines which point, P1 or P2, is the closest to P0
    dp1 = abs(m*(P1.x-P0.x)+P0.y-P1.y)/math.sqrt(m**2+1)
    dp2 = abs(m*(P2.x-P0.x)+P0.y-P2.y)/math.sqrt(m**2+1)
    if dp1 > dp2:
        return 2
    else:
        return 1

def check_in(v1,v2,v_check):
    quadrant = get_quadrant(v1,v2)
    if v2.x != v1.x:
        y_check = (v2.y-v1.y)/(v2.x-v1.x)*(v_check.x-v1.x)+v1.y
        if quadrant ==1 or quadrant ==4:
            if v_check.y <= y_check:
                return True
            else:
                return False
        else:
            if v_check.y >= y_check:
                return True
            else:
                return False 
    else:
        if quadrant ==1:
            if v_check.x >= v1.x:
                return True
            else:
                return False
        else:
            if v_check.x <= v1.x:
                return True
            else:
                return False   

def choose_sense(n, index_i, index_f): #n is the number of points in the face.
    if index_f > index_i:
        watch_sense = index_f-index_i
        anti_watch_sense = n-(index_f-index_i)
    else:
        watch_sense = n-(index_i-index_f)
        anti_watch_sense = index_i-index_f
    if watch_sense <= anti_watch_sense:
        return 1
    else:
        return -1

def get_flux_tP(input_flux,list_edges,t,edges_sense):
    edge_field = input_flux[str(t)].iF.field
    edge_field_u = input_flux[str(t)].iF.u
    flux = np.zeros((1,len(list_edges)))
    if edge_field_u == 'uniform':
        flux = np.ones((1,len(list_edges)))*edge_field[0]
    else:
        for i in range(len(list_edges)):
            flux[0,i] = edge_field[list_edges[i]]*edges_sense[i]
    return flux

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def get_point_t(point_i, point_f, dx):
    point_m = (point_i+point_f)*0.5
    quadrant = get_quadrant(point_i,point_f)
    if point_f.x != point_i.x:
        m = (point_f.y-point_i.y)/(point_f.x-point_i.x)
        if m != 0:
            if quadrant == 1 or quadrant == 2:
                dist = 0.001*dx
            else:
                dist = -0.001*dx
            point_t = vector(1,-1/m,0)/math.sqrt(1+(1/m)**2)*dist+point_m
        else:
            if quadrant == 1 or quadrant == 4:
                dist = 0.001*dx
            else:
                dist = -0.001*dx
            point_t = vector(0,-1,0)*dist+point_m
    else:
        if quadrant == 1 or quadrant == 2:
            dist = 0.001*dx
        else:
            dist = -0.001*dx
        point_t = vector(1,0,0)*dist+point_m
        
    return point_t
                      
def get_quadrant(v1,v2):
    if v2.x >= v1.x and v2.y >= v1.y:
        quadrant = 1
    elif v2.x < v1.x and v2.y >= v1.y:
        quadrant = 2
    elif v2.x < v1.x and v2.y < v1.y:
        quadrant = 3
    elif v2.x >= v1.x and v2.y < v1.y:
        quadrant = 4
    return quadrant

def sign(x):
    if float(x) > 0:
        return 1
    elif float(x) < 0:
        return -1
    elif float(x) == 0:
        return 0


########################################  main function  #########################################################

def clean_dir(path, dir, rank = 0):
    if rank == 0:
        listdir = os.listdir(path)
        if dir in listdir:
            path += '/'+dir
            for filename in os.listdir(path):
                file_path = os.path.join(path, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print('Failed to delete %s. Reason: %s' % (file_path, e))
        
            os.chdir(path)
            os.rmdir(path)

def make_dir(path, dir, rank = 0):
    if rank == 0:
        listdir = os.listdir(path)
        if not(dir in listdir):
            os.mkdir(path+'/'+dir)

def date_format(x):
    if x<10:
        return '0'+str(x)
    else:
        return str(x)

def date(t):
    yr = date_format(t[0])
    mo = date_format(t[1])
    dy = date_format(t[2])
    hr = date_format(t[3])
    mi = date_format(t[4])
    sc = date_format(t[5])
    
    return dy+'/'+mo+'/'+yr+'  '+hr+':'+mi+':'+sc    
    
def write_summary(path, lp, tp, dx, dy, rx, ry, alpha, niter, alphafield, niterfield, dist, ntp, rank):
    if rank == 0:
        import time
        t = time.gmtime()
    
        filename = 'Summary.dat'
        fi = open(path+'/'+filename, 'w')
        dt = date(t)
        
        fi.write("Summary file written on "+ dt +"\n")
        fi.write("\n")
        fi.write("With dx = "+ str(dx) + " and dy = "+ str(dy) + "\n")
        fi.write("With rx = "+ str(rx) + " and ry = "+ str(ry) + "\n")
        fi.write("With alpha = "+ str(alpha) + " and n° iterations = "+ str(niter) +"\n")
        fi.write("With alpha field = "+ str(alphafield) + " and n° iterations field = "+ str(niterfield) +"\n")
        fi.write("With transversal profiles width = "+ str(dist) + " and n° of points = "+ str(ntp) + "\n")
        fi.write("\n")
        fi.close()    

def main(argv):
    from mpi4py import MPI
    import time
    import shutil
    
    parser = argparse.ArgumentParser(description='Reading the outputs of a debrisfaSavageHutterFoam simulation')

    parser.add_argument('-path', type=str, help='input simulationPath', default='')
    parser.add_argument('-o', type=str, help='output folder', default='Results')
    parser.add_argument('-h_flag', help='Activate to read h field', action="store_true")
    parser.add_argument('-pb_flag', help='Activate to read pb field', action="store_true")
    parser.add_argument('-Cv_flag', help='Activate to read Cv field', action="store_true")
    parser.add_argument('-deltaz0_flag', help='Activate to read deltaz0 field', action="store_true")
    parser.add_argument('-Us_flag', help='Activate to read Us field', action="store_true")
    parser.add_argument('-tau_flag', help='Activate to read tau field', action="store_true")
    parser.add_argument('-phi2s_flag', help='Activate to read phi2s field', action="store_true")
    parser.add_argument('-Q_flag', help='Activate to read Q field', action="store_true")
    parser.add_argument('-c_flag', help='Activate to read c field', action="store_true")
    parser.add_argument('-n_flag', help='Activate to read n field', action="store_true")
    parser.add_argument('-he_flag', help='Activate to read he field', action="store_true")
    parser.add_argument('-A_flag', help='Activate to read A field', action="store_true")
        
    parser.add_argument('-lp', type=str, help='longitudinal profile', default='')
    parser.add_argument('-tp', type=str, help='transversal profile', default='')
    parser.add_argument('-dx', type=float, help='distance between elements in the x direction')
    parser.add_argument('-dy', type=float, help='distance between elements in the y direction')
    parser.add_argument('-rx', type=float, help='distance max used for the square inverse distance method')
    parser.add_argument('-ry', type=float, help='distance max used for the square inverse distance method')
    parser.add_argument('-alpha', type=float, help='alpha value used to smooth the interpolated z field', default=0)
    parser.add_argument('-niter', type=int, help='number of iterations used in the z field smoothing algorithm', default=0)
    parser.add_argument('-alphafield', type=float, help='alpha value used to smooth the interpolated fields', default=0)
    parser.add_argument('-niterfield', type=int, help='number of iterations used in the fields smoothing algorithm', default=0)
    parser.add_argument('-dist', type=float, help='wide of the transversal profiels', default=0)
    parser.add_argument('-ntp', type=int, help='number of points in every transversal profile', default=100)
    parser.add_argument('-Sm_flag', help='Activate Sm calculation', action="store_true")
    parser.add_argument('-rho_flag', help='Activate rho calculation', action="store_true")
    parser.add_argument('-rcg_flag', help='Activate rcg calculation', action="store_true")
    parser.add_argument('-M_flag', help='Activate M calculation', action="store_true")
    parser.add_argument('-V_flag', help='Activate V calculation', action="store_true")
    parser.add_argument('-Vsed_flag', help='Activate Vsed calculation', action="store_true")
    
    args = parser.parse_args()

    input_path = args.path
    
    h  = args.h_flag
    pb = args.pb_flag
    Cv = args.Cv_flag
    deltaz0 = args.deltaz0_flag
    Us  = args.Us_flag
    tau = args.tau_flag
    phi2s = args.phi2s_flag
    Q  = args.Q_flag
    c  = args.c_flag
    n  = args.n_flag
    he = args.he_flag
    A  = args.A_flag
        
    lp = args.lp
    tp = args.tp
    dx = args.dx
    dy = args.dy
    rx = args.rx
    ry = args.ry
    alpha = args.alpha
    niter = args.niter
    alphafield = args.alphafield
    niterfield = args.niterfield
    dist = args.dist
    ntp = args.ntp
    Sm  = args.Sm_flag
    rho = args.rho_flag
    rcg = args.rcg_flag
    M = args.M_flag
    V = args.V_flag
    Vsed = args.Vsed_flag
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print('Rank = ' + str(rank)+'  started at = ' + str(date(time.gmtime())))

    if input_path == '':
        input_path = os.getcwd()

    output_path = input_path+'/'+args.o

    clean_dir(input_path, args.o, rank)
    make_dir(input_path, args.o, rank)

    if rank == 0:
        shutil.copy(input_path+'/constant/transportProperties', output_path+'/transportProperties')
        #Be careful with the following line, because if the file does not exist the code will fail //AG
        shutil.copy(input_path+'/log.debrisfaSavageHutterFoam', output_path+'/log.debrisfaSavageHutterFoam')
    write_summary(output_path,lp,tp,dx,dy,rx,ry,alpha,niter,alphafield,niterfield,dist,ntp,rank)

    if rank == 0:
        r = runCase(input_path, lp, tp, h, pb, Cv, deltaz0, Us, tau, phi2s, Q, c, n, he, A)
        r.report_Case(output_path+'/'+'Summary.dat')
    else:
        r = None
        
    r = comm.bcast(r, root=0)
    
    o = outPut(r,dx,dy,rx,ry,alpha,niter,alphafield,niterfield,dist,ntp,rank,Sm,rho,rcg,M,Vsed,V)
    o.write_output(output_path, size)

    print('Rank = ' + str(o.rank)+'  ended at = ' + str(date(time.gmtime())))

if __name__ == "__main__":
  main(sys.argv[1:]) 