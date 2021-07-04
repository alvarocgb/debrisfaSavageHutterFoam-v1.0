#!/usr/bin/env python3

'''
Description
    Code used to reconstruct a parallel simulation of the debrisfaSavageHutterFoam. 
    With this code there is no need to reconstruct some fields, like the ones generated 
    in a simulation where the mesh is modified. It also can be run in parallel, something that reconstructPar cannot.

    The code is not capable of recognising when an internalfield or a boundarypatchfield is uniform.
    It will always write a list of values (even if they are all the same) and say that it is a nonuniform list.
    This problem can be easily solved calculating the max() and min() in every list before writing //AG
Author
    Álvaro González Bilbao alvaro.gonzalez.bilbao@gmail.com

@version: 1.10 (21/06/21)
'''

import sys
import os
import argparse
import shutil

##########################################  Definitions  ##############################################################

def header(objectName, className, locationName):
    return ("/*--------------------------------*- C++ -*----------------------------------*\\\n"
            "| =========                 |                                                 |\n"
            "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
            "|  \\\\    /   O peration     | Version:  v1812                                 |\n"
            "|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n"
            "|    \\\\/     M anipulation  |                                                 |\n"
            "\*---------------------------------------------------------------------------*/\n"
            "FoamFile\n"
            "{}\n"
            "    version     2.0;\n"
            "    format      ascii;\n"
            "    class       {};\n"
            "    location    \"{}\";\n"
            "    object      {};\n"
            "{}\n"
            "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
            ).format("{", className, locationName, objectName, "}");

def footer():
    return ("\n"
            "// ************************************************************************* //\n")
            
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

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False   

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

def correct_EdgeLabels(list):
    new_list = []
    for i in range(len(list)):
        if (list[i][0] == 'edgeLabels' and len(list[i])>2):
            if (list[i][1] == 'List<label>' and list[i][2][0] == '0'):
                new_list.append([list[i][0], list[i][2]])            
            else:
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

def create_files_list(path):
    listdir = os.listdir(path)
    k = 0
    while True:
        if listdir[k]=='polyMesh' or listdir[k]=='uniform' or listdir[k]=='meshPhi':
            listdir.remove(listdir[k])
        else:
            k += 1
        if k == len(listdir):
            break
    return listdir

def read_proc_number_faces(path, n_proc, output):
    for p in range(n_proc):
        list = []
        get_input_number(path+'/processor'+str(p)+'/constant/faMesh/faceLabels',list) 
        output['processor'+str(p)] = int(list[0][0])

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
                                
        def __eq__(self, other):
                if self.x == other.x and self.y == other.y and self.z == other.z:
                    return True
                else:
                    return False

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
                        
def read_areaField_proc(path, time, name, type, output, n_processors, proc_fB, proc_nF):
    dict = {}
    for p in range(n_processors):
        list = []
        get_input(path+'/processor'+str(p)+'/'+str(time)+'/'+name,list)
        list = correct_bFEdgeFields(list,type)
        field = {}
        read_dictionary(list,field)
        dimension = dimensions(field['dimensions'])
        boundaryfield = boundaryField(field['boundaryField'],type,proc_fB['processor'+str(p)])
        internalfield = internalField(field,type,proc_nF['processor'+str(p)])
        dict['processor'+str(p)] = areaField(type,internalfield,boundaryfield,dimension)
    output[str(time)] = dict

def correct_bFEdgeFields(list,type): #This is incredibly inefficient. It must be improved in a future version. //AG
    # Since this is useful only for simulations generated by the newest versions of OpenFOAM, it might be useful to get the version inside the code
    # and apply the correction only if it is completely necessary. //AG
    new_list = []
    key = 'List<'+type+'>'
    for i in range(len(list)):
        if (list[i][0] == 'value' and len(list[i])>3):
            if list[i][2] == key:
                pop = list[i].pop(2)
        new_list.append(list[i])
    return new_list

def read_edgeField_proc(path, time, name, type, output, n_processors, proc_fB, proc_nE):
    dict = {}
    for p in range(n_processors):
        list = []
        get_input(path+'/processor'+str(p)+'/'+str(time)+'/'+name,list)
        list = correct_bFEdgeFields(list,type)
        field = {}
        read_dictionary(list,field)
        dimension = dimensions(field['dimensions'])
        boundaryfield = boundaryField(field['boundaryField'],type,proc_fB['processor'+str(p)])
        internalfield = internalField(field,type,proc_nE['processor'+str(p)]['internal'])
        dict['processor'+str(p)] = areaField(type,internalfield,boundaryfield,dimension)
    output[str(time)] = dict

def read_edgeField(path,time,name,type,output,number_edges,faBoundary):
    list = []
    get_input(path+'/'+str(time)+'/'+name,list)
    field = {}
    read_dictionary(list,field)
    dimension = dimensions(field['dimensions'])
    boundaryfield = boundaryField(field['boundaryField'],type,faBoundary)
    internalfield = internalField(field,type,number_edges)
    output[str(time)] = areaField(type,internalfield,boundaryfield,dimension)

def read_ec(path,output, number_edges,faBoundary):
    read_edgeField(path, 0, 'ec','vector',output, number_edges,faBoundary)

def read_ec_proc(path,output,n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, 0, 'ec','vector',output, n_processors, proc_fB, proc_nE)
       
def read_h_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'h','scalar',output, n_processors, proc_fB, proc_nF)

def read_pb_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'pb','scalar',output, n_processors, proc_fB, proc_nF)
                
def read_Cv_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'Cv','scalar',output, n_processors, proc_fB, proc_nF)
                
def read_deltac0_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'deltac0','scalar',output, n_processors, proc_fB, proc_nF)

def read_deltah0_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'deltah0','scalar',output, n_processors, proc_fB, proc_nF)

def read_Us_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'Us','vector',output, n_processors, proc_fB, proc_nF)

def read_tau_proc(path, time, output, n_processors, proc_fB, proc_nF):
    read_areaField_proc(path, time, 'tau','vector',output, n_processors, proc_fB, proc_nF)

def read_Q_proc(path, time, output, n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, time, 'Q','scalar',output, n_processors, proc_fB, proc_nE)

def read_phi2s_proc(path, time, output, n_processors, proc_fB, proc_nE):
    read_edgeField_proc(path, time, 'phi2s','scalar',output, n_processors, proc_fB, proc_nE)

##########################################  runCase class  ##############################################################

class runCase:
        def __init__(self,path, h='on',Cv='on',deltaz0='on',Us='on',Q='on',pb='on',tau='on', phi2s = 'on', rank = 0):
                self.t= []
                self.h_proc  = {}               
                self.Cv_proc = {}
                self.pb_proc = {}
                self.deltah0_proc = {}
                self.deltac0_proc = {}
                self.Us_proc  = {}                               
                self.tau_proc = {}
                self.Q_proc = {}
                self.phi2s_proc = {}
                
                self.h  = {}               
                self.Cv = {}
                self.pb = {}
                self.deltah0 = {}
                self.deltac0 = {}
                self.Us  = {}               
                self.tau = {}
                self.Q = {}
                self.phi2s = {}
                
                self.p  = path
                self.np = number_processors(self.p)
                self.fB = {}  #faBoundary
                self.nF = get_number_faces(self.p)
                self.nE = {}  #number Edges             
                self.ec = {}
                self.proc_fB = {}  #processors faBoundary
                self.proc_nF = {}  #processors number of faces
                self.proc_nE = {}  #processors number of edges
                self.proc_ec = {}
                self.proc_faceaddr = {}
                self.proc_edgeaddr = {}
                self.t_size = []  # time list used when runnning in parallel                
                self.rank = rank  #processor number
               
                self.h_flag  = h
                self.Cv_flag = Cv
                self.pb_flag = pb
                self.deltaz0_flag = deltaz0
                self.deltac0_flag = False
                self.deltah0_flag = False
                self.Us_flag  = Us
                self.tau_flag = tau
                self.Q_flag = Q              
                self.phi2s_flag = phi2s
                
                self.get_case()

        def get_case(self):
                self.get_time()
                self.get_files_names()
                self.get_proc_nF()
                self.get_faBoundary()
                self.get_nEdges()
                self.get_proc_faBoundary()
                self.get_proc_nE()
                self.get_proc_faceaddr()
                self.get_proc_edgeaddr()
                self.get_ec()
                self.get_proc_ec()
                self.order_proc_faBoundary()

        def get_time(self):
                if self.rank == 0:
                    print ('Reading times ...')
                create_time(self.p+'/processor0', self.t)

        def get_files_names(self):
                if self.rank == 0:
                    print ('Reading files names ...')
                files_list = create_files_list(self.p+'/processor0'+'/'+str(self.t[-1]))

                if (self.h_flag == 'on' or  self.h_flag == True) and not('h' in files_list):
                    self.h_flag = 'off'
                if (self.Cv_flag == 'on' or  self.Cv_flag == True) and not('Cv' in files_list):
                    self.Cv_flag = 'off'
                if (self.pb_flag == 'on' or  self.pb_flag == True) and not('pb' in files_list):
                    self.pb_flag = 'off'
                if (self.Us_flag == 'on' or  self.Us_flag == True) and not('Us' in files_list):
                    self.Us_flag = 'off'
                if (self.tau_flag == 'on' or  self.tau_flag == True) and not('tau' in files_list):
                    self.tau_flag = 'off'
                if (self.Q_flag == 'on' or  self.Q_flag == True) and not('Q' in files_list):
                    self.Q_flag = 'off'              
                if (self.phi2s_flag == 'on' or  self.phi2s_flag == True) and not('phi2s' in files_list):
                    self.phi2s_flag = 'off'

                if (self.deltaz0_flag == 'on' or  self.deltaz0_flag == True) and not('deltac0' in files_list or 'deltah0' in files_list):
                    self.deltaz0_flag = 'off'
                elif (self.deltaz0_flag == 'on' or  self.deltaz0_flag == True) and ('deltac0' in files_list):
                    self.deltac0_flag = 'on'
                elif (self.deltaz0_flag == 'on' or  self.deltaz0_flag == True) and ('deltah0' in files_list):
                    self.deltah0_flag = 'on' 

        def get_proc_nF(self):
                read_proc_number_faces(self.p, self.np, self.proc_nF)

        def get_proc_nE(self):
                read_proc_number_edges(self.p, self.np, self.proc_nE, self.proc_fB)

        def get_faBoundary(self):
                if self.rank == 0:
                    print ('Reading faBoundary ...')
                read_faBoundary(self.p,self.fB)

        def get_nEdges(self):
                self.nE['internal'] = get_number_edges(self.p)
                sum = 0
                for key in self.fB['boundaries']:
                    sum += len(self.fB['boundaries'][key])
                self.nE['Total'] = self.nE['internal'] + sum

        def get_proc_faBoundary(self):
                if self.rank == 0:
                    print ('Reading processors faBoundary ...')            
                read_proc_faBoundary(self.p, self.np, self.proc_fB)

        def get_proc_faceaddr(self):
                if self.rank == 0:
                    print ('Reading processors faceaddressing ...')            
                read_proc_faceaddr(self.p, self.np, self.proc_faceaddr)

        def get_proc_edgeaddr(self):
                if self.rank == 0:
                    print ('Reading processors edgeaddressing ...')            
                read_proc_edgeaddr(self.p, self.np, self.proc_edgeaddr)

        def create_t_size(self, size, withZero):
                for i in range(size):
                    sub_list = []
                    for j in range(int(len(self.t)/size)):
                        sub_list.append(self.t[j*size+i])
                    self.t_size.append(sub_list)
                if withZero == False:
                    self.t_size[0].pop(0)
                    for i in range(len(self.t)-int(len(self.t)/size)*size):
                        self.t_size[i].append(self.t[int(len(self.t)/size)*size+i])
                else:
                    for i in range(len(self.t)-int(len(self.t)/size)*size):
                        self.t_size[i+1].append(self.t[int(len(self.t)/size)*size+i])                

        def get_ec(self):
                if self.rank == 0:
                    print('Reading ec field ...')
                read_ec(self.p,self.ec, self.nE['internal'], self.fB)

        def get_proc_ec(self):
                if self.rank == 0:
                    print('Reading processor ec field ...')            
                read_ec_proc(self.p,self.proc_ec,self.np, self.proc_fB, self.proc_nE)

        def get_h(self, report = True, time = -1):                
                if (self.h_flag == 'on' or self.h_flag == 'yes' or self.h_flag == True):
                    if report == True:
                        print('Reading h field ...')
                    if time == -1:
                        for t in self.t:
                            read_h_proc(self.p, t, self.h_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_h_proc(self.p, time, self.h_proc, self.np, self.proc_fB, self.proc_nF)

        def get_Cv(self, report = True, time = -1):                
                if (self.Cv_flag == 'on' or self.Cv_flag == 'yes' or self.Cv_flag == True):
                    if report == True:
                        print('Reading Cv field ...')
                    if time == -1:
                        for t in self.t:
                            read_Cv_proc(self.p, t, self.Cv_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_Cv_proc(self.p, time, self.Cv_proc, self.np, self.proc_fB, self.proc_nF)

        def get_pb(self, report = True, time = -1):                
                if (self.pb_flag == 'on' or self.pb_flag == 'yes' or self.pb_flag == True):
                    if report == True:
                        print('Reading pb field ...')
                    if time == -1:
                        for t in self.t:
                            read_pb_proc(self.p, t, self.pb_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_pb_proc(self.p, time, self.pb_proc, self.np, self.proc_fB, self.proc_nF)

        def get_deltaz0(self, report = True, time = -1):                
                if (self.deltac0_flag == 'on' or self.deltac0_flag == 'yes' or self.deltac0_flag == True):
                    if report == True:
                        print('Reading deltac0 field ...')
                    if time == -1:
                        for t in self.t:
                            read_deltac0_proc(self.p, t, self.deltac0_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_deltac0_proc(self.p, time, self.deltac0_proc, self.np, self.proc_fB, self.proc_nF)
                elif (self.deltah0_flag == 'on' or self.deltah0_flag == 'yes' or self.deltah0_flag == True):
                    if report == True:
                        print('Reading deltah0 field ...')
                    if time == -1:
                        for t in self.t:
                            read_deltah0_proc(self.p, t, self.deltah0_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_deltah0_proc(self.p, time, self.deltah0_proc, self.np, self.proc_fB, self.proc_nF)

        def get_Us(self, report = True, time = -1):                
                if (self.Us_flag == 'on' or self.Us_flag == 'yes' or self.Us_flag == True):
                    if report == True:
                        print('Reading Us field ...')
                    if time == -1:
                        for t in self.t:
                            read_Us_proc(self.p, t, self.Us_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_Us_proc(self.p, time, self.Us_proc, self.np, self.proc_fB, self.proc_nF)

        def get_tau(self, report = True, time = -1):                
                if (self.tau_flag == 'on' or self.tau_flag == 'yes' or self.tau_flag == True):
                    if report == True:
                        print('Reading tau field ...')
                    if time == -1:
                        for t in self.t:
                            read_tau_proc(self.p, t, self.tau_proc, self.np, self.proc_fB, self.proc_nF)
                    else:
                        read_tau_proc(self.p, time, self.tau_proc, self.np, self.proc_fB, self.proc_nF)

        def get_Q(self, report = True, time = -1):                
                if (self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True):
                    if report == True:
                        print('Reading Q field ...')
                    if time == -1:
                        for t in self.t:
                            read_Q_proc(self.p, t, self.Q_proc, self.np, self.proc_fB, self.proc_nE)
                    else:
                        read_Q_proc(self.p, time, self.Q_proc, self.np, self.proc_fB, self.proc_nE)

        def get_phi2s(self, report = True, time = -1):                
                if (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True):
                    if report == True:
                        print('Reading phi2s field ...')
                    if time == -1:
                        for t in self.t:
                            read_phi2s_proc(self.p, t, self.phi2s_proc, self.np, self.proc_fB, self.proc_nE)
                    else:
                        read_phi2s_proc(self.p, time, self.phi2s_proc, self.np, self.proc_fB, self.proc_nE)

        def create_h(self, report = True, time = -1):                
                if (self.h_flag == 'on' or self.h_flag == 'yes' or self.h_flag == True):
                    if report == True:
                        print('Creating h field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.h_proc, self.get_h, self.h, 'scalar')
                    else:
                        self.create_areaField(time, self.h_proc, self.get_h, self.h, 'scalar')

        def create_Cv(self, report = True, time = -1):                
                if (self.Cv_flag == 'on' or self.Cv_flag == 'yes' or self.Cv_flag == True):
                    if report == True:
                        print('Creating Cv field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.Cv_proc, self.get_Cv, self.Cv, 'scalar')
                    else:
                        self.create_areaField(time, self.Cv_proc, self.get_Cv, self.Cv, 'scalar')

        def create_pb(self, report = True, time = -1):                
                if (self.pb_flag == 'on' or self.pb_flag == 'yes' or self.pb_flag == True):
                    if report == True:
                        print('Creating pb field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.pb_proc, self.get_pb, self.pb, 'scalar')
                    else:
                        self.create_areaField(time, self.pb_proc, self.get_pb, self.pb, 'scalar')

        def create_deltaz0(self, report = True, time = -1):                
                if (self.deltac0_flag == 'on' or self.deltac0_flag == 'yes' or self.deltac0_flag == True):
                    if report == True:
                        print('Creating deltac0 field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.deltac0_proc, self.get_deltaz0, self.deltac0, 'scalar')
                    else:
                        self.create_areaField(time, self.deltac0_proc, self.get_deltaz0, self.deltac0, 'scalar')                
                elif (self.deltah0_flag == 'on' or self.deltah0_flag == 'yes' or self.deltah0_flag == True):
                    if report == True:
                        print('Creating deltah0 field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.deltah0_proc, self.get_deltaz0, self.deltah0, 'scalar')
                    else:
                        self.create_areaField(time, self.deltah0_proc, self.get_deltaz0, self.deltah0, 'scalar') 

        def create_Us(self, report = True, time = -1):                
                if (self.Us_flag == 'on' or self.Us_flag == 'yes' or self.Us_flag == True):
                    if report == True:
                        print('Creating Us field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.Us_proc, self.get_Us, self.Us, 'vector')
                    else:
                        self.create_areaField(time, self.Us_proc, self.get_Us, self.Us, 'vector')

        def create_tau(self, report = True, time = -1):                
                if (self.tau_flag == 'on' or self.tau_flag == 'yes' or self.tau_flag == True):
                    if report == True:
                        print('Creating tau field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_areaField(t, self.tau_proc, self.get_tau, self.tau, 'vector')
                    else:
                        self.create_areaField(time, self.tau_proc, self.get_tau, self.tau, 'vector')

        def create_Q(self, report = True, time = -1):                
                if (self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True):
                    if report == True:
                        print('Creating Q field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_edgeField(t, self.Q_proc, self.get_Q, self.Q, 'scalar')
                    else:
                        self.create_edgeField(time, self.Q_proc, self.get_Q, self.Q, 'scalar')

        def create_phi2s(self, report = True, time = -1):                
                if (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True):
                    if report == True:
                        print('Creating phi2s field ...')
                    if time == -1:
                        for t in self.t:
                            self.create_edgeField(t, self.phi2s_proc, self.get_phi2s, self.phi2s, 'scalar')
                    else:
                        self.create_edgeField(time, self.phi2s_proc, self.get_phi2s, self.phi2s, 'scalar')

        def create_areaField(self, time, field_proc, field_func, field, type):
                if (str(time) in field_proc.keys()) == False:
                    field_func(False, time)
                list = [-9999]*self.nF
                for p in range(self.np):
                    for i in range(self.proc_nF['processor'+str(p)]):
                        if field_proc[str(time)]['processor'+str(p)].iF.u == 'uniform':
                            list[self.proc_faceaddr['processor'+str(p)][i]] = field_proc[str(time)]['processor'+str(p)].iF.field[0]
                        else:
                            list[self.proc_faceaddr['processor'+str(p)][i]] = field_proc[str(time)]['processor'+str(p)].iF.field[i]
                dict_iF = {}                
                dict_iF['field'] = list
                dict_iF['u'] = 'nonuniform'
                dimension = field_proc[str(time)]['processor0'].dim
                internalfield = internalField(dict_iF,type,0,False)
                
                dict_bF = {}
                for key in self.fB['boundaries'].keys():
                    dict_bF[key] = {}
                for proc_key in field_proc[str(time)].keys():
                    for patch_key in field_proc[str(time)][proc_key].bF.d.keys():
                        if field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'zeroGradient':
                            dict_bF[patch_key]['type'] = 'zeroGradient'
                        elif field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'fixedValue' or field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'calculated':
                            dict_bF[patch_key]['type'] = field_proc[str(time)][proc_key].bF.d[patch_key]['type']
                            value_dict = {}
                            value_dict['uniform'] = 'nonuniform'
                            value_dict['List'] = [-9999]*len(self.fB['boundaries'][patch_key])
                            dict_bF[patch_key]['value'] = value_dict
                            
                for proc_key in field_proc[str(time)].keys():
                    for patch_key in field_proc[str(time)][proc_key].bF.d.keys():
                        if field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'fixedValue' or field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'calculated':
                            proc_fB = self.proc_fB[proc_key]['boundaries'][patch_key]
                            fB = self.fB['boundaries'][patch_key]
                            values = field_proc[str(time)][proc_key].bF.d[patch_key]['value']['List']
                            if field_proc[str(time)][proc_key].bF.d[patch_key]['value']['uniform'] == 'yes':
                                for i in range(len(self.proc_fB[proc_key]['boundaries'][patch_key])):
                                    index = fB.index(self.proc_edgeaddr[proc_key][proc_fB[i]])
                                    dict_bF[patch_key]['value']['List'][index] = values[0]
                            else:
                                for i in range(len(values)):
                                    index = fB.index(self.proc_edgeaddr[proc_key][proc_fB[i]])
                                    dict_bF[patch_key]['value']['List'][index] = values[i]                                
                
                boundaryfield = boundaryField(dict_bF,type,0,False)                
                field[str(time)] = areaField(type,internalfield,boundaryfield,dimension)

        def create_edgeField(self, time, field_proc, field_func, field, type):
                if (str(time) in field_proc.keys()) == False:
                    field_func(False, time)
                list = [-9999]*self.nE['Total']

                for p in range(self.np):
                    if field_proc[str(time)]['processor'+str(p)].iF.u == 'uniform':
                        for i in range(self.proc_nE['processor'+str(p)]['internal']):
                            list[self.proc_edgeaddr['processor'+str(p)][i]] = field_proc[str(time)]['processor'+str(p)].iF.field[0]
                    else:
                        for i in range(self.proc_nE['processor'+str(p)]['internal']):
                            list[self.proc_edgeaddr['processor'+str(p)][i]] = field_proc[str(time)]['processor'+str(p)].iF.field[i]
                    for key in self.proc_fB['processor'+str(p)]['boundaries'].keys():
                        if field_proc[str(time)]['processor'+str(p)].bF.d[key]['value']['uniform'] == 'yes':
                            for i in range(len(self.proc_fB['processor'+str(p)]['boundaries'][key])):
                                index = self.proc_fB['processor'+str(p)]['boundaries'][key][i]
                                list[self.proc_edgeaddr['processor'+str(p)][index]] = field_proc[str(time)]['processor'+str(p)].bF.d[key]['value']['List'][0]
                        else:
                            for i in range(len(self.proc_fB['processor'+str(p)]['boundaries'][key])):
                                index = self.proc_fB['processor'+str(p)]['boundaries'][key][i]
                                list[self.proc_edgeaddr['processor'+str(p)][index]] = field_proc[str(time)]['processor'+str(p)].bF.d[key]['value']['List'][i]
                dict_iF = {}
                dict_iF['field'] = list[:self.nE['internal']]
                dict_iF['u'] = 'nonuniform'
                dimension = field_proc[str(time)]['processor0'].dim
                internalfield = internalField(dict_iF,type,0,False)
                
                dict_bF = {}
                for key in self.fB['boundaries'].keys():
                    dict_bF[key] = {}
                for proc_key in field_proc[str(time)].keys():
                    for patch_key in field_proc[str(time)][proc_key].bF.d.keys():
                        if field_proc[str(time)][proc_key].bF.d[patch_key]['type'] == 'calculated':
                            dict_bF[patch_key]['type'] = 'calculated'
                            value_dict = {}
                            value_dict['uniform'] = 'nonuniform'
                            value_dict['List'] = [-9999]*len(self.fB['boundaries'][patch_key])
                            dict_bF[patch_key]['value'] = value_dict
                            
                for patch_key in dict_bF.keys():
                    for i in range(len(self.fB['boundaries'][patch_key])):
                        dict_bF[patch_key]['value']['List'][i] = list[self.fB['boundaries'][patch_key][i]]

                boundaryfield = boundaryField(dict_bF,type,0,False)
                field[str(time)] = areaField(type,internalfield,boundaryfield,dimension)

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
                        
        def write_output(self,path,size, time = -1, withZero = False):                    
                if self.rank == 0:
                    print('Writing output files ...')
                self.create_dir(path, size, time, withZero)
                self.write_fields(path, size, time) 
             
        def create_dir(self, path, size, time, withZero):
                if size == 1:
                    if time == -1:
                        for t in self.t:
                            os.mkdir(path+'/'+str(t))
                    else:
                        os.mkdir(path+'/'+str(time))
                else:
                    if len(self.t_size) == 0:
                        self.create_t_size(size, withZero)
                        
                    if time == -1:
                        for t in self.t_size[self.rank]:
                            os.mkdir(path+'/'+str(t))

        def write_fields(self, path, size, time):
                import time as timex
                if time == -1:
                    if size == 1:
                        for t in self.t:
                            print('   Reconstructing fields for time = '+ str(t)+ ' at time = ' + date(timex.gmtime()))
                            self.write_areaFields(path,t)
                            self.write_edgeFields(path,t)
                            self.clean_fields(False, t)
                    else:
                        for t in self.t_size[self.rank]:
                            print('   Reconstructing fields for time = '+ str(t)+ ' at time = ' + date(timex.gmtime())+ ' with rank = ' + str(self.rank))
                            self.write_areaFields(path,t)
                            self.write_edgeFields(path,t)
                            self.clean_fields(False, t)
                else:
                    print('   Reconstructing fields for time = '+ str(time)+ ' at time = ' + date(timex.gmtime()))
                    self.write_areaFields(path,time)
                    self.write_edgeFields(path,time)
                    self.clean_fields(False, time)
                
        def write_areaFields(self, path, t):
                self.write_h(path,t)
                self.write_Cv(path,t)
                self.write_deltac0(path,t)
                self.write_deltah0(path,t)
                self.write_Us(path,t)
                self.write_tau(path,t)
                self.write_pb(path,t)

        def write_edgeFields(self, path, t):
                self.write_Q(path,t)
                self.write_phi2s(path,t)

        def write_h(self, path, time):
                if (self.h_flag == 'on' or self.h_flag == 'yes' or self.h_flag == True):
                    self.create_h(False, time)
                    self.write_scalarfield(path, 'h', self.h, time, 'areaScalarField')

        def write_Cv(self, path, time):
                if (self.Cv_flag == 'on' or self.Cv_flag == 'yes' or self.Cv_flag == True):
                    self.create_Cv(False, time)
                    self.write_scalarfield(path, 'Cv', self.Cv, time, 'areaScalarField')

        def write_deltah0(self, path, time):
                if (self.deltah0_flag == 'on' or self.deltah0_flag == 'yes' or self.deltah0_flag == True):
                    self.create_deltaz0(False, time)
                    self.write_scalarfield(path, 'deltah0', self.deltah0, time, 'areaScalarField') 

        def write_deltac0(self, path, time):
                if (self.deltac0_flag == 'on' or self.deltac0_flag == 'yes' or self.deltac0_flag == True):
                    self.create_deltaz0(False, time)
                    self.write_scalarfield(path, 'deltac0', self.deltac0, time, 'areaScalarField') 

        def write_Us(self, path, time):
                if (self.Us_flag == 'on' or self.Us_flag == 'yes' or self.Us_flag == True):
                    self.create_Us(False, time)
                    self.write_scalarfield(path, 'Us', self.Us, time, 'areaVectorField')

        def write_tau(self, path, time):
                if (self.tau_flag == 'on' or self.tau_flag == 'yes' or self.tau_flag == True):
                    self.create_tau(False, time)
                    self.write_scalarfield(path, 'tau', self.tau, time, 'areaVectorField')

        def write_pb(self, path, time):
                if (self.pb_flag == 'on' or self.pb_flag == 'yes' or self.pb_flag == True):
                    self.create_pb(False, time)
                    self.write_scalarfield(path, 'pb', self.pb, time, 'areaScalarField')
                    
        def write_Q(self, path, time):
                if (self.Q_flag == 'on' or self.Q_flag == 'yes' or self.Q_flag == True):
                    self.create_Q(False, time)
                    self.write_scalarfield(path, 'Q', self.Q, time, 'edgeScalarField')

        def write_phi2s(self, path, time):
                if (self.phi2s_flag == 'on' or self.phi2s_flag == 'yes' or self.phi2s_flag == True):
                    self.create_phi2s(False, time)
                    self.write_scalarfield(path, 'phi2s', self.phi2s, time, 'edgeScalarField')

        def write_scalarfield(self, path, name, field, time, header_class):
                if time == -1:
                    for t in range(len(self.t)):
                        fi = open(path+'/'+str(self.t[t])+'/'+name, 'w')
                        fi.write(header(name, header_class, time))
                        self.write_dims(fi, field[str(self.t[t])].dim) 
                        self.write_field(fi, field[str(self.t[t])].iF.field, field[str(self.t[t])].ty)
                        self.write_bfield(fi, field[str(self.t[t])].bF.d, field[str(self.t[t])].ty)                        
                        fi.write(footer())                           
                        fi.close()
                else:
                    fi = open(path+'/'+str(time)+'/'+name, 'w')
                    fi.write(header(name, header_class, time))
                    self.write_dims(fi, field[str(time)].dim) 
                    self.write_field(fi, field[str(time)].iF.field, field[str(time)].ty)
                    self.write_bfield(fi, field[str(time)].bF.d, field[str(time)].ty)
                    fi.write(footer())                       
                    fi.close()                    
            
        def write_field(self, fi, field, type):
                fi.write("\n")
                fi.write('internalField   nonuniform List<'+type+'>'+"\n")
                fi.write(str(len(field))+"\n")
                fi.write('('+"\n")
                if type == 'scalar':
                    for i in range(len(field)):
                        fi.write(str(field[i])+"\n") 
                elif type == 'vector':
                    for i in range(len(field)):
                        fi.write('('+str(field[i].x)+' '+str(field[i].y)+' '+str(field[i].z)+')'+"\n") 
                fi.write(')'+"\n")
                fi.write(';'+"\n")
                fi.write("\n") 

        def write_dims(self, fi, dim):
                fi.write("\n")
                fi.write('dimensions'+'      [')
                List = dim.List()
                for i in range(len(List)-1):
                    fi.write(str(int(List[i]))+' ')
                fi.write(str(int(List[-1])))
                fi.write('];'+"\n")
                fi.write("\n") 

        def write_bfield(self, fi, bfield, type):
                fi.write('boundaryField'+"\n")
                fi.write('{'+"\n")
                
                for patch_key in bfield.keys():
                    fi.write('    '+str(patch_key)+"\n")
                    fi.write('    '+'{'+"\n")
                    if bfield[patch_key]['type']  == 'zeroGradient':
                        fi.write('        type            '+'zeroGradient'+';'+"\n")                             
                    elif bfield[patch_key]['type']  == 'fixedValue':
                        fi.write('        type            '+'fixedValue'+';'+"\n")
                        fi.write('        value           nonuniform List<'+type+'> '+"\n")
                        fi.write(str(len(bfield[patch_key]['value']['List']))+"\n")
                        fi.write('('+"\n")
                        if type == 'scalar':
                            for i in range(len(bfield[patch_key]['value']['List'])):
                                fi.write(str(bfield[patch_key]['value']['List'][i])+"\n")                             
                        elif type == 'vector':
                            for i in range(len(bfield[patch_key]['value']['List'])):
                                v = bfield[patch_key]['value']['List'][i]
                                fi.write('('+str(v.x)+' '+str(v.y)+' '+str(v.z)+')'+"\n") 
                        fi.write(')'+"\n")
                        fi.write(';'+"\n")
                    elif bfield[patch_key]['type']  == 'calculated':
                        fi.write('        type            '+'calculated'+';'+"\n") 
                        fi.write('        value           nonuniform List<'+type+'> '+"\n")
                        fi.write(str(len(bfield[patch_key]['value']['List']))+"\n")
                        fi.write('('+"\n")
                        if type == 'scalar':
                            for i in range(len(bfield[patch_key]['value']['List'])):
                                fi.write(str(bfield[patch_key]['value']['List'][i])+"\n")                             
                        elif type == 'vector':
                            for i in range(len(bfield[patch_key]['value']['List'])):
                                v = bfield[patch_key]['value']['List'][i]
                                fi.write('('+str(v.x)+' '+str(v.y)+' '+str(v.z)+')'+"\n") 
                        fi.write(')'+"\n")
                        fi.write(';'+"\n")                          
                    fi.write('    '+'}'+"\n")                               
                fi.write('}'+"\n")
                fi.write("\n") 

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
                self.h_proc.pop(str(t), None)
                self.Cv_proc.pop(str(t), None)
                self.deltac0_proc.pop(str(t), None)
                self.deltah0_proc.pop(str(t), None)
                self.Us_proc.pop(str(t), None)
                self.pb_proc.pop(str(t), None)
                self.tau_proc.pop(str(t), None)                
                self.h.pop(str(t), None)
                self.Cv.pop(str(t), None)
                self.deltac0.pop(str(t), None)
                self.deltah0.pop(str(t), None)
                self.Us.pop(str(t), None)
                self.pb.pop(str(t), None)
                self.tau.pop(str(t), None) 
                
        def clean_edgeFields(self, t):
                self.Q_proc.pop(str(t), None)
                self.phi2s_proc.pop(str(t), None)
                self.Q.pop(str(t), None)
                self.phi2s.pop(str(t), None)


 ##########################################  Main function  ##############################################################

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
                          
def main(argv):
    from mpi4py import MPI
    import time
    
    parser = argparse.ArgumentParser(description='Reconstructing the outputs of a debrisfaSavageHutterFoam simulation')

    parser.add_argument('-path', type=str, help='input simulationPath',default='')
    parser.add_argument('-o', type=str, help='output folder',default='')
    parser.add_argument('-h_flag', help='Activate to reconstruct h field', action="store_true")
    parser.add_argument('-Cv_flag', help='Activate to reconstruct Cv field ', action="store_true")
    parser.add_argument('-deltaz0_flag', help='Activate to reconstruct deltaz0 field ', action="store_true")
    parser.add_argument('-Us_flag', help='Activate to reconstruct Us field ', action="store_true")
    parser.add_argument('-Q_flag', help='Activate to reconstruct Q field ', action="store_true")
    parser.add_argument('-tau_flag', help='Activate to reconstruct tau field ', action="store_true")
    parser.add_argument('-pb_flag', help='Activate to reconstruct pb field ', action="store_true")
    parser.add_argument('-phi2s_flag', help='Activate to reconstruct phi2s field ', action="store_true")
    parser.add_argument('-withZero', help='Activate to reconstruct fields for time 0', action="store_true")
        
    args = parser.parse_args()

    path = args.path
    o = args.o
    h = args.h_flag
    Cv = args.Cv_flag
    pb = args.pb_flag
    deltaz0 = args.deltaz0_flag
    Us  = args.Us_flag
    tau = args.tau_flag
    Q = args.Q_flag   
    phi2s = args.phi2s_flag

    withZero = args.withZero
            
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print('Rank = ' + str(rank)+'  started at = ' + str(date(time.gmtime())))

    if path == '':
        path = os.getcwd()

    if o != '':        
        clean_dir(path, o, rank)
        make_dir(path, o, rank)
        output = path+'/'+o
    else:
        output = path
        
    if rank == 0:
        r = runCase(path, h, Cv, deltaz0, Us, Q, pb, tau, phi2s, rank)
    else:
        r = None
        
    r = comm.bcast(r, root=0)
    r.rank = rank
    r.write_output(output, size, -1, withZero)
    
    print('Rank = ' + str(r.rank)+'  ended at = ' + str(date(time.gmtime())))

if __name__ == "__main__":
    main(sys.argv[1:]) 