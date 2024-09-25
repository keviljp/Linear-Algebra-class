#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 15:46:41 2024

@author: joshkevil
"""

import json
import numpy as np

with open('ProjectOne.txt', 'w')as proj:
    proj.write("")
    
def read_resistances_json(file_name):
    resistances = {}
    with open(file_name, 'r') as f:
        data = json.load(f)
        for entry in data:
            node1 = entry["node1"]
            node2 = entry["node2"]
            resistance = entry["resistance"]
            resistances[(node1, node2)] = resistance
    return resistances

def read_fixed_voltages_json(file_name):
    fixed_voltages = {}
    with open(file_name, 'r') as f:
        data = json.load(f)
        for entry in data:
            node = entry["node"]
            voltage = entry["voltage"]
            fixed_voltages[node] = voltage
    return fixed_voltages

def read_resistances_txt(file_name):
    resistances = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.split(" ")
            node1 = line[0]
            node2 = line[1]
            resistance = line[2]
            resistances[(node1, node2)] = resistance
    return resistances

def read_fixed_voltages_txt(file_name):
    fixed_voltages = {}
    with open(file_name, 'r') as f:
        for line in f:
            line = line.split(" ")
            node = line[0]
            voltage = line[1]
            fixed_voltages[node] = voltage
    return fixed_voltages


fileResistances = input("please enter the name of the resistances you would like to analyze")
fileVoltages = input("please enter the name of the fixed voltages you would like to analyze")
fileType = input("is this a space separated txt file or a json file? (enter txt or json): ")

def setUpNeeded(fileVoltages, fileResistances,fileType):
    if fileType == 'txt':
        fixed_voltages = read_fixed_voltages_txt(fileVoltages)
        resistances = read_resistances_txt(fileResistances)
    elif fileType == 'json':
        fixed_voltages = read_fixed_voltages_json(fileVoltages)
        resistances = read_resistances_json(fileResistances)
    else:
        print("invalid filetype")
        
    return (fixed_voltages, resistances)

Needed = setUpNeeded(fileResistances, fileVoltages, fileType)
fixed_voltages = Needed[0]
resistances = Needed[1]

howMany = max(resistances.keys())
maxMany = max(howMany)
A = np.zeros((maxMany,maxMany))
b = np.zeros(maxMany)

def findNeighbors(i):
    Neighbors = {}
    values = []
    for value in resistances.keys():
        if i in value:
            if value[0]==i:
                values.append(value[1])
            else:
                values.append(value[0])
        Neighbors[i] = values
    return Neighbors
    
def useSetVolts(i):
    A[i-1,i-1] = 1
    b[i-1] = fixed_voltages[i]

def makeA(i):
    sum_resistances = 0
    Neighbors = findNeighbors(i)
    for neighbor in Neighbors[i]:
        if (i, neighbor) in resistances:
            resistance = resistances[(i, neighbor)]
        else:
            resistance = resistances[(neighbor, i)]
        if resistance > 0:
            A[i-1, neighbor-1] = -1 / resistance
            sum_resistances += 1 / resistance
        A[i-1, i-1] = sum_resistances
        
for i in range(1,maxMany+1):
    if i in fixed_voltages.keys():
        useSetVolts(i)
    else:
        findNeighbors(i)
        makeA(i)
A = np.matrix(A)

def LUFactorization (A, maxMany):
    A_kmin1 = A
    L = np.matrix(np.identity(maxMany))
    #for k in range(0,maxMany):
    for k in range(0,maxMany):
        E_k = np.zeros((maxMany,maxMany))
        for l in range (0,maxMany):
            for j in range (0,maxMany):
                if j == l:
                    E_k[l,j] = 1
                elif (j == k) and (j<l):
                    E_k[l,j] = -(A_kmin1[l,k])/(A_kmin1[k,k])
                else:
                    E_k[l,j] = 0
        A_kmin1 = E_k*A_kmin1
        L = E_k*L
    return(A_kmin1, L)

Sol = LUFactorization(A,maxMany)
B = Sol[0]
D = Sol[1]
C = np.linalg.inv(D)

def writeMatrix(title,B,file,n,m):
    with open(file, 'a') as proj:
        proj.write("Matrix "+title+":\n")
        for i in range(n):
            for j in range(m):
                proj.write(str(round(B[i,j],2))+'\t')
            proj.write('\n')
        proj.write('\n')
            
writeMatrix('U',B,'ProjectOne.txt',25,25)
writeMatrix('L',C,'ProjectOne.txt',25,25)



def forewardSubstitution (L,b):
    n = len(b)
    sols = [0]*n
    for i in range(0,n):
        knowns = 0
        for m in range(i):
            knowns += sols[m]*L[i,m]
        unk = b[i]-knowns
        nextsol = unk/L[i,i]
        sols[i] = nextsol
    return (sols)

def backwardSubstitution (U,y):
    n = len(y)
    V = [0]*n
    for i in range(0,n):
        knowns = 0
        for m in range(i):
            knowns += V[n-1-m]*U[n-1-i,n-1-m]
        unk = y[n-1-i]-knowns
        nextsol = unk/U[n-1-i,n-1-i]
        V[n-1-i] = nextsol
    return (V)

y = forewardSubstitution(C,b)
V = backwardSubstitution(B,y)

def compute_currents(V, resistances):
    currents = {}
    for (node1, node2), resistance in resistances.items():
        current = (V[node1 - 1] - V[node2 - 1]) / resistance
        currents[(node1, node2)] = current
    return currents

currents = compute_currents(V, resistances)

print("\nCurrents through each link:")
for (node1, node2), current in currents.items():
    print(f"Current from Node {node1} to Node {node2}: {current:.6f} A")
    with open('ProjectOne.txt', 'a') as proj:
        proj.write(f"Current from Node {node1} to Node {node2}: {current:.6f} A \n")

for i in range(len(V)):
    print('Voltage at node '+str(i+1)+" is "+str(V[i]))
    with open('ProjectOne.txt', 'a') as proj:
        proj.write('Voltage at node '+str(i+1)+" is "+str(V[i])+" V\n")
        






    
