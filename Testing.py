#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 20:37:56 2024

@author: joshkevil
"""

import json
from sympy import Matrix, symbols, solve
import numpy as np
from Project_code import LUFactorization
from Project_code import forewardSubstitution
from Project_code import backwardSubstitution
#from Project_code import testBack

A = np.matrix([[1,-2,1],[-2,1,-3],[4,-3,1]])

Sol = LUFactorization(A,3)
B = Sol[0]
D = Sol[1]
C = np.linalg.inv(D)
print(B)
print(C)
print(C*B)
b = [1,2,1]

y = forewardSubstitution(C,b)
print(y)
#print(y[2])
print(backwardSubstitution(B,y))
