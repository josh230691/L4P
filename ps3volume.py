#! /usr/bin/env python

#python stuff that may or may not be useful
from __future__ import with_statement
import math, cmath, string, os, sys, fileinput, pprint
from optparse import OptionParser
import random as rnd

#comment out the following if not using matplotlib and numpy
import mpmath as mp
import numpy as np
import pylab as pl
from scipy import interpolate, signal
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import vegas

mhiggs = 125
mb = 4.18

# As in the PDG notes write up
def m232plus(m122):

    # Define constants in the expression
    a = mhiggs**4 / (4 * m122)
    b = m122**2 / mhiggs**4
    c = mhiggs**2 / m122
    d = (4 * mb**2) / m122

    return a * (1 - b * ((-1)*(c - 1) + math.sqrt(1 - d))**2)

def m232minus(m122):

    # Define constants in the expression
    a = mhiggs**4 / (4 * m122)
    b = m122**2 / mhiggs**4
    c = mhiggs**2 / m122
    d = (4 * mb**2) / m122

    return a * (1 - b * ((c - 1) + math.sqrt(1 - d))**2)

def deltam232(m122):

    return m232plus(m122) - m232minus(m122)

# The integrand restricted to the region defined by m232plus and m232minus
def f(x):

    m122 = x[0]
    m232 = x[1]

    # If out of region return 0
    if m232 > m232plus(m122) or m232 < m232minus(m122):
        return 0
    # Else return the integrand
    else:
        return (1 / 256) * (1 / math.pi**3) * (1 / mhiggs**3)

m122max = mhiggs**2
m122min = 4 * mb**2
deltam122 = m122max - m122min

# Number of evaluations for Monte Carlo
N = 100000
sum = 0

# Monte Carlo integration of the integrand

for i in range(0, N):

    # Random m122
    m122_i = m122min + (rnd.random() * deltam122)
    # Random m232 in range allowed by given m122
    m232_i = m232minus(m122_i) + (rnd.random() * deltam232(m122_i))

    sum += deltam122 * deltam232(m122_i) * f([m122_i, m232_i])

result = sum / N
print(result)

# VEGAS implementation of the integral

# Define region to include all of the region defined by m232plus and m232minus
integ = vegas.Integrator([[m122min, m122max], [mb**2, (mhiggs - mb)**2]])
result = integ(f, nitn=10, neval=100000)
print(result.summary())
print('total CS(pb) = %s    Q = %.2f' % (result, result.Q))
                                  
# Code to generate Dalitz plot

# Number of points for each of m232plus and minus
n = 200
# Store generated values of m232minus
#data1 = [[]] * n
# Store generated values of m232plus
#data2 = [[]] * n

#for i in range(0, n):

    # Generate n values for m122_i
    #m122_i = m122min + (i * deltam122) / n
    # Add new data point for m232minus
    #data1[i].append([m122_i, m232minus(m122_i)])
    # Add new data point for m232plus
    #data2[i].append([m122_i, m232plus(m122_i)])

# Convert to numpy array
#data1 = np.array(data1)
#data2 = np.array(data2)

#a, b = data1.T
#c, d = data2.T

#plt.scatter(a, b)
#plt.scatter(c, d)
#plt.show()




    
