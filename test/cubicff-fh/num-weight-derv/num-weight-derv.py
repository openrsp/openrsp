#!/usr/bin/env python

from pylab import *
import sys
import string
from optparse import OptionParser

usage = '''
  utility to numerically approximate numerical grid weight derivatives
  using polynomial fits
  requires python-pylab

  example: ./%prog -p 2 -d 1 --grids="grid_+0.001 grid_+0.000 grid_-0.001" > grid_1'''

parser = OptionParser(usage)

parser.add_option('-p',
                  dest = 'polynomial_order',
                  help = 'order of the polynomial (2 for quadratic)',
                  metavar='N')
parser.add_option('-d',
                  dest = 'derivative',
                  help = 'derivative (2 for second derivative)',
                  metavar='N')
parser.add_option('--grids',
                  dest = 'grids',
                  help = 'grid files to polyfit (space separated in "")',
                  metavar='"grid_-0.001 grid_+0.000 grid_+0.001 ..."')

(options, args) = parser.parse_args()

if len(sys.argv) == 1:
    # user has given no arguments: print help and exit
    print parser.format_help().strip()
    sys.exit()

#-------------------------------------------------------------------------------

def fact(x):
    return (1 if x==0 else x * fact(x-1))

def get_derivative(x_l, y_l, p_order, d_order):
    c_l = polyfit(x_l, y_l, p_order)
    d_l = []
    for i in range(d_order + 1):
        d_l.append(fact(i)*c_l[len(c_l) - i - 1])
    return d_l[d_order]

def readlines_strip(file_name):
    f = open(file_name, 'r')
    l = [line.strip('\n') for line in f.readlines()]
    f.close()
    return l

#-------------------------------------------------------------------------------

class grid:
    def __init__(self, file_name):
        self.x = []
        self.y = []
        self.z = []
        self.w = []
        for line in readlines_strip(file_name):
            if len(line.split()) > 1:
                self.x.append(float(string.replace(line.split()[0], 'D', 'e')))
                self.y.append(float(string.replace(line.split()[1], 'D', 'e')))
                self.z.append(float(string.replace(line.split()[2], 'D', 'e')))
                self.w.append(float(string.replace(line.split()[3], 'D', 'e')))
        self.nr_points = len(self.x)
        self.f = float(file_name.split('_')[-1])

#-------------------------------------------------------------------------------

grid_d = {}
f_l    = []
undiff_file_name = None
for file_name in options.grids.split():
    grid_d[file_name] = grid(file_name)
    f = grid_d[file_name].f
    f_l.append(f)
    if abs(f) < 1e-10:
        undiff_file_name = file_name

polynomial_order = int(options.polynomial_order)
derivative       = int(options.derivative)

nr_points = grid_d[undiff_file_name].nr_points

print nr_points

for i in range(nr_points):
    x = grid_d[undiff_file_name].x[i]
    y = grid_d[undiff_file_name].y[i]
    z = grid_d[undiff_file_name].z[i]
    w_l = []
    for file_name in options.grids.split():
        w_l.append(grid_d[file_name].w[i])
    print x, y, z, get_derivative(f_l, w_l, polynomial_order, derivative)

print -1
