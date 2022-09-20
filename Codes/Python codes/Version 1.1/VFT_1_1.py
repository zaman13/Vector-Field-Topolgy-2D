#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:28:07 2022

@author: Mohammad Asif Zaman

VFT: Finding non-trivial zeros 
"""

from __future__ import print_function   

import os
import numpy as np
import pylab as py





from local_minima import *
from critical_points_lines import *

# from utilities_definition import my_contourf, my_scatter, color_distinct
from diff_matrices import Diff_mat_1D, Diff_mat_2D   



py.close('all')

os.system('clear')


# ====================================================================================================
# Parameters
# ====================================================================================================
# Npoints = 4  # number of critical points to select
Navg_x = 12  # number of points taken when calculating local average (one sided along x direction)
Navg_y = 12  # number of points taken when calculating local average (one sided along y direction)
Npoints_th = .6
th_thin = 0.1  # Thereshold for thinning critical points (i.e. two points with Euclidean dist < th_thin are counted only once)



# ====================================================================================================
# ====================================================================================================



# ====================================================================================================
# Read csv files where multiple datasets are present in the same file as blocks of rows
# ====================================================================================================

def block_csv_read(block_no,fname):
    path = '/home/asif/pCloudDrive/Data/C-aperture enhancement polarization/Current data wavelength sweep/'
    # fname = 'c60_90deg_1064nm_z=-150nm_Jx_fset2.csv'

    f_path =  path + fname

    Nlines = 15464
    Nheader = 13
    block_length = 301
    offset = 2

    # block_no = 7

    Nblock1 = Nheader + block_length

    head_skip = Nheader + block_no*(block_length + offset)
    foot_skip = Nlines - head_skip - block_length    
        

    Mdata = np.genfromtxt(f_path,delimiter=',',skip_header = head_skip, skip_footer = foot_skip)
    return Mdata


# ====================================================================================================
# ====================================================================================================





# ====================================================================================================
# Define various test functions
# ====================================================================================================
def test_vector(ftype):
    
    Nx = 400                         # No. of grid points along x direction
    Ny = 300                         # No. of grid points along y direction
    x = np.linspace(-2,2,Nx)         # x variables in 1D
    y = np.linspace(-2,2,Ny)         # y variable in 1D
    
    
    X,Y = np.meshgrid(x,y)          # 2D meshgrid
    
    
    if ftype == 1:
        Fx = np.exp(-X**2 - Y**2)*(-2*X + Y)
        Fy = np.exp(-X**2 - Y**2)*(3*Y - X)
        
    if ftype == 2:
        Fx = 4*np.exp(-2*X**2 - Y**2)*np.cos(3*X)
        Fy = 2*np.exp(-2*X**2 - Y**2)*np.cos(3.5*Y)
    
    if ftype == 3:
        Fx = np.sin(X+Y) * np.exp(-0.5*X**2 - 0.3*Y**2)
        Fy = np.cos(X-Y) * np.exp(-0.5*X**2 - 0.3*Y**2)
    
    if ftype == 4:
        Fx = (np.sin(X))**2+ Y
        Fy = np.cos(X+Y**2)
    
    if ftype == 5:
        Fx = np.cos(X**2 + Y)
        Fy = X - Y**2 + 1
    
    if ftype == 6:
        Fx = -2*X -3 + Y**2
        Fy = .2*Y+.8*X
        
    if ftype == 7:  # Monkey Saddle
        Fx = 6*X*Y
        Fy = 3*X**2 - 3*Y**2
        
    if ftype == 8:                         
        Fx = Y*np.exp(-X**2 - Y**2)
        Fy = -X*np.exp(-X**2 - Y**2)




    
    # return np.transpose(x), np.transpose(y), np.transpose(Fx), np.transpose(Fy)
    return x,y,Fx,Fy
# ====================================================================================================







# ====================================================================================================
# Setup 
# ====================================================================================================
ftype = 3                               # select function 

x, y, Fx, Fy = test_vector(ftype)       # the spatial variables and the vector function
X,Y = np.meshgrid(x,y)                  # Meshgrid
F =  np.sqrt(Fx**2 + Fy**2)             # Vector field magnitude 

Nx = len(x)
Ny = len(y)

# # Loading finite difference matrix operators
# #================================================================================================
# Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton


# dx = x[1]-x[0]
# dy = y[1]-y[0]

# Dx_2d = Dx_2d/(2*dx)
# Dy_2d = Dy_2d/(2*dy)


# Fx = (Dx_2d*F.ravel()).reshape(Ny,Nx)
# Fy = (Dy_2d*F.ravel()).reshape(Ny,Nx)



py.contourf(x,y,F,cmap = 'inferno')
py.streamplot(x,y,Fx,Fy,color = [.7, .7 , .7],density = 1.2, linewidth = 0.4)

# ====================================================================================================





# ====================================================================================================
# ====================================================================================================
# Find the indices that correspond to the critical points
# ====================================================================================================
# ind_lm = detect_local_minima(F)

ind_1 = adj_local_minima(x,y,F)                   # Find indices of the local minimas. Correct for corner points
c_type = critical_point_type(x, y, Fx, Fy, ind_1)
print('Adjusted Local Minima function call... \n')
print_list(x,y,Fx,Fy,ind_1,Navg_x,Navg_y)

print('Thin Local Minima function call... \n')
ind_2 = thin_local_minima(x,y,ind_1,th_thin)    # Ignore points that are close to each other (select only one)

print_list(x,y,Fx,Fy,ind_2,Navg_x,Navg_y)



# Local average based ranking and selection (not currently used)
# ====================================================================================================
print('Local average selection... \n')
avF = avgF_at_local_minima(F,ind_2,Navg_x,Navg_y)            # Calculate local averages
# ind_sort = np.flip(np.argsort(avF))                         # sort from maximum average value to minimum average value
ind_sort = np.argsort(avF)                                   # sort from minimum average value to maximum average value
Npoints = Npoint_estimate(avF[ind_sort], Npoints_th,1)       # Estimate the number of critical points
ind_3 = adj_with_rank(ind_2, ind_sort, Npoints)              # Clip to the estimated number of critical points

# print_list(x,y,Fx,Fy,ind_3,Navg_x,Navg_y)

# ====================================================================================================
# ====================================================================================================



# print_list(x,y,Fx,Fy,ind_2,Navg_x,Navg_y)



# ====================================================================================================
# Select adjusted critical points
# ====================================================================================================

ind_adj = ind_3                       # Local average raknin algorithm needs more work. It is skipped for now
print('Final selection of CPs \n')
print_list(x,y,Fx,Fy,ind_adj,Navg_x,Navg_y)
# ====================================================================================================



py.figure()

if Npoints > 0:
    stream_lines(x,y,Fx,Fy,ind_adj)
    
py.xlabel('x')
py.ylabel('y')
py.savefig('out1.svg')

py.streamplot(x,y,Fx,Fy,color = [.7, .7 , .7],density = 1.2, linewidth = 0.4)
py.savefig('out2.svg')


py.figure()
py.contourf(x,y,F,cmap = 'Greys')
py.streamplot(x,y,Fx,Fy,color = [.7, .7 , .7],density = 1.2, linewidth = 0.4)

if Npoints > 0:
    stream_lines(x,y,Fx,Fy,ind_adj)
    
py.xlabel('x')
py.ylabel('y')
py.savefig('out3.svg')


