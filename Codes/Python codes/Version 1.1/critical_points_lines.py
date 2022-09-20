#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  7 12:22:42 2022

@author: Mohammad Asif Zaman
"""


from __future__ import print_function   


import numpy as np
import pylab as py


from diff_matrices import Diff_mat_1D, Diff_mat_2D   




def Jacobian_mat(x,y,Fx,Fy):
    
    Nx = len(x)
    Ny = len(y)

        
    # Loading finite difference matrix operators
    #================================================================================================
    Dx_2d, Dy_2d, D2x_2d, D2y_2d = Diff_mat_2D(Nx,Ny)   # Calling 2D matrix operators from funciton

    
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    
    Dx_2d = Dx_2d/(2*dx)
    Dy_2d = Dy_2d/(2*dy)
    
    
    D2x_2d = D2x_2d/(dx**2)
    D2y_2d = D2y_2d/(dy**2)
    
    #================================================================================================
    
    
    Jxx = (Dx_2d*Fx.ravel()).reshape(Ny,Nx)
    Jxy = (Dy_2d*Fx.ravel()).reshape(Ny,Nx)
    Jyx = (Dx_2d*Fy.ravel()).reshape(Ny,Nx)
    Jyy = (Dy_2d*Fy.ravel()).reshape(Ny,Nx)

    return Jxx, Jxy, Jyx, Jyy    


def critical_point_type(x,y,Fx,Fy,ind_adj):
    
    th = 1e-5   # Threshold value for considering the real part of the Eigen value to be zero
    Jxx, Jxy, Jyx, Jyy = Jacobian_mat(x,y,Fx,Fy)
    
    cr_type = np.zeros(len(ind_adj[0]),dtype=np.uint32) + 23
    
    for ki in range(len(ind_adj[0])):
            
        
        xc = np.array([x[ind_adj[1][ki]]])
        yc = np.array([y[ind_adj[0][ki]]])
        
        
                
        
        Jo = np.zeros([2,2])
        
        Jo[0,0] = Jxx[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[0,1] = Jxy[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[1,0] = Jyx[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[1,1] = Jyy[ind_adj[0][ki], ind_adj[1][ki]]
        
        eig_val,eig_vect = np.linalg.eig(Jo)
        
        eig_i = np.imag(eig_val)
        eig_r = np.real(eig_val)
        
        if np.sqrt(eig_i[0]**2 + eig_i[1]**2) < th:      # if imaginary part of both eigen values are zero
            if (eig_r[0] > 0) & (eig_r[1] > 0):
                cr_type[ki] = 0                     # repelling node
            if (eig_r[0] < 0) & (eig_r[1] < 0):
                cr_type[ki] = 1                     # attracting node
            if eig_r[0]*eig_r[1] < 0:
                cr_type[ki] = 2                     # saddle node
        
        else:
            
            if (eig_r[0] > 0) & (eig_r[1] > 0):
                cr_type[ki] = 3                     # repelling focus
            if (eig_r[0] < 0) & (eig_r[1] < 0):
                cr_type[ki] = 4                     # attracting focus
            if (eig_r[0] == 0) & (eig_r[1] == 0):
                cr_type[ki] = 5                     # center
            
            # adjustments for non-ideal cases
            if (abs(eig_r[0]) < th) & (abs(eig_r[1]) < th):
                cr_type[ki] = 5                     # center
            
                
                
    return cr_type
        
        
        
def stream_lines(x,y,Fx,Fy,ind_adj):
    
    F =  np.sqrt(Fx**2 + Fy**2)
    text_map = ['Repelling Node', 'Attracting Node', 'Saddle Node', 'Repelling Focus','Attracting Focus', 'Center']
    color_map = ['#5c9090', '#ff756d', '#0f4c81',  '#e9b666','#f5b19c', '#85de77']
    
    
    Jxx, Jxy, Jyx, Jyy = Jacobian_mat(x,y,Fx,Fy)
    
    c_type = critical_point_type(x, y, Fx, Fy, ind_adj)
    
    
    for ki in range(len(ind_adj[0])):
        
            
        
        xc = np.array([x[ind_adj[1][ki]]])
        yc = np.array([y[ind_adj[0][ki]]])
        
        # print('Critical Point no. %i at (%1.2f,%1.2f) Type = %s' % (ki,xc,yc,text_map[c_type[ki]]) )
        # 
        py.plot(xc,yc,'o',color = color_map[c_type[ki]])
        
        
        Jo = np.zeros([2,2])
        
        Jo[0,0] = Jxx[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[0,1] = Jxy[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[1,0] = Jyx[ind_adj[0][ki], ind_adj[1][ki]]
        Jo[1,1] = Jyy[ind_adj[0][ki], ind_adj[1][ki]]
        
        w,v = np.linalg.eig(Jo)
        
        v = np.real(v)

        # seed_points = np.array([[-2, -1, 0, 1, 2, -1], [-2, -1,  0, 1, 2, 2]])
        
        del_t = .1 if c_type[ki] != 5 else .5
        
        seed_points = np.array([ xc + v[0][0]*del_t , yc + v[1][0]*del_t])
        py.streamplot(x, y, Fx, Fy, color=F, linewidth=1, cmap='inferno', start_points=seed_points.T)
        # py.streamplot(x, y, Fx, Fy,  start_points=seed_points.T)
        
        seed_points = np.array([ xc - v[0][0]*del_t , yc - v[1][0]*del_t])
        py.streamplot(x, y, Fx, Fy, color=F, linewidth=1, cmap='inferno', start_points=seed_points.T)
        
        
        seed_points = np.array([ xc + v[0][1]*del_t , yc + v[1][1]*del_t])
        py.streamplot(x, y, Fx, Fy, color=F, linewidth=1, cmap='inferno', start_points=seed_points.T)
        
        seed_points = np.array([ xc - v[0][1]*del_t , yc - v[1][1]*del_t])
        py.streamplot(x, y, Fx, Fy, color=F, linewidth=1, cmap='inferno', start_points=seed_points.T)
        
        
    return 0




        
        
 
    
 
    
 
    
 
    
 
    
 
    
 
    
 
    