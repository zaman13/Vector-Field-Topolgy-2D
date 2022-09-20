#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:14:38 2022

@author: Mohammad Asif Zaman
"""

import numpy as np
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology


from critical_points_lines import critical_point_type

def detect_local_minima(arr):
    # Code taken from the following sources:
    # https://stackoverflow.com/questions/3986345/how-to-find-the-local-minima-of-a-smooth-multidimensional-array-in-numpy-efficie
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    
    eroded_background = morphology.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min ^ eroded_background
    
    
    
    return np.where(detected_minima) 




# ====================================================================================================
# Correct for corner points
# ====================================================================================================
def adj_local_minima(x,y,F):
    
    ind_lm = detect_local_minima(F)
    #  ind_lm[1][..] is the x indices and ind_lm[0][...] are the y indices
    
    # Adjust to see if the indices do not point to corner points
    # =========================================================================
    xL = min(x)
    xH = max(x)
    yL = min(y)
    yH = max(y)
    
    adj_ind = np.where( (x[ind_lm[1]] != xL) & (x[ind_lm[1]] != xH) & (y[ind_lm[0]] != yH) & (y[ind_lm[0]] != yL) )
    
    # =========================================================================
    
    #ind_out is an array of arrrays. ind_out[1][..] is the x indices and ind_out[0][...] are the y indices
    
    ind_out = ind_lm[0][adj_ind],ind_lm[1][adj_ind]   
    
    return ind_out


# ====================================================================================================



# ====================================================================================================
# Correct for multiple critical points close to each other
# ====================================================================================================
def thin_local_minima(x,y,ind,th):
    
    counter = 0
    ox = []
    oy = []
    
    # Loop over all the points except the last one
    for m in range(len(ind[0])-1):
        flag = 0
        x1 = x[ind[1][m]]
        y1 = y[ind[0][m]]
        
        # For each point, loop over the other right hand points in the list 
        for n in range(m+1,len(ind[0])):
                
            x2 = x[ind[1][n]]
            y2 = y[ind[0][n]]
            
            r2 = (x2-x1)**2 + (y2-y1)**2
            
            if r2 < th**2:   # check if any of the n loop points are close to the m loop point
                flag = 1     # if so, we will ignore the m point
                break        # break with the flag value 1
            
        if flag == 0:        # if the m point was found to be sufficiently far away from all the n points (i.e. flag == 0), we store it
            ox.append( ind[1][m] )
            oy.append( ind[0][m] )
            
    
    # Note that the loop starts eliminating points from the left of the list. The right most point (end of the list) is thus
    # always kept
    
    # Add the right most  point to the list
    ox.append( ind[1][-1] )
    oy.append( ind[0][-1] )
            
    
    return np.array(oy), np.array(ox)
# ====================================================================================================


# ====================================================================================================
# Find the average of the function value near a potential critical point
# ====================================================================================================
def local_average(F,ix,iy,Nx,Ny):
    
    
    
    # Find the local average of F around F[iy,ix] taking into account Nx and Ny neigboring points
    
    # Note, np.size(F,1) = len(x) and np.size(F,0) = len(y)
    
    
    # The region over which the averaging will be performed. Checking for boundary crossings and setting the 
    # limits appropriately
    ix_low = ix - Nx if (ix-Nx) >=0 else 0
    iy_low = iy - Ny if (iy-Ny) >=0 else 0
    
    ix_high = ix + Nx if (ix+Nx) < np.size(F,1) else np.size(F,1)-1
    iy_high = iy + Ny if (iy+Ny) < np.size(F,0) else np.size(F,0)-1
    
    # F_mean_x = np.zeros(iy_high-iy_low)
    counter = 0
    F_mean = 0
    for m in range(iy_low,iy_high+1):
        F_mean =  F_mean + np.mean(F[m][ix_low:ix_high])
        counter = counter + 1
    
    
    # return ix_low, ix_high, iy_low, iy_high,counter, F_mean
    return abs(F_mean/counter)
# ====================================================================================================
    


# ====================================================================================================
# Calculate local average for all the potential critical points    
# ====================================================================================================    
def avgF_at_local_minima(F,ind,Nx,Ny):
    
    avgF =  np.zeros(len(ind[0]))
    
    for m in range(len(ind[0])):
        avgF[m] = local_average( F, ind[1][m], ind[0][m] ,Nx,Ny)
    
        
    return avgF
# ====================================================================================================    
    



# ====================================================================================================   
# Estimate the number of critical points from local averages
# ====================================================================================================    
def Npoint_estimate(avgF_sort,th,order):
    # Estimate the number of critical points
    # For order == 0
    # Select the points with high local average values, meaning it was surrounded by zeros.     
    # As avF_sort is sorted, the [0] element is the maximum value. For points with local average th times
    # lower than avF, we disregard them.
    # For order == 1, the process is reversed (high local averages are discarded). The [0] element is the minimum value
    Npoints = 0
    for m in range(len(avgF_sort)):
        temp = avgF_sort[m]/avgF_sort[0]
        
        
        flg = (temp < th) if order == 0 else (temp > 1/th)
        
        print('m = %i, ratio = %1.3f, flag = %i \n' % (m,temp,flg))
        
        if flg == 1:
            break
        
        Npoints = Npoints + 1
  
    
    return Npoints
        
# ====================================================================================================       


# ====================================================================================================
# Sort the critical points according to the indices ind_sort
# ====================================================================================================    
def adj_with_rank(ind_adj, ind_sort, Npoints):
    ind_adj2_x = np.zeros(Npoints, dtype = int)
    ind_adj2_y = np.zeros(Npoints, dtype = int)

    for m in range(Npoints):
        ind_adj2_x[m] = int( ind_adj[1][ind_sort[m]] )
        ind_adj2_y[m] = int( ind_adj[0][ind_sort[m]] )
        
    ind_adj2 = ind_adj2_y, ind_adj2_x
    
    return ind_adj2
# ====================================================================================================


# ====================================================================================================
# Print data for the critical points
# ====================================================================================================
def print_list(x,y,Fx,Fy,ind_adj,Navg_x,Navg_y):
    
    F =  np.sqrt(Fx**2 + Fy**2)                             # Vector field magnitude 
    avF = avgF_at_local_minima(F,ind_adj,Navg_x,Navg_y)
    
    c_type = critical_point_type(x, y, Fx, Fy, ind_adj)
    text_map = ['Repl N', 'Attr N', 'Saddle', 'Repl F','Attr F', 'Center']
    
    print('Number of potential CP = %2i \n' % len(ind_adj[0]))
    
    for n in range(len(ind_adj[0])):
        
        xc = np.array([x[ind_adj[1][n]]])
        yc = np.array([y[ind_adj[0][n]]])
        
        txt = text_map[c_type[n]]
        
        print('CP %i: (%+1.2f,%+1.2f), type = %s, LA=%1.2e '  % (n,xc,yc,txt,avF[n]) )
       

    print('\n') 
    
    return 0
# ====================================================================================================
    
    
    
    
    
    
    
    
    