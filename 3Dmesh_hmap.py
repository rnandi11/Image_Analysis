#!/usr/bin/env python

import binarymeshformat as bmf
import vedo.mesh
import vedo
from vedo import *
import pandas as pd
import numpy as np
from skimage import io
import tifffile
import matplotlib.path as mplPath
from matplotlib.path import Path
import math
from matplotlib.image import imread
from matplotlib.image import imsave
from scipy.signal import convolve2d as conv2
from matplotlib.patches import Polygon
import scipy as sp
from scipy import ndimage
import cv2
import sys
import time

## Convert bmf scales to movie scales
def convert(x,y,z):
    px=[1944,1024,72]
    dv=[0.2075665,0.2075665,0.4794330]
    width = px[0]*dv[0]
    height = px[1]*dv[1]
    depth = px[2]*dv[2]
    scale = max(width, height, depth)
    offset = [0.5*width / scale, 0.5*height/scale, 0.5*depth/scale ]
    return (
        (x+ offset[0])*scale, 
        (y+ offset[1])*scale, 
        (z+ offset[2])*scale
        )

###########################################################################

## get list of points and triangles
def getVedoMesh( bmfmesh ):
    points = [];
    triangles = []
    p = bmfmesh.positions
    t = bmfmesh.triangles
    
    for i in range( len(bmfmesh.positions)//3):
        x=p[3*i]
        y=p[3*i+1]
        z=p[3*i+2]
        pt=convert(x,y,z)
        #pt=(x,y,z)
        point = [ int(pt[0]/0.2075665), int(pt[1]/0.2075665), int(pt[2]/0.4794330)]
        points.append(point)
    for i in range( len(bmfmesh.triangles) // 3 ):
        triangle = [ t[i*3] , t[i*3+1], t[i*3 + 2] ]
        triangles.append( triangle );
    
    return vedo.mesh.Mesh( [ points, triangles ] )

###########################################################################

# calculate centroid of vmesh object
def calculate_mesh_centroid(vmesh):
    total_points = len(vmesh.points())
    centroid_sum = [0.0, 0.0, 0.0]

    for point in vmesh.points():
        centroid_sum[0] += point[0]
        centroid_sum[1] += point[1]
        centroid_sum[2] += point[2]

    centroid = [
        centroid_sum[0] / total_points,
        centroid_sum[1] / total_points,
        centroid_sum[2] / total_points
    ]
    
    return centroid    

###########################################################################

#Calculates volume inside each mesh
def FillMeshes(bmf_meshes, frame, name):
    image = np.zeros((72, 1024, 1944), dtype=np.uint8)
    dv=[0.2075665,0.2075665,0.4794330]

    for bmesh in bmf_meshes:
        
        volcount = 0
        vmesh = getVedoMesh(bmesh)
        CENTROID = calculate_mesh_centroid(vmesh)
        boundar = np.array(vmesh.bounds(), dtype=int)

        x = np.arange(boundar[0], boundar[1] + 1)
        y = np.arange(boundar[2], boundar[3] + 1)
        z = np.arange(boundar[4], boundar[5] + 1)
        points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

        # Get the points inside the `vmesh`
        inside_points = vmesh.inside_points(points)

        # Extract points and IDs
        if isinstance(inside_points, vedo.pointcloud.Points):
            inside_points_array = inside_points.points()  # Extract coordinates of inside points
            is_inside = inside_points.pointdata["IsInside"]  # Extract the 'IsInside' field if available

            '''
            # Debugging output
            print(f"Type of inside_points_array: {type(inside_points_array)}")
            print(f"Shape of inside_points_array: {inside_points_array.shape}")
            print(f"Content of inside_points_array: {inside_points_array}")
            '''
            if is_inside is not None:
                print(f"IsInside field: {is_inside}")

            for point in inside_points_array:
                x, y, z = point.astype(int)
                image[z, y, x] = 255
                volcount += 1
            exactvol.append([name,frame,CENTROID[0] * dv[0],CENTROID[1] * dv[1],CENTROID[2] * dv[2],
                volcount * (dv[0] * dv[1] * dv[2])])

        else:
            print(f"Unexpected type for inside_points: {type(inside_points)}")
            continue

    #create image file for each frame
    #tifffile.imwrite(f'image_folder/image_GFO_{frame}.tif', image)

    # Create the sum projection along the z-axis
    sum_image = np.sum(image, axis=0).astype(np.float32)
    tifffile.imwrite(f'image_folder/sum_image_GFO_{frame}.tif', sum_image)

    return exactvol

###########################################################################


##Calculate mesh volume from the 2D height map
def Heightmap(bmf_meshes, frame,name):

    dv=[0.2075665,0.2075665,0.4794330]
    summed_zstack=Image.open('image_folder/SUM_image_GFO_'+str(frame).zfill(1)+'.tif')
    Zstack_array=np.array(summed_zstack)
    '''
    height, width = Zstack_array.shape
    output=np.zeros((height,width,3), dtype=float)
    #print(height,width)
    '''
    zplane=[]
    for bmesh in bmf_meshes:
        vmesh = getVedoMesh(bmesh)
        boundar = np.array(vmesh.bounds(), dtype=int)
        zplane.append([boundar[4]])
    min_z=int(np.min(np.array(zplane)))
    #print(min_z)   
    for bmesh in bmf_meshes:
        volcount=0
        vol=0
        vmesh = getVedoMesh(bmesh)
        CENTROID=calculate_mesh_centroid(vmesh)
        boundar = np.array(vmesh.bounds(), dtype=int)
        x = np.arange(boundar[0], boundar[1] + 1)
        y = np.arange(boundar[2], boundar[3] + 1)

        #This criteria is somewhat arbitrary: chosen all points in 3 z-planes such that 
        #there are enough points to mimic the continuous shape of the apical layer.
    
    
        z = np.arange(min_z+12,min_z+15)                            
        points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)
        point_list=[]
        # Get the points inside the `vmesh`
        inside_points = vmesh.inside_points(points)

        # Extract points and IDs
        if isinstance(inside_points, vedo.pointcloud.Points):
            inside_points_array = inside_points.points()  # Extract coordinates of inside points
            is_inside = inside_points.pointdata["IsInside"]  # Extract the 'IsInside' field if available
        for point in inside_points_array:
            x, y, z = point.astype(int)
            if([x,y] not in point_list):
                point_list.append([x,y])
                vol += (Zstack_array[y, x])
                volcount+=1
                    #output[y,x,:]=[255,0,0]
        if (volcount!=0):    
            heightmap.append([name,frame,CENTROID[0]*dv[0],CENTROID[1]*dv[1],(vol*dv[0]*dv[1]*dv[2]/255),
                (vol*dv[2]/(255*volcount))])

    return heightmap

if __name__ == "__main__":
    tracks = bmf.loadMeshTracks('/path/to/meshfile.bmf')
    print("loaded %s tracks"%len(tracks))
    '''
    # Print attributes of the first track in the tracks list
    if len(tracks) > 0:
        print("Attributes of a track object:")
        print(dir(tracks[0]))
    '''
    exactvol=[] 
    heightmap=[]
    frames = set()

    for track in tracks:
        for k in track.meshes:
            frames.add(k)
    
    for frame in frames:
        print(frame)
        to_see = [(track.meshes[frame], track.name) for track in tracks if frame in track.meshes.keys()]
        for mesh, name in to_see:
            FillMeshes([mesh], frame, name)
            Heightmap([mesh], frame, name)
        #FillMeshes(to_see,frame)  

    np.savetxt('exactvol_movie_name.csv',np.array(exactvol),delimiter =",", fmt ='% s')
    np.savetxt('heightmap_movie_name.csv',np.array(heightmap),delimiter=',',fmt='%s')




