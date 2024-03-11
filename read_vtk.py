import os
import csv
import pandas as pd
import numpy as np

input_Path  = '.....'  #insert the location of the *.vtk files


folders = os.listdir(input_Path)

nSamples=len(folders)

def read_vtkfiles(path):
    f=open(path, "r")
    a=""
    temp = False
    for line in f:
        if "POINTS" in line:
            temp = True
            continue
        elif "POLYGONS" in line:
            temp = False
            continue
        elif temp:
            a=a+line
    f.close()
    
    b=a.strip()
    c=b.replace("\n","")
    d=c.split(" ")
    x = np.array(d)
    points = x.astype(np.float)

    return(points)
    
    
    
all_meshes_points=[]
all_file_names=[]
srepNumber=1
for i in folders:
    tempFilesList = os.listdir(input_Path+"/"+i)
    
    for j in tempFilesList:
        if ".vtk" in j:
            print(j)
            upSpokeFileLocation=input_Path+"/"+i+"/"+j
            points=read_vtkfiles(path=upSpokeFileLocation)
            
            all_meshes_points.append(points)
            all_file_names.append(j)

    srepNumber=srepNumber+1

data={}
df = pd.DataFrame(data) 

for i in range(len(all_file_names)):
    df[all_file_names[i]]=all_meshes_points[i]


df.to_csv("output.csv", sep=',',index=False)  # output.csv is the name of the output file