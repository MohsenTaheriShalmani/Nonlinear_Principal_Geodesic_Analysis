import xml.etree.ElementTree as ET
import os
import csv
import pandas as pd
import numpy as np

input_Path  = '............' # insert the location of the *.vtp files

def readXML_data(vtpAddress,spokeType,srepNumber,file):
    
    root = ET.parse(vtpAddress)
    
    Spoke = spokeType
    srepNumber = srepNumber
    file=file
    
    for i in root.iter('DataArray'):
        if(i.attrib["Name"]=="spokeDirection"):
            directions=i.text
        if(i.attrib["Name"]=="spokeLength"):
            lengths=i.text
        if(i.attrib["Name"]=="Points"):
            positions=i.text
        
        
    Directions=(" ".join(directions.split())).split()
    Lengths=(" ".join(lengths.split())).split()
    Positions=(" ".join(positions.split())).split()

    Ux=[]
    Uy=[]
    Uz=[]
    for i in np.arange(0,len(Directions),3):
        Ux.append(float(Directions[i]))
        Uy.append(float(Directions[i+1]))
        Uz.append(float(Directions[i+2]))
    

    Radii=[]
    for i in np.arange(0,len(Lengths)):
        Radii.append(float(Lengths[i]))
    
    xPos=[]
    yPos=[]
    zPos=[]
    for i in np.arange(0,len(Positions),3):
        xPos.append(float(Positions[i]))
        yPos.append(float(Positions[i+1]))
        zPos.append(float(Positions[i+2]))
    
    numberOfSpokes=len(Radii)
    SpokesNumber=np.arange(1,numberOfSpokes+1)
        
    result={"file":file,"srepNumber":srepNumber,"Spoke":Spoke,"Ux":Ux,"Uy":Uy,"Uz":Uz,
            "Radii":Radii,"xPos":xPos,"yPos":yPos,"zPos":zPos, "SpokesNumber":SpokesNumber}
    
    return(result)



folders = os.listdir(input_Path)

nSamples=len(folders)

allSreps=[]
srepNumber=1
for i in folders:
    tempFilesList = os.listdir(input_Path+"/"+i+"/"+i.replace("refined_",""))
    k=0
    for j in tempFilesList:
        if "refined_up_" in j:
            upSpokeFileLocation=input_Path+"/"+i+"/"+i.replace("refined_","")+"/"+j
            spokeType="up"
            upSpokeData=readXML_data(upSpokeFileLocation,spokeType,srepNumber,file=i.replace("refined_",""))
            k=k+1
        elif "refined_down_" in j:
            downSpokeFileLocation=input_Path+"/"+i+"/"+i.replace("refined_","")+"/"+j
            spokeType="down"
            downSpokeData=readXML_data(downSpokeFileLocation,spokeType,srepNumber,file=i.replace("refined_",""))
            k=k+1
        elif "refined_crest_" in j:
            crestSpokeFileLocation=input_Path+"/"+i+"/"+i.replace("refined_","")+"/"+j
            spokeType="crest"
            crestSpokeData=readXML_data(crestSpokeFileLocation,spokeType,srepNumber,file=i.replace("refined_",""))
            k=k+1
    if(k!=3):
        print("Error!!!!!")
        print("The files are not encoded by slicer! Or There is a folder without up or down or crest spokes data!")
        break
    
    srepNumber=srepNumber+1
    
    allSreps.append([upSpokeData,downSpokeData,crestSpokeData])
    
    
    
dataFrameList=[]
for i in range(nSamples):
    for j in range(3):
        d= {"file":allSreps[i][j]["file"],"srepNumber":allSreps[i][j]["srepNumber"],"SpokesNumber":allSreps[i][j]["SpokesNumber"],
            "Spoke":allSreps[i][j]["Spoke"],
            "xPos":allSreps[i][j]["xPos"],"yPos":allSreps[i][j]["yPos"],"zPos":allSreps[i][j]["zPos"],
            "radii":allSreps[i][j]["Radii"],
            "Ux":allSreps[i][j]["Ux"],"Uy":allSreps[i][j]["Uy"],"Uz":allSreps[i][j]["Uz"]}
        df=pd.DataFrame(data=d)
        dataFrameList.append(df)
    
    
new_df = pd.concat(dataFrameList)

new_df.to_csv("output2.csv", sep=',',index=False)  # output2.csv is the name of the output file

print("Done!")