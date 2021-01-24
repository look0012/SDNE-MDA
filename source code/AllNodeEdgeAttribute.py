from numpy import *
import numpy as np
import random
import math
import os
import time
import pandas as pd
import csv
import math
import random

def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return


LDAllLncDisease = []
ReadMyCsv(LDAllLncDisease, "LDAllLncDisease.csv")
print('LDAllLncDisease[0]', LDAllLncDisease[0])
LMSNPLncMi = []
ReadMyCsv(LMSNPLncMi, "LMSNPLncMi.csv")
print('LMSNPLncMi[0]', LMSNPLncMi[0])
LPLncRNA2TargetLncProtein3 = []
ReadMyCsv(LPLncRNA2TargetLncProtein3, "LPLncRNA2TargetLncProtein3.csv")
print('LPLncRNA2TargetLncProtein3[0]', LPLncRNA2TargetLncProtein3[0])
MDCuiMiDisease = []
ReadMyCsv(MDCuiMiDisease, "MDCuiMiDisease.csv")
print('MDCuiMiDisease[0]', MDCuiMiDisease[0])
MPmiRTarBaseMiProtein9 = []
ReadMyCsv(MPmiRTarBaseMiProtein9, "MPmiRTarBaseMiProtein5.csv")
print('MPmiRTarBaseMiProtein9[0]', MPmiRTarBaseMiProtein9[0])
PDDisGeNETProteinDisease15 = []
ReadMyCsv(PDDisGeNETProteinDisease15, "PDDisGeNETProteinDisease20.csv")
print('PDDisGeNETProteinDisease15[0]', PDDisGeNETProteinDisease15[0])
PPI = []
ReadMyCsv(PPI, "PPI.csv")
print('PPI[0]', PPI[0])
AllEdge = []
AllEdge.extend(LDAllLncDisease)
AllEdge.extend(LMSNPLncMi)
AllEdge.extend(MPmiRTarBaseMiProtein9)
AllEdge.extend(LPLncRNA2TargetLncProtein3)
AllEdge.extend(PDDisGeNETProteinDisease15)
AllEdge.extend(PPI)
print(len(AllEdge))
print(AllEdge[0])
StorFile(AllEdge, 'AllEdge.csv')
FinalDiseaseFeature = []
ReadMyCsv(FinalDiseaseFeature, "FinalDiseaseFeature.csv")
FinalAllDisease = np.array(FinalDiseaseFeature)[:, 0]
print('len(FinalAllDisease)', len(FinalAllDisease))
print('FinalAllDisease[0]', FinalAllDisease[0])
FinalLncKmer = []
ReadMyCsv(FinalLncKmer, "FinalLncKmer.csv")
FinalLnc = np.array(FinalLncKmer)[:, 0]
print('len(FinalLnc)', len(FinalLnc))
print('FinalLnc[0]', FinalLnc[0])
FinalMiKmer = []
ReadMyCsv(FinalMiKmer, "FinalMiKmer.csv")
FinalMi = np.array(FinalMiKmer)
FinalMi = FinalMi[:, 0]
print('len(FinalMi)', len(FinalMi))
print('FinalMi[0]', FinalMi[0])
FinalProteinFeature = []
ReadMyCsv(FinalProteinFeature, "FinalProteinFeature.csv")
FinalProtein = np.array(FinalProteinFeature)[:, 0]
print('len(FinalProtein)', len(FinalProtein))
print('FinalProtein[0]', FinalProtein[0])
AllNode = []
AllNode.extend(FinalMi)
AllNode.extend(FinalLnc)
AllNode.extend(FinalProtein)
AllNode.extend(FinalAllDisease)
print('len(AllNode)', len(AllNode))
counter = 0
while counter < len(AllNode):
    pair = []
    pair.append(AllNode[counter])
    AllNode[counter] = pair
    counter = counter + 1
print(AllNode[0])
StorFile(AllNode, 'AllNode.csv')
AllNodeAttribute = []
AllNodeAttribute.extend(FinalMiKmer)
AllNodeAttribute.extend(FinalLncKmer)
AllNodeAttribute.extend(FinalProteinFeature)
AllNodeAttribute.extend(FinalDiseaseFeature)
StorFile(AllNodeAttribute, 'AllNodeAttribute.csv')
