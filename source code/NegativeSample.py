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
import Tool

def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:  
        SaveList.append(row)
    return

def ReadMyCsv2(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        counter = 0
        while counter < len(row):
            row[counter] = int(row[counter])      
            counter = counter + 1
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return




MDCuiMiDiseaseNum = []
ReadMyCsv(MDCuiMiDiseaseNum, "MDCuiMiDiseaseNum.csv")
print('MDCuiMiDiseaseNum[0]', MDCuiMiDiseaseNum[0])
print('len(MDCuiMiDiseaseNum)', len(MDCuiMiDiseaseNum))

AllDiseaseNum = []
ReadMyCsv(AllDiseaseNum, "AllDiseaseNum.csv")         

AllMiNum = []
ReadMyCsv(AllMiNum, "AllMiNum.csv")         

print(AllMiNum[0])

import Tool
NegativeSample = Tool.NegativeGenerate(MDCuiMiDiseaseNum, AllDiseaseNum, AllMiNum)
Tool.StorFile(NegativeSample, 'NegativeSample.csv')


