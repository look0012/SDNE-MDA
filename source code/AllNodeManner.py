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
import copy

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

AllNode = []
ReadMyCsv(AllNode, 'AllNode.csv')
FinalDrugFeatureNum = []
ReadMyCsv(FinalDrugFeatureNum, 'FinalDiseaseFeatureNum.csv')
FinalDiseaseFeatureNum = []
ReadMyCsv(FinalDiseaseFeatureNum, 'FinalMiKmerNum.csv')
LineEmbeddingName1 = 'vec_all0.txt'
LineEmbedding1 = np.loadtxt(LineEmbeddingName1, dtype=str, skiprows=1)
StorFile(LineEmbedding1, 'LineEmbedding1.csv')
counterP = 0
while counterP < 5:
    LineEmbeddingName = 'vec_all' + str(counterP) + '.txt'
    LineEmbedding = np.loadtxt(LineEmbeddingName, dtype=str, skiprows=1)
    AllNodeMannerNum = []
    counter = 0
    while counter < len(AllNode):
        pair = []
        counter1 = 0
        while counter1 < len(LineEmbedding[0]) - 1:
            pair.append(0)
            counter1 = counter1 + 1
        AllNodeMannerNum.append(pair)
        counter = counter + 1
    counter = 0
    while counter < len(LineEmbedding):
        num = int(LineEmbedding[counter][0])
        AllNodeMannerNum[num] = list(LineEmbedding[counter][1:])
        counter = counter + 1
    print(np.array(AllNodeMannerNum).shape)
    AllNodeMannerNumName = 'AllNodeMannerNum' + str(counterP) + '.csv'

    num1 = 0
    counter = 0
    while counter < len(FinalDrugFeatureNum):
        AllNodeMannerNum[int(FinalDrugFeatureNum[counter][0])].extend(FinalDrugFeatureNum[counter][1:])
        num1 = num1 + 1
        counter = counter + 1
    print(num1)
    num2 = 0
    counter = 0
    while counter < len(FinalDiseaseFeatureNum):
        AllNodeMannerNum[int(FinalDiseaseFeatureNum[counter][0])].extend(FinalDiseaseFeatureNum[counter][1:])
        num2 = num2 + 1
        counter = counter + 1
    print(num2)
    print(np.array(AllNodeMannerNum).shape)
    AllNodeFeatureNumName = 'AllNodeAttributeMannerNum' + str(counterP) + '.csv'
    counterP = counterP + 1


