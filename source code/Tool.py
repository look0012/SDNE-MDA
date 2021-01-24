import numpy as np
from sklearn.decomposition import NMF
from pylab import *
import matplotlib.pyplot as plt
import csv
import random
from scipy.sparse.linalg import svds
from scipy import sparse
from numpy import *

def ReadMyCsv(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        for i in range(len(row)):
            row[i] = float(row[i])
        SaveList.append(row)
    return

def ReadMyCsv2(SaveList, fileName):
    csv_reader = csv.reader(open(fileName))
    for row in csv_reader:
        SaveList.append(row)
    return

def StorFile(data, fileName):
    with open(fileName, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)
    return

def GenerateAllDrug(LncRNADiseaseAssociation):

    AllDrug = []
    counter1 = 0
    while counter1 < len(LncRNADiseaseAssociation):
        counter2 = 0
        flag = 0
        while counter2 < len(AllDrug):
            if LncRNADiseaseAssociation[counter1][0] != AllDrug[counter2][0]:
                counter2 = counter2 + 1
            elif LncRNADiseaseAssociation[counter1][0] == AllDrug[counter2][0]:
                flag = 1
                break
        if flag == 0:
            pair = []
            pair.append(LncRNADiseaseAssociation[counter1][0])
            AllDrug.append(pair)
        counter1 = counter1 + 1
        print(counter1)
    return AllDrug
def GenerateAllDisease(LncRNADiseaseAssociation):
    AllDisease = []
    counter1 = 0
    while counter1 < len(LncRNADiseaseAssociation): 
        counter2 = 0
        flag = 0
        while counter2 < len(AllDisease):  
            if LncRNADiseaseAssociation[counter1][1] != AllDisease[counter2][0]:  
                counter2 = counter2 + 1
            elif LncRNADiseaseAssociation[counter1][1] == AllDisease[counter2][0]:  
                flag = 1
                counter2 = counter2 + 1
        if flag == 0:
            pair = []
            pair.append(LncRNADiseaseAssociation[counter1][1])
            AllDisease.append(pair)
        counter1 = counter1 + 1
        print(counter1)
    return AllDisease
def LowerData(Data):
    counter = 0
    while counter < len(Data):
        Data[counter][0] = Data[counter][0].lower()
        Data[counter][1] = Data[counter][1].lower()
        counter = counter + 1
    return Data
def MySampleLabel(num):
    SampleLabel = []
    counter = 0
    while counter < num:
        SampleLabel.append(1)
        counter = counter + 1
    counter1 = 0
    while counter1 < num:
        SampleLabel.append(0)
        counter1 = counter1 + 1
    return SampleLabel
def MyAssociationMatix(AssociationMatrix, TrainList, AllDrug, AllDisease):
    counter = 0
    while counter < len(TrainList):
        Drug = TrainList[counter][0]
        Disease = TrainList[counter][1]

        flag1 = 0
        counter1 = 0
        while counter1 < len(AllDrug):
            if Drug == AllDrug[counter1]:
                flag1 = 1
                break
            counter1 = counter1 + 1

        flag2 = 0
        counter2 = 0
        while counter2 < len(AllDisease):
            if Disease == AllDisease[counter2]:
                flag2 = 1
                break
            counter2 = counter2 + 1

        if flag1 == 1 and flag2 == 1:
            AssociationMatrix[counter1][counter2] = 1

        print(counter)
        counter = counter + 1
    return AssociationMatrix
def PositiveGenerate(DiseaseAndRNABinaryOld, DGS, RGS):
    PositiveFeature = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 1:
                pair = []
                pair.extend(RGS[counter1])
                pair.extend(DGS[counter])
                PositiveFeature.append(pair)
            counter1 = counter1 + 1
        counter = counter + 1
    return PositiveFeature
def GenerateSampleFeature(Sample, LncFeature, DiseaseFeature):
    SampleFeature = []
    counter = 0
    while counter < len(Sample):
        FeaturePair = []
        lnc = Sample[counter][0]
        disease = Sample[counter][1]
        counter1 = 0
        while counter1 < len(LncFeature):
            if lnc == LncFeature[counter1][0]:
                FeaturePair.extend(LncFeature[counter1][1:])
                break
            counter1 = counter1 + 1

        counter2 = 0
        while counter2 < len(DiseaseFeature):
            if disease == DiseaseFeature[counter2][0]:
                FeaturePair.extend(DiseaseFeature[counter2][1:])
                break
            counter2 = counter2 + 1
        SampleFeature.append(FeaturePair)
        counter = counter + 1
    return SampleFeature
def TTFeature(AllSampleFeature, MyList):
    SampleFeature = []
    counter = 0
    while counter < len(MyList):
        SampleFeature.append(AllSampleFeature[MyList[counter]])
        counter = counter + 1
    return SampleFeature
def NegativeCandidateGenerate(DiseaseAndRNABinaryOld, RNAGaussianOld, DiseaseGaussianOld):
    NegativeFeatureAll = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 0:
                pair = []
                pair.extend(RNAGaussianOld[counter1])
                pair.extend(DiseaseGaussianOld[counter])
                NegativeFeatureAll.append(pair)
            counter1 = counter1 + 1
        counter = counter + 1
        print(counter)
    return NegativeFeatureAll

def NegativeGenerate(LncDisease, AllRNA,AllDisease):

    import random
    NegativeSample = []
    counterN = 0
    while counterN < len(LncDisease):
        counterR = random.randint(0, len(AllRNA) - 1)
        counterD = random.randint(0, len(AllDisease) - 1)
        DiseaseAndRnaPair = []
        DiseaseAndRnaPair.append(AllRNA[counterR])
        DiseaseAndRnaPair.append(AllDisease[counterD])
        flag1 = 0
        counter = 0
        while counter < len(LncDisease):
            if DiseaseAndRnaPair == LncDisease[counter]:
                flag1 = 1
                break
            counter = counter + 1
        if flag1 == 1:
            continue
        flag2 = 0
        counter1 = 0
        while counter1 < len(NegativeSample):
            if DiseaseAndRnaPair == NegativeSample[counter1]:
                flag2 = 1
                break
            counter1 = counter1 + 1
        if flag2 == 1:
            continue
        if (flag1 == 0 & flag2 == 0):
            NamePair = []
            NamePair.append(AllRNA[counterR])
            NamePair.append(AllDisease[counterD])
            NegativeSample.append(NamePair)

            counterN = counterN + 1
        print(counterN)
    return NegativeSample

def NegativeGenerate2(RNAFeatureDAG, DiseaseFeatureDAG,RNAFeatureNMFDAG, DiseaseFeatureNMFDAG):

    LncDisease = []
    ReadMyCsv2(LncDisease, 'LncDisease.csv')
    AllDisease = []
    ReadMyCsv2(AllDisease, 'AllDisease.csv')
    AllRNA = []
    ReadMyCsv2(AllRNA, 'AllRNA.csv')
    import random
    NegativeSample = []
    NegativeSampleFeature = []
    NegativeSampleFeatureNMF = []
    counterN = 0
    while counterN < len(LncDisease):
        counterD = random.randint(0, len(AllDisease) - 1)
        counterR = random.randint(0, len(AllRNA) - 1)
        DiseaseAndRnaPair = []
        DiseaseAndRnaPair.append(AllRNA[counterR])
        DiseaseAndRnaPair.append(AllDisease[counterD])
        flag1 = 0
        counter = 0
        while counter < len(LncDisease):
            if DiseaseAndRnaPair == LncDisease[counter]:
                flag1 = 1
                break
            counter = counter + 1
        if flag1 == 1:
            continue
        flag2 = 0
        counter1 = 0
        while counter1 < len(NegativeSample):
            if DiseaseAndRnaPair == NegativeSample[counter1]:
                flag2 = 1
                break
            counter1 = counter1 + 1
        if flag2 == 1:
            continue
        if (flag1 == 0 & flag2 == 0):
            NamePair = []
            NamePair.append(AllRNA[counterR][0])
            NamePair.append(AllDisease[counterD][0])
            NegativeSample.append(NamePair)

            FeaturePair0 = []
            FeaturePair0.extend(RNAFeatureDAG[counterR])
            FeaturePair0.extend(DiseaseFeatureDAG[counterD])
            NegativeSampleFeature.append(FeaturePair0)

            FeaturePair1 = []
            FeaturePair1.extend(RNAFeatureNMFDAG[counterR])
            FeaturePair1.extend(DiseaseFeatureNMFDAG[counterD])
            NegativeSampleFeatureNMF.append(FeaturePair0)

            counterN = counterN + 1
    return NegativeSampleFeature, NegativeSampleFeatureNMF, NegativeSample

def StrongNegativeGenerate(DiseaseAndRNABinaryOld, RNAGaussianOld, DiseaseGaussianOld,LncRNADiseaseAssociationOld):

    NegativeFeatureAll = []
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld[counter]):
            if DiseaseAndRNABinaryOld[counter][counter1] == 0:
                pairFeature = []
                pairFeature.extend(RNAGaussianOld[counter1])
                pairFeature.extend(DiseaseGaussianOld[counter])
                NegativeFeatureAll.append(pairFeature)
            counter1 = counter1 + 1
        counter = counter + 1
    from sklearn.ensemble import IsolationForest
    clf = IsolationForest(contamination=0.1)
    clf.fit(NegativeFeatureAll)
    scores_pred = clf.decision_function(NegativeFeatureAll)

    PredictionScoreNum = []
    counter = 0
    while counter < len(scores_pred):
        pair = []
        pair.append(scores_pred[counter])
        pair.append(counter)
        PredictionScoreNum.append(pair)
        counter = counter + 1
    SerialNumber = 0
    MaxScoreNum = []
    counter = 0
    while counter < len(LncRNADiseaseAssociationOld):
        max = PredictionScoreNum[0][0]
        counter1 = 0
        while counter1 < len(PredictionScoreNum):
            if max < PredictionScoreNum[counter1][0]:
                max = PredictionScoreNum[counter1][0]
                SerialNumber = counter1
            counter1 = counter1 + 1
        MaxScoreNum.append(PredictionScoreNum[SerialNumber][1])
        del PredictionScoreNum[SerialNumber]
        counter = counter + 1
    NegativeFeature = []
    counter = 0
    while counter < len(MaxScoreNum):
        NegativeFeature.append(NegativeFeatureAll[MaxScoreNum[counter]])
        counter = counter + 1
    return NegativeFeature, NegativeFeatureAll
def MyEvaluate(prediction, prediction_proba,TestSample):
    import math
    num = 0
    SumRMSE = 0
    SumMAE = 0
    counter = 0
    while counter < len(prediction):
        SumRMSE = SumRMSE + math.pow((1 - prediction_proba[counter][1]), 2)
        SumMAE = SumMAE + abs(1 - prediction_proba[counter][1])
        if prediction[counter] == 1:
            num = num + 1
        counter = counter + 1
    RMSE = math.sqrt(SumRMSE / len(TestSample))
    MAE = SumMAE / len(TestSample)
    print('TrueNum ?/243: ', num)
    print('RMSE:', RMSE)
    print('MAE:', MAE)
    MyResult = []
    MyResult.append(num)
    MyResult.append(RMSE)
    MyResult.append(MAE)
    return MyResult

def MyPrediction(SampleFeature,SampleLabel,TestSample):
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=100)
    model.fit(SampleFeature, SampleLabel)
    prediction = model.predict(TestSample)
    prediction_proba = model.predict_proba(TestSample)
    print('RandomForestClassifier!')
    result = MyEvaluate(prediction, prediction_proba, TestSample)
    return result
def MyPredictionAndMatrixCompletion(SampleFeature,SampleLabel,NegativeFeatureAll,DiseaseAndRNABinaryOld1,DiseaseAndRNABinaryOld2,TestSample):
    from sklearn.ensemble import RandomForestClassifier
    model = RandomForestClassifier(n_estimators=100)
    model.fit(SampleFeature, SampleLabel)
    prediction = model.predict(TestSample)
    prediction_proba = model.predict_proba(TestSample)
    print('RandomForestClassifier!')
    result = MyEvaluate(prediction, prediction_proba, TestSample)

    prediction_proba_all = model.predict_proba(NegativeFeatureAll)
    num = 0
    counter = 0
    while counter < len(DiseaseAndRNABinaryOld1):
        counter1 = 0
        while counter1 < len(DiseaseAndRNABinaryOld1[counter]):
            if DiseaseAndRNABinaryOld1[counter][counter1] == 0:
                DiseaseAndRNABinaryOld2[counter][counter1] = prediction_proba_all[num][1]
                num = num + 1
            counter1 = counter1 + 1
        counter = counter + 1
    return DiseaseAndRNABinaryOld2, result

def LncRNAGaussianKernel(DiseaseAndRNABinary):

    counter1 = 0
    sum1 = 0
    while counter1 < (len(DiseaseAndRNABinary)):
        counter2 = 0
        while counter2 < (len(DiseaseAndRNABinary[counter1])):
            sum1 = sum1 + pow((DiseaseAndRNABinary[counter1][counter2]), 2)
            counter2 = counter2 + 1
        counter1 = counter1 + 1

    Ak = sum1
    Nd = len(DiseaseAndRNABinary)
    rdpie = 0.5
    rd = rdpie * Nd / Ak
    DiseaseGaussian = []
    counter1 = 0
    while counter1 < len(DiseaseAndRNABinary):
        counter2 = 0
        DiseaseGaussianRow = []
        while counter2 < len(DiseaseAndRNABinary):
            AiMinusBj = 0
            sum2 = 0
            counter3 = 0
            AsimilarityB = 0
            while counter3 < len(DiseaseAndRNABinary[counter2]):
                sum2 = pow((DiseaseAndRNABinary[counter1][counter3] - DiseaseAndRNABinary[counter2][counter3]),2)  
                AiMinusBj = AiMinusBj + sum2
                counter3 = counter3 + 1
            AsimilarityB = math.exp(- (AiMinusBj / rd))
            DiseaseGaussianRow.append(AsimilarityB)
            counter2 = counter2 + 1
        DiseaseGaussian.append(DiseaseGaussianRow)
        counter1 = counter1 + 1
        print(counter1)
    return DiseaseGaussian

def DiseaseGaussianKernel(DiseaseAndRNABinary):
    MDiseaseAndRNABinary = np.array(DiseaseAndRNABinary)
    RNAAndDiseaseBinary = MDiseaseAndRNABinary.T
    RNAGaussian = []
    counter1 = 0
    sum1 = 0
    while counter1 < (len(RNAAndDiseaseBinary)):
        counter2 = 0
        while counter2 < (len(RNAAndDiseaseBinary[counter1])):
            sum1 = sum1 + pow((RNAAndDiseaseBinary[counter1][counter2]), 2)
            counter2 = counter2 + 1
        counter1 = counter1 + 1
    Ak = sum1
    Nm = len(RNAAndDiseaseBinary)
    rdpie = 0.5
    rd = rdpie * Nm / Ak
    counter1 = 0
    while counter1 < len(RNAAndDiseaseBinary):
        counter2 = 0
        RNAGaussianRow = []
        while counter2 < len(RNAAndDiseaseBinary):
            AiMinusBj = 0
            sum2 = 0
            counter3 = 0
            AsimilarityB = 0
            while counter3 < len(RNAAndDiseaseBinary[counter2]):
                sum2 = pow((RNAAndDiseaseBinary[counter1][counter3] - RNAAndDiseaseBinary[counter2][counter3]),2)  
                AiMinusBj = AiMinusBj + sum2
                counter3 = counter3 + 1
            AsimilarityB = math.exp(- (AiMinusBj / rd))
            RNAGaussianRow.append(AsimilarityB)
            counter2 = counter2 + 1
        RNAGaussian.append(RNAGaussianRow)
        counter1 = counter1 + 1
        print(counter1)
    return RNAGaussian

def TestSampleFeatureGenerate(LncRNADiseaseAssociationNew, AllDiseaseOld, AllRNAOld, RNAGaussianOld, DiseaseGaussianOld, DiseaseAndRNABinaryOld):

    ExtraPairNum = []
    ExtraPairName = []
    TestSampleFeature = []
    counter = 0
    while counter < len(LncRNADiseaseAssociationNew):
        rna = LncRNADiseaseAssociationNew[counter][0]
        disease = LncRNADiseaseAssociationNew[counter][1]
        counter1 = 0
        while counter1 < len(AllDiseaseOld):
            if disease == AllDiseaseOld[counter1]:
                counter2 = 0
                while counter2 < len(AllRNAOld):
                    if rna == AllRNAOld[counter2]:
                        if DiseaseAndRNABinaryOld[counter1][counter2] == 0:
                            pairNum = []
                            pairNum.append(counter1)
                            pairNum.append(counter2)
                            ExtraPairNum.append(pairNum)
                            pairName = []
                            pairName.append(AllDiseaseOld[counter1])
                            pairName.append(AllRNAOld[counter2])
                            ExtraPairName.append(pairName)
                            pairFeature = []
                            pairFeature.extend(RNAGaussianOld[counter2])
                            pairFeature.extend(DiseaseGaussianOld[counter1])
                            TestSampleFeature.append(pairFeature)
                        break
                    counter2 = counter2 + 1
                break
            counter1 = counter1 + 1
        counter = counter + 1
    return TestSampleFeature