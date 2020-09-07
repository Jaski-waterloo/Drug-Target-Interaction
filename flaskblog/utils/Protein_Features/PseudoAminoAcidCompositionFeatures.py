#!/usr/bin/env python
# coding: utf-8

# In[1]:


import string
import math


# In[2]:


AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={"A":0.62,"R":-2.53,"N":-0.78,"D":-0.90,"C":0.29,"Q":-0.85,"E":-0.74,"G":0.48,"H":-0.40,"I":1.38,"L":1.06,"K":-1.50,"M":0.64,"F":1.19,"P":0.12,"S":-0.18,"T":-0.05,"W":0.81,"Y":0.26,"V":1.08}

_hydrophilicity={"A":-0.5,"R":3.0,"N":0.2,"D":3.0,"C":-1.0,"Q":0.2,"E":3.0,"G":0.0,"H":-0.5,"I":-1.8,"L":-1.8,"K":3.0,"M":-1.3,"F":-2.5,"P":0.0,"S":0.3,"T":-0.4,"W":-3.4,"Y":-2.3,"V":-1.5}

_residuemass={"A":15.0,"R":101.0,"N":58.0,"D":59.0,"C":47.0,"Q":72.0,"E":73.0,"G":1.000,"H":82.0,"I":57.0,"L":57.0,"K":73.0,"M":75.0,"F":91.0,"P":42.0,"S":31.0,"T":45.0,"W":130.0,"Y":107.0,"V":43.0}

_pK1={"A":2.35,"C":1.71,"D":1.88,"E":2.19,"F":2.58,"G":2.34,"H":1.78,"I":2.32,"K":2.20,"L":2.36,"M":2.28,"N":2.18,"P":1.99,"Q":2.17,"R":2.18,"S":2.21,"T":2.15,"V":2.29,"W":2.38,"Y":2.20}

_pK2={"A":9.87,"C":10.78,"D":9.60,"E":9.67,"F":9.24,"G":9.60,"H":8.97,"I":9.76,"K":8.90,"L":9.60,"M":9.21,"N":9.09,"P":10.6,"Q":9.13,"R":9.09,"S":9.15,"T":9.12,"V":9.74,"W":9.39,"Y":9.11}

_pI={"A":6.11,"C":5.02,"D":2.98,"E":3.08,"F":5.91,"G":6.06,"H":7.64,"I":6.04,"K":9.47,"L":6.04,"M":5.74,"N":10.76,"P":6.30,"Q":5.65,"R":10.76,"S":5.68,"T":5.60,"V":6.02,"W":5.88,"Y":5.63}


# In[3]:


def _mean(listvalue):
    return sum(listvalue)/len(listvalue)


# In[4]:


def _std(listvalue,ddof=1):
    mean=_mean(listvalue)
    temp=[math.pow(i-mean,2) for i in listvalue]
    res=math.sqrt(sum(temp)/(len(listvalue)-ddof))
    return res


# In[5]:


def NormalizeEachAAP(AAP):
    if len(AAP.values())!=20:
        print ('You can not input the correct number of properities of Amino acids!')
    else:
        Result={}
        for i,j in AAP.items():
            Result[i]=(j-_mean(AAP.values()))/_std(AAP.values(),ddof=0)

    return Result


# In[6]:


#TYPE 1


# In[7]:


def _GetCorrelationFunction(Ri='S',Rj='D',AAP=[_Hydrophobicity,_hydrophilicity,_residuemass]):
    Hydrophobicity=NormalizeEachAAP(AAP[0])
    hydrophilicity=NormalizeEachAAP(AAP[1])
    residuemass=NormalizeEachAAP(AAP[2])
    theta1=math.pow(Hydrophobicity[Ri]-Hydrophobicity[Rj],2)
    theta2=math.pow(hydrophilicity[Ri]-hydrophilicity[Rj],2)
    theta3=math.pow(residuemass[Ri]-residuemass[Rj],2)
    theta=round((theta1+theta2+theta3)/3.0,3)
    return theta


# In[8]:


def _GetSequenceOrderCorrelationFactor(ProteinSequence,k=1):
    LengthSequence=len(ProteinSequence)
    res=[]
    for i in range(LengthSequence-k):
        AA1=ProteinSequence[i]
        AA2=ProteinSequence[i+k]
        res.append(_GetCorrelationFunction(AA1,AA2))
    result=round(sum(res)/(LengthSequence-k),3)
    return result


# In[9]:


def GetAAComposition(ProteinSequence):
    LengthSequence=len(ProteinSequence)
    Result={}
    for i in AALetter:
        Result[i]=round(float(ProteinSequence.count(i))/LengthSequence*100,3)
    return Result


# In[10]:


def _GetPseudoAAC1(ProteinSequence,lamda=10,weight=0.05):
    rightpart=0.0
    for i in range(lamda):
        rightpart=rightpart+_GetSequenceOrderCorrelationFactor(ProteinSequence,k=i+1)
    AAC=GetAAComposition(ProteinSequence)
    
    result={}
    temp=1+weight*rightpart
    for index,i in enumerate(AALetter):
        result['PAAC'+str(index+1)]=round(AAC[i]/temp,3)
    
    return result


# In[11]:


def _GetPseudoAAC2(ProteinSequence,lamda=10,weight=0.05):
    rightpart=[]
    for i in range(lamda):
        rightpart.append(_GetSequenceOrderCorrelationFactor(ProteinSequence,k=i+1))
    
    result={}
    temp=1+weight*sum(rightpart)
    for index in range(20,20+lamda):
        result['PAAC'+str(index+1)]=round(weight*rightpart[index-20]/temp*100,3)
    
    return result


# In[12]:


def _GetPseudoAAC(ProteinSequence,lamda=10,weight=0.05):
    res={}
    res.update(_GetPseudoAAC1(ProteinSequence,lamda=lamda,weight=weight))
    res.update(_GetPseudoAAC2(ProteinSequence,lamda=lamda,weight=weight))
    return res


# In[13]:


#TYPE 2


# In[14]:


def _GetCorrelationFunctionForAPAAC(Ri='S',Rj='D',AAP=[_Hydrophobicity,_hydrophilicity]):
    Hydrophobicity=NormalizeEachAAP(AAP[0])
    hydrophilicity=NormalizeEachAAP(AAP[1])
    theta1=round(Hydrophobicity[Ri]*Hydrophobicity[Rj],3)
    theta2=round(hydrophilicity[Ri]*hydrophilicity[Rj],3)

    return theta1,theta2


# In[15]:


def GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence,k=1):
    LengthSequence=len(ProteinSequence)
    resHydrophobicity=[]
    reshydrophilicity=[]
    for i in range(LengthSequence-k):
        AA1=ProteinSequence[i]
        AA2=ProteinSequence[i+k]
        temp=_GetCorrelationFunctionForAPAAC(AA1,AA2)
        resHydrophobicity.append(temp[0])
        reshydrophilicity.append(temp[1])
    result=[]
    result.append(round(sum(resHydrophobicity)/(LengthSequence-k),3))
    result.append(round(sum(reshydrophilicity)/(LengthSequence-k),3))
    return result


# In[16]:


def GetAPseudoAAC1(ProteinSequence,lamda=30,weight=0.5):
    rightpart=0.0
    for i in range(lamda):
        rightpart=rightpart+sum(GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence,k=i+1))
    AAC=GetAAComposition(ProteinSequence)
    
    result={}
    temp=1+weight*rightpart
    for index,i in enumerate(AALetter):
        result['APAAC'+str(index+1)]=round(AAC[i]/temp,3)
    
    return result


# In[17]:


def GetAPseudoAAC2(ProteinSequence,lamda=30,weight=0.5):
    rightpart=[]
    for i in range(lamda):
        temp=GetSequenceOrderCorrelationFactorForAPAAC(ProteinSequence,k=i+1)
        rightpart.append(temp[0])
        rightpart.append(temp[1])
        
    
    result={}
    temp=1+weight*sum(rightpart)
    for index in range(20,20+2*lamda):
        result['PAAC'+str(index+1)]=round(weight*rightpart[index-20]/temp*100,3)
    
    return result


# In[18]:


def GetAPseudoAAC(ProteinSequence,lamda=30,weight=0.5):
    res={}
    res.update(GetAPseudoAAC1(ProteinSequence,lamda=lamda,weight=weight))
    res.update(GetAPseudoAAC2(ProteinSequence,lamda=lamda,weight=weight))
    return res


# In[19]:


#TYPE I


# In[20]:


def GetCorrelationFunction(Ri='S',Rj='D',AAP=[]):
    NumAAP=len(AAP)
    theta=0.0
    for i in range(NumAAP):
        temp=NormalizeEachAAP(AAP[i])
        theta=theta+math.pow(temp[Ri]-temp[Rj],2)
    result=round(theta/NumAAP,3)
    return result


# In[21]:


def GetSequenceOrderCorrelationFactor(ProteinSequence,k=1,AAP=[]):
    LengthSequence=len(ProteinSequence)
    res=[]
    for i in range(LengthSequence-k):
        AA1=ProteinSequence[i]
        AA2=ProteinSequence[i+k]
        res.append(GetCorrelationFunction(AA1,AA2,AAP))
    result=round(sum(res)/(LengthSequence-k),3)
    return result


# In[22]:


def GetPseudoAAC1(ProteinSequence,lamda=30,weight=0.05,AAP=[]):
    rightpart=0.0
    for i in range(lamda):
        rightpart=rightpart+GetSequenceOrderCorrelationFactor(ProteinSequence,i+1,AAP)
    AAC=GetAAComposition(ProteinSequence)
    
    result={}
    temp=1+weight*rightpart
    for index,i in enumerate(AALetter):
        result['PAAC'+str(index+1)]=round(AAC[i]/temp,3)
    
    return result


# In[23]:


def GetPseudoAAC2(ProteinSequence,lamda=30,weight=0.05,AAP=[]):
    rightpart=[]
    for i in range(lamda):
        rightpart.append(GetSequenceOrderCorrelationFactor(ProteinSequence,i+1,AAP))
    
    result={}
    temp=1+weight*sum(rightpart)
    for index in range(20,20+lamda):
        result['PAAC'+str(index+1)]=round(weight*rightpart[index-20]/temp*100,3)
    
    return result


# In[24]:


def GetPseudoAAC(ProteinSequence,lamda=30,weight=0.05):
    
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    _Hydrophobicity={"A":0.62,"R":-2.53,"N":-0.78,"D":-0.90,"C":0.29,"Q":-0.85,"E":-0.74,"G":0.48,"H":-0.40,"I":1.38,"L":1.06,"K":-1.50,"M":0.64,"F":1.19,"P":0.12,"S":-0.18,"T":-0.05,"W":0.81,"Y":0.26,"V":1.08}

    _hydrophilicity={"A":-0.5,"R":3.0,"N":0.2,"D":3.0,"C":-1.0,"Q":0.2,"E":3.0,"G":0.0,"H":-0.5,"I":-1.8,"L":-1.8,"K":3.0,"M":-1.3,"F":-2.5,"P":0.0,"S":0.3,"T":-0.4,"W":-3.4,"Y":-2.3,"V":-1.5}

    _residuemass={"A":15.0,"R":101.0,"N":58.0,"D":59.0,"C":47.0,"Q":72.0,"E":73.0,"G":1.000,"H":82.0,"I":57.0,"L":57.0,"K":73.0,"M":75.0,"F":91.0,"P":42.0,"S":31.0,"T":45.0,"W":130.0,"Y":107.0,"V":43.0}

    _pK1={"A":2.35,"C":1.71,"D":1.88,"E":2.19,"F":2.58,"G":2.34,"H":1.78,"I":2.32,"K":2.20,"L":2.36,"M":2.28,"N":2.18,"P":1.99,"Q":2.17,"R":2.18,"S":2.21,"T":2.15,"V":2.29,"W":2.38,"Y":2.20}

    _pK2={"A":9.87,"C":10.78,"D":9.60,"E":9.67,"F":9.24,"G":9.60,"H":8.97,"I":9.76,"K":8.90,"L":9.60,"M":9.21,"N":9.09,"P":10.6,"Q":9.13,"R":9.09,"S":9.15,"T":9.12,"V":9.74,"W":9.39,"Y":9.11}

    _pI={"A":6.11,"C":5.02,"D":2.98,"E":3.08,"F":5.91,"G":6.06,"H":7.64,"I":6.04,"K":9.47,"L":6.04,"M":5.74,"N":10.76,"P":6.30,"Q":5.65,"R":10.76,"S":5.68,"T":5.60,"V":6.02,"W":5.88,"Y":5.63}
    AAP=[_Hydrophobicity,_hydrophilicity, _pK1, _pK2, _pI]
    res={}
    res.update(GetPseudoAAC1(ProteinSequence,lamda,weight,AAP))
    res.update(GetPseudoAAC2(ProteinSequence,lamda,weight,AAP))
    return res


# In[25]:


def GetPseudoAminoAcidCompositionFeatures(ProteinSequence):
    result = {}
    result.update(_GetPseudoAAC(ProteinSequence))
    result.update(GetAPseudoAAC(ProteinSequence))
    result.update(GetPseudoAAC(ProteinSequence, AAP=[_Hydrophobicity,_hydrophilicity, _pK1, _pK2, _pI]))
    
    return result


# In[ ]:




