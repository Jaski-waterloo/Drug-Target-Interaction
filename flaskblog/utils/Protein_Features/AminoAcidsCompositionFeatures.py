#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re


# In[2]:


AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]


# In[3]:


def CalculateAAComposition(ProteinSequence):
    LengthSequence=len(ProteinSequence)
    Result={}
    for i in AALetter:
        Result[i]=round(float(ProteinSequence.count(i))/LengthSequence*100,3)
    return Result


# In[4]:


def CalculateDipeptideComposition(ProteinSequence):
    LengthSequence=len(ProteinSequence)
    Result={}
    for i in AALetter:
        for j in AALetter:
            Dipeptide=i+j
            Result[Dipeptide]=round(float(ProteinSequence.count(Dipeptide))/(LengthSequence-1)*100,2)
    return Result


# In[5]:


def Getkmers():
    kmers=list()
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i+j+k)
    return kmers


# In[6]:


def GetSpectrumDict(proteinsequence):
    result={}
    kmers=Getkmers()
    for i in kmers:
        result[i]=len(re.findall(i,proteinsequence))
    return result


# In[7]:


def CalculateAADipeptideComposition(ProteinSequence):
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    result={}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))
    result.update(GetSpectrumDict(ProteinSequence))
    
    return result


# In[ ]:




