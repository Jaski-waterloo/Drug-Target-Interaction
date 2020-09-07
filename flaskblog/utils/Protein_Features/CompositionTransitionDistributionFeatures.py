#!/usr/bin/env python
# coding: utf-8

# In[1]:


import string, math, copy


# In[2]:


AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'} 
#'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_NormalizedVDWV={'1':'GASTPD','2':'NVEQIL','3':'MHKFRYW'}
#'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
#'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'}
#'1'stand for Helix; '2'stand for Strand, '3' stand for coil

_SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
#'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

_Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
#'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)


# In[3]:


_AATProperty=(_Hydrophobicity,_NormalizedVDWV,_Polarity,_Charge,_SecondaryStr,_SolventAccessibility,_Polarizability)

_AATPropertyName=('_Hydrophobicity','_NormalizedVDWV','_Polarity','_Charge','_SecondaryStr','_SolventAccessibility','_Polarizability')


# In[4]:


def StringtoNum(ProteinSequence,AAProperty):
    hardProteinSequence=copy.deepcopy(ProteinSequence)
    for k,m in AAProperty.items():
        for index in m:
            hardProteinSequence=str.replace(hardProteinSequence,index,k)
    TProteinSequence=hardProteinSequence

    return TProteinSequence


# In[5]:


def CalculateComposition(ProteinSequence,AAProperty,AAPName):
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    Result[AAPName+'C'+'1']=round(float(TProteinSequence.count('1'))/Num,3)
    Result[AAPName+'C'+'2']=round(float(TProteinSequence.count('2'))/Num,3)
    Result[AAPName+'C'+'3']=round(float(TProteinSequence.count('3'))/Num,3)
    return Result


# In[6]:


def CalculateTransition(ProteinSequence,AAProperty,AAPName):
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    CTD=TProteinSequence
    Result[AAPName+'T'+'12']=round(float(CTD.count('12')+CTD.count('21'))/(Num-1),3)
    Result[AAPName+'T'+'13']=round(float(CTD.count('13')+CTD.count('31'))/(Num-1),3)
    Result[AAPName+'T'+'23']=round(float(CTD.count('23')+CTD.count('32'))/(Num-1),3)
    return Result


# In[7]:


def CalculateDistribution(ProteinSequence,AAProperty,AAPName):
    TProteinSequence=StringtoNum(ProteinSequence,AAProperty)
    Result={}
    Num=len(TProteinSequence)
    temp=('1','2','3')
    for i in temp:
        num=TProteinSequence.count(i)
        ink=1
        indexk=0
        cds=[]
        while ink<=num:
            indexk=str.find(TProteinSequence,i,indexk)+1
            cds.append(indexk)
            ink=ink+1
                
        if cds==[]:
            Result[AAPName+'D'+i+'001']=0
            Result[AAPName+'D'+i+'025']=0
            Result[AAPName+'D'+i+'050']=0
            Result[AAPName+'D'+i+'075']=0
            Result[AAPName+'D'+i+'100']=0
        else:
                
            Result[AAPName+'D'+i+'001']=round(float(cds[0])/Num*100,3)
            Result[AAPName+'D'+i+'025']=round(float(cds[int(math.floor(num*0.25))-1])/Num*100,3)
            Result[AAPName+'D'+i+'050']=round(float(cds[int(math.floor(num*0.5))-1])/Num*100,3)
            Result[AAPName+'D'+i+'075']=round(float(cds[int(math.floor(num*0.75))-1])/Num*100,3)
            Result[AAPName+'D'+i+'100']=round(float(cds[-1])/Num*100,3)

    return Result


# In[8]:


def CalculateCompositionHydrophobicity(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[9]:


def CalculateCompositionNormalizedVDWV(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result


# In[10]:


def CalculateCompositionPolarity(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_Polarity,'_Polarity')
    return result


# In[11]:


def CalculateCompositionCharge(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_Charge,'_Charge')
    return result


# In[12]:


def CalculateCompositionSecondaryStr(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result


# In[13]:


def CalculateCompositionSolventAccessibility(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result


# In[14]:


def CalculateCompositionPolarizability(ProteinSequence):
    result=CalculateComposition(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[15]:


def CalculateTransitionHydrophobicity(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[16]:


def CalculateTransitionNormalizedVDWV(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result


# In[17]:


def CalculateTransitionPolarity(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_Polarity,'_Polarity')
    return result


# In[18]:


def CalculateTransitionCharge(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_Charge,'_Charge')
    return result


# In[19]:


def CalculateTransitionSecondaryStr(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result


# In[20]:


def CalculateTransitionSolventAccessibility(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result


# In[21]:


def CalculateTransitionPolarizability(ProteinSequence):
    result=CalculateTransition(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[22]:


def CalculateDistributionHydrophobicity(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[23]:


def CalculateDistributionNormalizedVDWV(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_NormalizedVDWV,'_NormalizedVDWV')
    return result


# In[24]:


def CalculateDistributionPolarity(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_Polarity,'_Polarity')
    return result


# In[25]:


def CalculateDistributionCharge(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_Charge,'_Charge')
    return result


# In[26]:


def CalculateDistributionSecondaryStr(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_SecondaryStr,'_SecondaryStr')
    return result


# In[27]:


def CalculateDistributionSolventAccessibility(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_SolventAccessibility,'_SolventAccessibility')
    return result


# In[28]:


def CalculateDistributionPolarizability(ProteinSequence):
    result=CalculateDistribution(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[29]:


def CalculateC(ProteinSequence):
    result={}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    #result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    return result


# In[30]:


def CalculateT(ProteinSequence):
    result={}
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    #result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    return result


# In[31]:


def CalculateD(ProteinSequence):
    result={}
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    #result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result


# In[32]:


def CalculateCompositionTransitionDistribution(ProteinSequence):
    AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    _Hydrophobicity={'1':'RKEDQN','2':'GASTPHY','3':'CLVIMFW'} 
    #'1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

    _NormalizedVDWV={'1':'GASTPD','2':'NVEQIL','3':'MHKFRYW'}
    #'1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

    _Polarity={'1':'LIFWCMVY','2':'CPNVEQIL','3':'KMHFRYW'}
    #'1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

    _Charge={'1':'KR','2':'ANCQGHILMFPSTWYV','3':'DE'}
    #'1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

    _SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'}
    #'1'stand for Helix; '2'stand for Strand, '3' stand for coil

    _SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND','3':'MPSTHY'}
    #'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

    _Polarizability={'1':'GASDT','2':'CPNVEQIL','3':'KMHFRYW'}
    #'1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)
    result={}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    #result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    #result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    #result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result


# In[ ]:




