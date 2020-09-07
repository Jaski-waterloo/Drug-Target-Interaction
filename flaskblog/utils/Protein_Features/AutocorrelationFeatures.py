#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math,string


# In[2]:


AALetter=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

_Hydrophobicity={"A":0.02,"R":-0.42,"N":-0.77,"D":-1.04,"C":0.77,"Q":-1.10,"E":-1.14,"G":-0.80,"H":0.26,"I":1.81,"L":1.14,"K":-0.41,"M":1.00,"F":1.35,"P":-0.09,"S":-0.97,"T":-0.77,"W":1.71,"Y":1.11,"V":1.13}

_AvFlexibility={"A":0.357,"R":0.529,"N":0.463,"D":0.511,"C":0.346,"Q":0.493,"E":0.497,"G":0.544,"H":0.323,"I":0.462,"L":0.365,"K":0.466,"M":0.295,"F":0.314,"P":0.509,"S":0.507,"T":0.444,"W":0.305,"Y":0.420,"V":0.386}

_Polarizability={"A":0.046,"R":0.291,"N":0.134,"D":0.105,"C":0.128,"Q":0.180,"E":0.151,"G":0.000,"H":0.230,"I":0.186,"L":0.186,"K":0.219,"M":0.221,"F":0.290,"P":0.131,"S":0.062,"T":0.108,"W":0.409,"Y":0.298,"V":0.140}

_FreeEnergy={"A":-0.368,"R":-1.03,"N":0.0,"D":2.06,"C":4.53,"Q":0.731,"E":1.77,"G":-0.525,"H":0.0,"I":0.791,"L":1.07,"K":0.0,"M":0.656,"F":1.06,"P":-2.24,"S":-0.524,"T":0.0,"W":1.60,"Y":4.91,"V":0.401}

_ResidueASA={"A":115.0,"R":225.0,"N":160.0,"D":150.0,"C":135.0,"Q":180.0,"E":190.0,"G":75.0,"H":195.0,"I":175.0,"L":170.0,"K":200.0,"M":185.0,"F":210.0,"P":145.0,"S":115.0,"T":140.0,"W":255.0,"Y":230.0,"V":155.0}

_ResidueVol={"A":52.6,"R":109.1,"N":75.7,"D":68.4,"C":68.3,"Q":89.7,"E":84.7,"G":36.3,"H":91.9,"I":102.0,"L":102.0,"K":105.1,"M":97.7,"F":113.9,"P":73.6,"S":54.9,"T":71.2,"W":135.4,"Y":116.2,"V":85.1}

_Steric={"A":0.52,"R":0.68,"N":0.76,"D":0.76,"C":0.62,"Q":0.68,"E":0.68,"G":0.00,"H":0.70,"I":1.02,"L":0.98,"K":0.68,"M":0.78,"F":0.70,"P":0.36,"S":0.53,"T":0.50,"W":0.70,"Y":0.70,"V":0.76}

_Mutability={"A":100.0,"R":65.0,"N":134.0,"D":106.0,"C":20.0,"Q":93.0,"E":102.0,"G":49.0,"H":66.0,"I":96.0,"L":40.0,"K":-56.0,"M":94.0,"F":41.0,"P":56.0,"S":120.0,"T":97.0,"W":18.0,"Y":41.0,"V":74.0}


# In[3]:


_AAProperty=(_Hydrophobicity,_AvFlexibility,_Polarizability,_FreeEnergy,_ResidueASA,_ResidueVol,_Steric,_Mutability)

_AAPropertyName=('_Hydrophobicity','_AvFlexibility','_Polarizability','_FreeEnergy','_ResidueASA','_ResidueVol','_Steric','_Mutability')


# In[4]:


def _mean(listvalue):
    """
    The mean value of the list data.
    """
    return sum(listvalue)/len(listvalue)


# In[5]:


def _std(listvalue,ddof=1):
    """
    The standard deviation of the list data.
    """
    mean=_mean(listvalue)
    temp=[math.pow(i-mean,2) for i in listvalue]
    res=math.sqrt(sum(temp)/(len(listvalue)-ddof))
    return res


# In[6]:


def NormalizeEachAAP(AAP):
    Result = {}
    if len(AAP.values())!=20:
        print('You can not input the correct number of properities of Amino acids!')
    else:
        for i,j in AAP.items():
            Result[i]=(j-_mean(AAP.values()))/_std(AAP.values(),ddof=0)

    return Result


# In[7]:


def CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,AAP,AAPName):
    AAPdic=NormalizeEachAAP(AAP)

    Result={}
    for i in range(1,31):
        temp=0
        for j in range(len(ProteinSequence)-i):
            temp=temp+AAPdic[ProteinSequence[j]]*AAPdic[ProteinSequence[j+1]]
        if len(ProteinSequence)-i==0:
            Result['MoreauBrotoAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)),3)
        else:
            Result['MoreauBrotoAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)-i),3)

    return Result


# In[8]:


def CalculateEachMoranAuto(ProteinSequence,AAP,AAPName):
    
    AAPdic=NormalizeEachAAP(AAP)
    cds=0
    for i in AALetter:
        cds=cds+(ProteinSequence.count(i))*(AAPdic[i])
    Pmean=cds/len(ProteinSequence)

    cc=[]
    for i in ProteinSequence:
        cc.append(AAPdic[i])

    K=(_std(cc,ddof=0))**2

    Result={}
    for i in range(1,31):
        temp=0
        for j in range(len(ProteinSequence)-i):
                
            temp=temp+(AAPdic[ProteinSequence[j]]-Pmean)*(AAPdic[ProteinSequence[j+i]]-Pmean)
        if len(ProteinSequence)-i==0:
            Result['MoranAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence))/K,3)
        else:
            Result['MoranAuto'+AAPName+str(i)]=round(temp/(len(ProteinSequence)-i)/K,3)

    return Result


# In[9]:


def CalculateEachGearyAuto(ProteinSequence,AAP,AAPName):

    AAPdic=NormalizeEachAAP(AAP)

    cc=[]
    for i in ProteinSequence:
        cc.append(AAPdic[i])

    K=((_std(cc))**2)*len(ProteinSequence)/(len(ProteinSequence)-1)
    Result={}
    for i in range(1,31):
        temp=0
        for j in range(len(ProteinSequence)-i):
                
            temp=temp+(AAPdic[ProteinSequence[j]]-AAPdic[ProteinSequence[j+i]])**2
        if len(ProteinSequence)-i==0:
            Result['GearyAuto'+AAPName+str(i)]=round(temp/(2*(len(ProteinSequence)))/K,3)
        else:
            Result['GearyAuto'+AAPName+str(i)]=round(temp/(2*(len(ProteinSequence)-i))/K,3)
    return Result


# In[10]:


def CalculateNormalizedMoreauBrotoAuto(ProteinSequence,AAProperty,AAPropertyName):
    Result={}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]]=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])


    return Result


# In[11]:


def CalculateMoranAuto(ProteinSequence,AAProperty,AAPropertyName):
    Result={}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]]=CalculateEachMoranAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])

    return Result


# In[12]:


def CalculateGearyAuto(ProteinSequence,AAProperty,AAPropertyName):
    Result={}
    for i in range(len(AAProperty)):
        Result[AAPropertyName[i]]=CalculateEachGearyAuto(ProteinSequence,AAProperty[i],AAPropertyName[i])

    return Result


# In[13]:


def CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[14]:


def CalculateNormalizedMoreauBrotoAutoAvFlexibility(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_AvFlexibility,'_AvFlexibility')
    return result


# In[15]:


def CalculateNormalizedMoreauBrotoAutoPolarizability(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[16]:


def CalculateNormalizedMoreauBrotoAutoFreeEnergy(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_FreeEnergy,'_FreeEnergy')
    return result


# In[17]:


def CalculateNormalizedMoreauBrotoAutoResidueASA(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
    return result


# In[18]:


def CalculateNormalizedMoreauBrotoAutoResidueVol(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_ResidueVol,'_ResidueVol')
    return result


# In[19]:


def CalculateNormalizedMoreauBrotoAutoSteric(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Steric,'_Steric')
    return result


# In[20]:


def CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence):
    result=CalculateEachNormalizedMoreauBrotoAuto(ProteinSequence,_Mutability,'_Mutability')
    return result


# In[21]:


def CalculateMoranAutoHydrophobicity(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[22]:


def CalculateMoranAutoAvFlexibility(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_AvFlexibility,'_AvFlexibility')
    return result


# In[23]:


def CalculateMoranAutoPolarizability(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[24]:


def CalculateMoranAutoFreeEnergy(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_FreeEnergy,'_FreeEnergy')
    return result


# In[25]:


def CalculateMoranAutoResidueASA(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
    return result


# In[26]:


def CalculateMoranAutoResidueVol(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_ResidueVol,'_ResidueVol')
    return result


# In[27]:


def CalculateMoranAutoSteric(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_Steric,'_Steric')
    return result


# In[28]:


def CalculateMoranAutoMutability(ProteinSequence):
    result=CalculateEachMoranAuto(ProteinSequence,_Mutability,'_Mutability')
    return result


# In[29]:


def CalculateGearyAutoHydrophobicity(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_Hydrophobicity,'_Hydrophobicity')
    return result


# In[30]:


def CalculateGearyAutoAvFlexibility(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_AvFlexibility,'_AvFlexibility')
    return result


# In[31]:


def CalculateGearyAutoPolarizability(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_Polarizability,'_Polarizability')
    return result


# In[32]:


def CalculateGearyAutoFreeEnergy(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_FreeEnergy,'_FreeEnergy')
    return result


# In[33]:


def CalculateGearyAutoResidueASA(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_ResidueASA,'_ResidueASA')
    return result


# In[34]:


def CalculateGearyAutoResidueVol(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_ResidueVol,'_ResidueVol')
    return result


# In[35]:


def CalculateGearyAutoSteric(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_Steric,'_Steric')
    return result


# In[36]:


def CalculateGearyAutoMutability(ProteinSequence):
    result=CalculateEachGearyAuto(ProteinSequence,_Mutability,'_Mutability')
    return result


# In[37]:


def CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence):
    result={}
    result.update(CalculateNormalizedMoreauBrotoAutoHydrophobicity(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoAvFlexibility(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoPolarizability(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoFreeEnergy(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoResidueASA(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoResidueVol(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoSteric(ProteinSequence))
    result.update(CalculateNormalizedMoreauBrotoAutoMutability(ProteinSequence))
    return result


# In[38]:


def CalculateMoranAutoTotal(ProteinSequence):
    result={}
    result.update(CalculateMoranAutoHydrophobicity(ProteinSequence))
    result.update(CalculateMoranAutoAvFlexibility(ProteinSequence))
    result.update(CalculateMoranAutoPolarizability(ProteinSequence))
    result.update(CalculateMoranAutoFreeEnergy(ProteinSequence))
    result.update(CalculateMoranAutoResidueASA(ProteinSequence))
    result.update(CalculateMoranAutoResidueVol(ProteinSequence))
    result.update(CalculateMoranAutoSteric(ProteinSequence))
    result.update(CalculateMoranAutoMutability(ProteinSequence))
    return result


# In[39]:


def CalculateGearyAutoTotal(ProteinSequence):
    result={}
    result.update(CalculateGearyAutoHydrophobicity(ProteinSequence))
    result.update(CalculateGearyAutoAvFlexibility(ProteinSequence))
    result.update(CalculateGearyAutoPolarizability(ProteinSequence))
    result.update(CalculateGearyAutoFreeEnergy(ProteinSequence))
    result.update(CalculateGearyAutoResidueASA(ProteinSequence))
    result.update(CalculateGearyAutoResidueVol(ProteinSequence))
    result.update(CalculateGearyAutoSteric(ProteinSequence))
    result.update(CalculateGearyAutoMutability(ProteinSequence))
    return result


# In[40]:


def CalculateAutoTotal(ProteinSequence):
    result={}
    result.update(CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence))
    result.update(CalculateMoranAutoTotal(ProteinSequence))
    result.update(CalculateGearyAutoTotal(ProteinSequence))
    return result


# In[41]:


def CalculateAutocorrelation(ProteinSequence):
    AALetter=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    _Hydrophobicity={"A":0.02,"R":-0.42,"N":-0.77,"D":-1.04,"C":0.77,"Q":-1.10,"E":-1.14,"G":-0.80,"H":0.26,"I":1.81,"L":1.14,"K":-0.41,"M":1.00,"F":1.35,"P":-0.09,"S":-0.97,"T":-0.77,"W":1.71,"Y":1.11,"V":1.13}

    _AvFlexibility={"A":0.357,"R":0.529,"N":0.463,"D":0.511,"C":0.346,"Q":0.493,"E":0.497,"G":0.544,"H":0.323,"I":0.462,"L":0.365,"K":0.466,"M":0.295,"F":0.314,"P":0.509,"S":0.507,"T":0.444,"W":0.305,"Y":0.420,"V":0.386}

    _Polarizability={"A":0.046,"R":0.291,"N":0.134,"D":0.105,"C":0.128,"Q":0.180,"E":0.151,"G":0.000,"H":0.230,"I":0.186,"L":0.186,"K":0.219,"M":0.221,"F":0.290,"P":0.131,"S":0.062,"T":0.108,"W":0.409,"Y":0.298,"V":0.140}

    _FreeEnergy={"A":-0.368,"R":-1.03,"N":0.0,"D":2.06,"C":4.53,"Q":0.731,"E":1.77,"G":-0.525,"H":0.0,"I":0.791,"L":1.07,"K":0.0,"M":0.656,"F":1.06,"P":-2.24,"S":-0.524,"T":0.0,"W":1.60,"Y":4.91,"V":0.401}

    _ResidueASA={"A":115.0,"R":225.0,"N":160.0,"D":150.0,"C":135.0,"Q":180.0,"E":190.0,"G":75.0,"H":195.0,"I":175.0,"L":170.0,"K":200.0,"M":185.0,"F":210.0,"P":145.0,"S":115.0,"T":140.0,"W":255.0,"Y":230.0,"V":155.0}

    _ResidueVol={"A":52.6,"R":109.1,"N":75.7,"D":68.4,"C":68.3,"Q":89.7,"E":84.7,"G":36.3,"H":91.9,"I":102.0,"L":102.0,"K":105.1,"M":97.7,"F":113.9,"P":73.6,"S":54.9,"T":71.2,"W":135.4,"Y":116.2,"V":85.1}

    _Steric={"A":0.52,"R":0.68,"N":0.76,"D":0.76,"C":0.62,"Q":0.68,"E":0.68,"G":0.00,"H":0.70,"I":1.02,"L":0.98,"K":0.68,"M":0.78,"F":0.70,"P":0.36,"S":0.53,"T":0.50,"W":0.70,"Y":0.70,"V":0.76}

    _Mutability={"A":100.0,"R":65.0,"N":134.0,"D":106.0,"C":20.0,"Q":93.0,"E":102.0,"G":49.0,"H":66.0,"I":96.0,"L":40.0,"K":-56.0,"M":94.0,"F":41.0,"P":56.0,"S":120.0,"T":97.0,"W":18.0,"Y":41.0,"V":74.0}
    result = {}
    result.update(CalculateNormalizedMoreauBrotoAutoTotal(ProteinSequence))
    result.update(CalculateMoranAutoTotal(ProteinSequence))
    result.update(CalculateGearyAutoTotal(ProteinSequence))
    result.update(CalculateAutoTotal(ProteinSequence))
    return result


# In[ ]:




