#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem import MolSurf as MOE 
from rdkit.Chem.EState import EState_VSA as EVSA


# In[2]:


def CalculateLabuteASA(mol):
    res={}
    temp=MOE.pyLabuteASA(mol,includeHs=1)
    res['LabuteASA']=round(temp,3)
    return res


# In[3]:


def CalculateTPSA(mol):
    res={}
    temp=MOE.TPSA(mol)
    res['TPSA']=round(temp,3)
    return res


# In[4]:


def CalculateSLOGPVSA(mol,bins=None):
    temp=MOE.SlogP_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['slogPVSA'+str(i)]=round(j,3)
    return res


# In[5]:


def CalculateSMRVSA(mol,bins=None):
    temp=MOE.SMR_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['MRVSA'+str(i)]=round(j,3)
    return res


# In[6]:


def CalculatePEOEVSA(mol,bins=None):
    temp=MOE.PEOE_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['PEOEVSA'+str(i)]=round(j,3)
    return res    


# In[7]:


def CalculateEstateVSA(mol,bins=None):
    temp=EVSA.EState_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['EstateVSA'+str(i)]=round(j,3)
    return res


# In[8]:


def CalculateVSAEstate(mol,bins=None):
    temp=EVSA.VSA_EState_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['VSAEstate'+str(i)]=round(j,3)
    return res


# In[9]:


def GetMOE(mol):
    result={}
    result.update(CalculateLabuteASA(mol))
    result.update(CalculateTPSA(mol))
    result.update(CalculateSLOGPVSA(mol,bins=None))
    result.update(CalculateSMRVSA(mol,bins=None))
    result.update(CalculatePEOEVSA(mol,bins=None))
    result.update(CalculateEstateVSA(mol,bins=None))
    result.update(CalculateVSAEstate(mol,bins=None))
    return result


# In[ ]:




