#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf as MS

import math
import numpy


# In[2]:


def CalculateMolLogP(mol):
    return round(Crippen.MolLogP(mol),3)


# In[3]:


def CalculateMolLogP2(mol):
    res=Crippen.MolLogP(mol)
    
    return round(res**2,3)


# In[4]:


def CalculateMolMR(mol):
    return round(Crippen.MolMR(mol),3)


# In[5]:


def CalculateTPSA(mol):
    return round(MS.TPSA(mol),3)


# In[6]:


def _CalculateBondNumber(mol,bondtype='SINGLE'):
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name==bondtype:
            i=i+1
            
    return i


# In[7]:


def CalculateUnsaturationIndex(mol):
    nd=_CalculateBondNumber(mol,bondtype='DOUBLE')
    nt=_CalculateBondNumber(mol,bondtype='TRIPLE')
    na=_CalculateBondNumber(mol,bondtype='AROMATIC')
    res=math.log((1+nd+nt+na),2)
    
    return round(res,3)


# In[8]:


def CalculateHydrophilicityFactor(mol):
    nheavy=mol.GetNumHeavyAtoms()
    nc=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==6:
            nc=nc+1
    nhy=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==7 or atom.GetAtomicNum()==8 or atom.GetAtomicNum()==16:
            atomn=atom.GetNeighbors()
            for i in atomn:
                if i.GetAtomicNum()==1:
                    nhy=nhy+1
                
    res=(1+nhy)*math.log((1+nhy),2)+nc*(1.0/nheavy*math.log(1.0/nheavy,2))+math.sqrt((nhy+0.0)/(nheavy^2))
    return round(res,3)


# In[9]:


def GetAllMoleculePropertyFeatures_ligand(mol):
    LogP = CalculateMolLogP(mol)
    LogP2 = CalculateMolLogP2(mol)
    MR = CalculateMolMR(mol)
    TPSA = CalculateTPSA(mol)
    Hy = CalculateHydrophilicityFactor(mol)
    UI = CalculateUnsaturationIndex(mol)
    
    return[
        LogP, 
        LogP2, 
        MR, 
        TPSA, 
        Hy, 
        UI
    ]


# In[10]:


def GetAllMoleculePropertyFeatures_pocket(mol):
    LogP = CalculateMolLogP(mol)
    LogP2 = CalculateMolLogP2(mol)
    MR = CalculateMolMR(mol)
    TPSA = CalculateTPSA(mol)
    Hy = CalculateHydrophilicityFactor(mol)
    UI = CalculateUnsaturationIndex(mol)
    
    return[
        LogP, 
        LogP2, 
        MR, 
        TPSA, 
        Hy, 
        UI
    ]


# In[11]:


def GetMolecularProperty(mol):
    MolecularProperty={'LogP':CalculateMolLogP,
                   'LogP2':CalculateMolLogP2,
                   'MR':CalculateMolMR,
                   'TPSA':CalculateTPSA,
                   'Hy':CalculateHydrophilicityFactor,
                   'UI':CalculateUnsaturationIndex
    }
    result={}
    for DesLabel in MolecularProperty.keys():
        try:
            result[DesLabel]=MolecularProperty[DesLabel](mol)
        except:
            result[DesLabel]=numpy.nan
    return result


# In[ ]:




