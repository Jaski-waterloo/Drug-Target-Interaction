#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem import rdPartialCharges as GMCharge
import numpy


# In[2]:


iter_step = 12


# In[3]:


def _CalculateElementMaxPCharge(mol,AtomicNum=6):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        return round(max(res),3)


# In[4]:


def _CalculateElementMaxNCharge(mol,AtomicNum=6):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(min(res),3)


# In[5]:


def CalculateHMaxPCharge(mol):
     return _CalculateElementMaxPCharge(mol,AtomicNum=1)


# In[6]:


def CalculateCMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol,AtomicNum=6)


# In[7]:


def CalculateNMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol,AtomicNum=7)


# In[8]:


def CalculateOMaxPCharge(mol):
    return _CalculateElementMaxPCharge(mol,AtomicNum=8)


# In[9]:


def CalculateHMaxNCharge(mol):
    return _CalculateElementMaxNCharge(mol,AtomicNum=1)


# In[10]:


def CalculateCMaxNCharge(mol):
    return _CalculateElementMaxNCharge(mol,AtomicNum=6)


# In[11]:


def CalculateNMaxNCharge(mol):
    return _CalculateElementMaxNCharge(mol,AtomicNum=7)


# In[12]:


def CalculateOMaxNCharge(mol):
    return _CalculateElementMaxNCharge(mol,AtomicNum=8)


# In[13]:


def CalculateAllMaxPCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(max(res),3)


# In[14]:


def CalculateAllMaxNCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(min(res),3)


# In[15]:


def CalculateTotalPCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(sum(cc[cc>0]),3)


# In[16]:


def CalculateMeanPCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(numpy.mean(cc[cc>0]),3)


# In[17]:


def CalculateTotalNCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(sum(cc[cc<0]),3)


# In[18]:


def CalculateMeanNCharge(mol):
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(numpy.mean(cc[cc<0]),3)


# In[19]:


def CalculateLocalDipoleIndex(mol):
    GMCharge.ComputeGasteigerCharges(mol,iter_step)
    res=[]
    for atom in mol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    cc = [numpy.absolute(res[x.GetBeginAtom().GetIdx()]-res[x.GetEndAtom().GetIdx()]) for x in mol.GetBonds()]
    B=len(mol.GetBonds())
    
    return round(sum(cc)/B,3)


# In[20]:


def CalculateSubmolPolarityPara(mol):
    return round(CalculateAllMaxPCharge(mol)-CalculateAllMaxNCharge(mol),3)


# In[21]:


def GetAllChargeProperties_pocket(mol):
        SPP = CalculateSubmolPolarityPara(mol)
        LDI = CalculateLocalDipoleIndex(mol)
        Mac = CalculateMeanAbsoulteCharge(mol)
        Tac = CalculateTotalAbsoulteCharge(mol)
        Mnc = CalculateMeanNCharge(mol)
        Tnc = CalculateTotalNCharge(mol)
        Mpc = CalculateMeanPCharge(mol)
        Tpc = CalculateTotalPCharge(mol)
        Qmin = CalculateAllMaxNCharge(mol)
        Qmax = CalculateAllMaxPCharge(mol)
        QOmin = CalculateOMaxNCharge(mol)
        QNmin = CalculateNMaxNCharge(mol)
        QCmin = CalculateCMaxNCharge(mol)
        QHmin = CalculateHMaxNCharge(mol)
        QOmax = CalculateOMaxPCharge(mol)
        QNmax = CalculateNMaxPCharge(mol)
        QCmax = CalculateCMaxPCharge(mol)
        QHmax = CalculateHMaxPCharge(mol)
        
        return [SPP, LDI,
                Mac, Tac, Mnc, Tnc,
                Mpc, Tpc, Qmin,
                Qmax, QOmin, QNmin, QCmin,
                QHmin, QOmax, QNmax, QCmax, QHmax]


# In[22]:


def GetAllChargeProperties_ligand(mol):
        SPP = CalculateSubmolPolarityPara(mol)
        LDI = CalculateLocalDipoleIndex(mol)
        #Mac = CalculateMeanAbsoulteCharge(mol)
        #Tac = CalculateTotalAbsoulteCharge(mol)
        Mnc = CalculateMeanNCharge(mol)
        Tnc = CalculateTotalNCharge(mol)
        Mpc = CalculateMeanPCharge(mol)
        Tpc = CalculateTotalPCharge(mol)
        Qmin = CalculateAllMaxNCharge(mol)
        Qmax = CalculateAllMaxPCharge(mol)
        QOmin = CalculateOMaxNCharge(mol)
        QNmin = CalculateNMaxNCharge(mol)
        QCmin = CalculateCMaxNCharge(mol)
        QHmin = CalculateHMaxNCharge(mol)
        QOmax = CalculateOMaxPCharge(mol)
        QNmax = CalculateNMaxPCharge(mol)
        QCmax = CalculateCMaxPCharge(mol)
        QHmax = CalculateHMaxPCharge(mol)
        
        return [SPP, LDI,
                Mnc, Tnc,
                Mpc, Tpc, Qmin,
                Qmax, QOmin, QNmin, QCmin,
                QHmin, QOmax, QNmax, QCmax, QHmax]


# In[23]:


def GetCharge(mol):
    _Charge={'SPP':CalculateSubmolPolarityPara,
        'LDI':CalculateLocalDipoleIndex,
        'Mnc':CalculateMeanNCharge,
        'Tnc':CalculateTotalNCharge,
        'Mpc':CalculateMeanPCharge,
        'Tpc':CalculateTotalPCharge,
        'Qmin':CalculateAllMaxNCharge,
        'Qmax':CalculateAllMaxPCharge,
        'QOmin':CalculateOMaxNCharge,
        'QNmin':CalculateNMaxNCharge,
        'QCmin':CalculateCMaxNCharge,
        'QHmin':CalculateHMaxNCharge,
        'QOmax':CalculateOMaxPCharge,
        'QNmax':CalculateNMaxPCharge,
        'QCmax':CalculateCMaxPCharge,
        'QHmax':CalculateHMaxPCharge,
    }
    result={}
    for DesLabel in _Charge.keys():
        try:
            result[DesLabel]=_Charge[DesLabel](mol)
        except:
            result[DesLabel]=numpy.nan
    return result


# In[ ]:




