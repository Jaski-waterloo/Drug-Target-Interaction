#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
from rdkit.Chem import rdchem
import numpy


# In[2]:


periodicTable = rdchem.GetPeriodicTable()


# In[3]:


def CalculateChi0(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum(numpy.sqrt(1./deltas))
    return res


# In[4]:


def CalculateChi1(mol):
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = sum(numpy.sqrt(1./cc))
    return res


# In[5]:


def CalculateMeanRandic(mol):
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = numpy.mean(numpy.sqrt(1./cc))
    
    return res


# In[6]:


def _CalculateChinp(mol,NumPath=2):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    for path in Chem.FindAllPathsOfLengthN(mol,NumPath+1,useBonds=0):
        cAccum=1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1./numpy.sqrt(cAccum)
    return accum


# In[7]:


def CalculateChi2(mol):
    return _CalculateChinp(mol,NumPath=2)


# In[8]:


def CalculateChi3p(mol):
    return _CalculateChinp(mol,NumPath=3)


# In[9]:


def CalculateChi4p(mol):
    return _CalculateChinp(mol,NumPath=4)


# In[10]:


def CalculateChi5p(mol):
    return _CalculateChinp(mol,NumPath=5)


# In[11]:


def CalculateChi6p(mol):
    return _CalculateChinp(mol,NumPath=6)


# In[12]:


def CalculateChi7p(mol):
    return _CalculateChinp(mol,NumPath=7)


# In[13]:


def CalculateChi8p(mol):
    return _CalculateChinp(mol,NumPath=8)


# In[14]:


def CalculateChi9p(mol):
    return _CalculateChinp(mol,NumPath=9)


# In[15]:


def CalculateChi10p(mol):
    return _CalculateChinp(mol,NumPath=10)


# In[16]:


def CalculateChi3c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[17]:


def CalculateChi4c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[18]:


def CalculateChi4pc(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[mol.GetAtomWithIdx(x).GetDegree() for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[19]:


def CalculateDeltaChi3c4pc(mol):
    return abs(CalculateChi3c(mol)-CalculateChi4pc(mol))


# In[20]:


def _CalculateChinch(mol,NumCycle=3):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    for tup in mol.GetRingInfo().AtomRings():
        cAccum=1.0
        if len(tup)==NumCycle:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1./numpy.sqrt(cAccum)

    return accum    


# In[21]:


def CalculateChi3ch(mol):
    return _CalculateChinch(mol,NumCycle=3)


# In[22]:


def CalculateChi4ch(mol):
    return _CalculateChinch(mol,NumCycle=4)


# In[23]:


def CalculateChi5ch(mol):
    return _CalculateChinch(mol,NumCycle=5)


# In[24]:


def CalculateChi6ch(mol):
    return _CalculateChinch(mol,NumCycle=6)


# In[25]:


def _HKDeltas(mol,skipHs=1):
    global periodicTable
    res=[]
    for atom in mol.GetAtoms():
        n=atom.GetAtomicNum()
        if n>1:
            nV=periodicTable.GetNOuterElecs(n)
            nHs=atom.GetTotalNumHs()
            if n<10:
                res.append(float(nV-nHs))
            else:
                res.append(float(nV-nHs)/float(n-nV-1))
        elif not skipHs:
            res.append(0.0)
    return res


# In[26]:


def CalculateChiv0(mol):
    deltas=_HKDeltas(mol,skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum(numpy.sqrt(1./deltas))
    return res


# In[27]:


def _CalculateChivnp(mol,NumPath=1):
    accum=0.0
    deltas=_HKDeltas(mol,skipHs=0)
    for path in Chem.FindAllPathsOfLengthN(mol,NumPath+1,useBonds=0):
        cAccum=1.0
        for idx in path:
            cAccum *= deltas[idx]
        if cAccum:
            accum += 1./numpy.sqrt(cAccum)
    return accum


# In[28]:


def CalculateChiv1(mol):
    return _CalculateChivnp(mol,NumPath=1)


# In[29]:


def CalculateChiv2(mol):
    return _CalculateChivnp(mol,NumPath=2)


# In[30]:


def CalculateChiv3p(mol):
    return _CalculateChivnp(mol,NumPath=3)


# In[31]:


def CalculateChiv4p(mol):
    return _CalculateChivnp(mol,NumPath=4)


# In[32]:


def CalculateChiv5p(mol):
    return _CalculateChivnp(mol,NumPath=5)


# In[33]:


def CalculateChiv6p(mol):
     return _CalculateChivnp(mol,NumPath=6)


# In[34]:


def CalculateChiv7p(mol):
    return _CalculateChivnp(mol,NumPath=7)


# In[35]:


def CalculateChiv8p(mol):
    return _CalculateChivnp(mol,NumPath=8)


# In[36]:


def CalculateChiv9p(mol):
    return _CalculateChivnp(mol,NumPath=9)


# In[37]:


def CalculateChiv10p(mol):
    return _CalculateChivnp(mol,NumPath=10)


# In[38]:


def CalculateDeltaChi0(mol):
    return abs(CalculateChiv0(mol)-CalculateChi0(mol))


# In[39]:


def CalculateDeltaChi1(mol):
    return abs(CalculateChiv1(mol)-CalculateChi1(mol))


# In[40]:


def CalculateDeltaChi2(mol):
    return abs(_CalculateChivnp(mol,NumPath=2)-_CalculateChinp(mol,NumPath=2))


# In[41]:


def CalculateDeltaChi3(mol):
    return abs(_CalculateChivnp(mol,NumPath=3)-_CalculateChinp(mol,NumPath=3))


# In[42]:


def CalculateDeltaChi4(mol):
    return abs(_CalculateChivnp(mol,NumPath=4)-_CalculateChinp(mol,NumPath=4))


# In[43]:


def _AtomHKDeltas(atom,skipHs=0):
    global periodicTable
    res=[]
    n=atom.GetAtomicNum()
    if n>1:
        nV=periodicTable.GetNOuterElecs(n)
        nHs=atom.GetTotalNumHs()
        if n<10:
            res.append(float(nV-nHs))
        else:
            res.append(float(nV-nHs)/float(n-nV-1))
    elif not skipHs:
        res.append(0.0)
    return res


# In[44]:


def CalculateChiv3c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[45]:


def CalculateChiv4c(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)(~*)~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[46]:


def CalculateChiv4pc(mol):
    accum=0.0
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    patt=Chem.MolFromSmarts('*~*(~*)~*~*')
    HPatt=mol.GetSubstructMatches(patt)
    #print HPatt
    for cluster in HPatt:
        deltas=[_AtomHKDeltas(mol.GetAtomWithIdx(x)) for x in cluster]
        while 0 in deltas:
            deltas.remove(0)
        if deltas!=[]:
            deltas1=numpy.array(deltas,numpy.float)
            accum=accum+1./numpy.sqrt(deltas1.prod())
    return accum


# In[47]:


def CalculateDeltaChiv3c4pc(mol):
    return abs(CalculateChiv3c(mol)-CalculateChiv4pc(mol))


# In[48]:


def _CalculateChivnch(mol,NumCyc=3):
    accum=0.0
    deltas=_HKDeltas(mol,skipHs=0)
    for tup in mol.GetRingInfo().AtomRings():
        cAccum=1.0
        if len(tup)==NumCyc:
            for idx in tup:
                cAccum *= deltas[idx]
            if cAccum:
                accum += 1./numpy.sqrt(cAccum)

    return accum


# In[49]:


def CalculateChiv3ch(mol):
    return _CalculateChivnch(mol,NumCyc=3)


# In[50]:


def CalculateChiv4ch(mol):
    return _CalculateChivnch(mol,NumCyc=4)


# In[51]:


def CalculateChiv5ch(mol):
    return _CalculateChivnch(mol,NumCyc=5)


# In[52]:


def CalculateChiv6ch(mol):
    return _CalculateChivnch(mol,NumCyc=6)


# In[53]:


def GetAllConnectivityProperties_pocket(mol):
    Chi0 = CalculateChi0(mol)
    Chi1 = CalculateChi1(mol)
    mChi1 = CalculateMeanRandic(mol)
    mChi1 = CalculateMeanRandic(mol)
    Chi2 = CalculateChi2(mol)
    Chi3 = CalculateChi3p(mol)
    Chi4 = CalculateChi4p(mol)
    Chi5 = CalculateChi5p(mol)
    Chi6 = CalculateChi6p(mol)
    Chi7 = CalculateChi7p(mol)
    Chi8 = CalculateChi8p(mol)
    Chi9 = CalculateChi9p(mol)
    Chi10 = CalculateChi10p(mol)
    Chi3c = CalculateChi3c(mol)
    Chi4c = CalculateChi4c(mol)
    Chi4pc = CalculateChi4pc(mol)
    Chi3ch = CalculateChi3ch(mol)
    Chi4ch = CalculateChi4ch(mol)
    Chi5ch = CalculateChi5ch(mol)
    Chi6ch = CalculateChi6ch(mol)
    knotp = CalculateDeltaChi3c4pc(mol)
    Chiv0 = CalculateChiv0(mol)
    Chiv1 = CalculateChiv1(mol)
    Chiv2 = CalculateChiv2(mol)
    Chiv3 = CalculateChiv3p(mol)
    Chiv4 = CalculateChiv4p(mol)
    Chiv5 = CalculateChiv5p(mol)
    Chiv6 = CalculateChiv6p(mol)
    Chiv7 = CalculateChiv7p(mol)
    Chiv8 = CalculateChiv8p(mol)
    Chiv9 = CalculateChiv9p(mol)
    Chiv10 = CalculateChiv10p(mol)
    dchi0 = CalculateDeltaChi0(mol)
    dchi1 = CalculateDeltaChi1(mol)
    dchi2 = CalculateDeltaChi2(mol)
    dchi3 = CalculateDeltaChi3(mol)
    dchi4 = CalculateDeltaChi4(mol)
    Chiv3c = CalculateChiv3c(mol)
    Chiv4c = CalculateChiv4c(mol)
    Chiv4pc = CalculateChiv4pc(mol)
    Chiv3ch = CalculateChiv3ch(mol)
    Chiv4ch = CalculateChiv4ch(mol)
    Chiv5ch = CalculateChiv5ch(mol)
    Chiv6ch = CalculateChiv6ch(mol)
    knotpv = CalculateDeltaChiv3c4pc(mol)
    
    return[Chi0,
    Chi1,
    mChi1,
    mChi1,
    Chi2,
    Chi3,
    Chi4,
    Chi5,
    Chi6,
    Chi7,
    Chi8,
    Chi9,
    Chi10,
    Chi3c,
    Chi4c,
    Chi4pc,
    Chi3ch,
    Chi4ch,
    Chi5ch,
    Chi6ch,
    knotp,
    Chiv0,
    Chiv1,
    Chiv2,
    Chiv3,
    Chiv4,
    Chiv5,
    Chiv6,
    Chiv7,
    Chiv8,
    Chiv9,
    Chiv10,
    dchi0,
    dchi1,
    dchi2,
    dchi3,
    dchi4,
    Chiv3c,
    Chiv4c,
    Chiv4pc,
    Chiv3ch,
    Chiv4ch,
    Chiv5ch,
    Chiv6ch,
    knotp]


# In[54]:


def GetAllConnectivityProperties_ligand(mol):
    Chi0 = CalculateChi0(mol)
    Chi1 = CalculateChi1(mol)
    mChi1 = CalculateMeanRandic(mol)
    Chi2 = CalculateChi2(mol)
    Chi3 = CalculateChi3p(mol)
    Chi4 = CalculateChi4p(mol)
    Chi5 = CalculateChi5p(mol)
    Chi6 = CalculateChi6p(mol)
    Chi7 = CalculateChi7p(mol)
    Chi8 = CalculateChi8p(mol)
    Chi9 = CalculateChi9p(mol)
    Chi10 = CalculateChi10p(mol)
    Chi3c = CalculateChi3c(mol)
    Chi4c = CalculateChi4c(mol)
    Chi4pc = CalculateChi4pc(mol)
    Chi3ch = CalculateChi3ch(mol)
    Chi4ch = CalculateChi4ch(mol)
    Chi5ch = CalculateChi5ch(mol)
    Chi6ch = CalculateChi6ch(mol)
    knotp = CalculateDeltaChi3c4pc(mol)
    Chiv0 = CalculateChiv0(mol)
    Chiv1 = CalculateChiv1(mol)
    Chiv2 = CalculateChiv2(mol)
    Chiv3 = CalculateChiv3p(mol)
    Chiv4 = CalculateChiv4p(mol)
    Chiv5 = CalculateChiv5p(mol)
    Chiv6 = CalculateChiv6p(mol)
    Chiv7 = CalculateChiv7p(mol)
    Chiv8 = CalculateChiv8p(mol)
    Chiv9 = CalculateChiv9p(mol)
    Chiv10 = CalculateChiv10p(mol)
    dchi0 = CalculateDeltaChi0(mol)
    dchi1 = CalculateDeltaChi1(mol)
    dchi2 = CalculateDeltaChi2(mol)
    dchi3 = CalculateDeltaChi3(mol)
    dchi4 = CalculateDeltaChi4(mol)
    Chiv3c = CalculateChiv3c(mol)
    Chiv4c = CalculateChiv4c(mol)
    Chiv4pc = CalculateChiv4pc(mol)
    Chiv3ch = CalculateChiv3ch(mol)
    Chiv4ch = CalculateChiv4ch(mol)
    Chiv5ch = CalculateChiv5ch(mol)
    Chiv6ch = CalculateChiv6ch(mol)
    knotpv = CalculateDeltaChiv3c4pc(mol)
    
    return[Chi0,
    Chi1,
    mChi1,
    Chi2,
    Chi3,
    Chi4,
    Chi5,
    Chi6,
    Chi7,
    Chi8,
    Chi9,
    Chi10,
    Chi3c,
    Chi4c,
    Chi4pc,
    Chi3ch,
    Chi4ch,
    Chi5ch,
    Chi6ch,
    knotp,
    Chiv0,
    Chiv1,
    Chiv2,
    Chiv3,
    Chiv4,
    Chiv5,
    Chiv6,
    Chiv7,
    Chiv8,
    Chiv9,
    Chiv10,
    dchi0,
    dchi1,
    dchi2,
    dchi3,
    dchi4,
    Chiv3c,
    Chiv4c,
    Chiv4pc,
    Chiv3ch,
    Chiv4ch,
    Chiv5ch,
    Chiv6ch,
    knotp]


# In[55]:


def GetConnectivity(mol):
    _connectivity={'Chi0':CalculateChi0,
                       'Chi1':CalculateChi1,
                       'mChi1':CalculateMeanRandic,
                       'Chi2':CalculateChi2,
                       'Chi3':CalculateChi3p,
                       'Chi4':CalculateChi4p,
                       'Chi5':CalculateChi5p,
                       'Chi6':CalculateChi6p,
                       'Chi7':CalculateChi7p,
                       'Chi8':CalculateChi8p,
                       'Chi9':CalculateChi9p,
                       'Chi10':CalculateChi10p,
                       'Chi3c':CalculateChi3c,
                       'Chi4c':CalculateChi4c,
                       'Chi4pc':CalculateChi4pc,
                       'Chi3ch':CalculateChi3ch,
                       'Chi4ch':CalculateChi4ch,
                       'Chi5ch':CalculateChi5ch,
                       'Chi6ch':CalculateChi6ch,
                       'knotp':CalculateDeltaChi3c4pc,
                       'Chiv0':CalculateChiv0,
                      'Chiv1':CalculateChiv1,
                      'Chiv2':CalculateChiv2,
                      'Chiv3':CalculateChiv3p,
                      'Chiv4':CalculateChiv4p,
                       'Chiv5':CalculateChiv5p,
                       'Chiv6':CalculateChiv6p,
                       'Chiv7':CalculateChiv7p,
                       'Chiv8':CalculateChiv8p,
                       'Chiv9':CalculateChiv9p,
                       'Chiv10':CalculateChiv10p,
                       'dchi0':CalculateDeltaChi0,
                       'dchi1':CalculateDeltaChi1,
                       'dchi2':CalculateDeltaChi2,
                       'dchi3':CalculateDeltaChi3,
                       'dchi4':CalculateDeltaChi4,
                       'Chiv3c':CalculateChiv3c,
                       'Chiv4c':CalculateChiv4c,
                       'Chiv4pc':CalculateChiv4pc,
                       'Chiv3ch':CalculateChiv3ch,
                       'Chiv4ch':CalculateChiv4ch,
                       'Chiv5ch':CalculateChiv5ch,
                       'Chiv6ch':CalculateChiv6ch,
                       'knotpv':CalculateDeltaChiv3c4pc
    }
    result={}
    for DesLabel in _connectivity.keys():
        try:
            result[DesLabel]=round(_connectivity[DesLabel](mol),3)
        except:
            result[DesLabel]=numpy.nan
    return result


# In[ ]:




