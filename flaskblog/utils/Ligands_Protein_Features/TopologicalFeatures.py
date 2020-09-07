#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy

from rdkit.Chem import GraphDescriptors as GD
from rdkit import Chem
from rdkit.Chem import rdchem


# In[2]:


def _GetPrincipleQuantumNumber(atNum):
    if atNum<=2:
        return 1
    elif atNum<=10:
        return 2
    elif atNum<=18:
        return 3
    elif atNum<=36:
        return 4
    elif atNum<=54:
        return 5
    elif atNum<=86:
        return 6
    else:
        return 7


# In[3]:


periodicTable = rdchem.GetPeriodicTable()


# In[4]:


def CalculateWeiner(mol):
    return 1.0/2*sum(sum(Chem.GetDistanceMatrix(mol)))


# In[5]:


def CalculateMeanWeiner(mol):
    N=mol.GetNumAtoms()
    WeinerNumber=CalculateWeiner(mol)
    return 2.0*WeinerNumber/(N*(N-1))


# In[6]:


def CalculateBalaban(mol):
    adjMat=Chem.GetAdjacencyMatrix(mol)
    Distance= Chem.GetDistanceMatrix(mol)
    Nbond=mol.GetNumBonds()
    Natom=mol.GetNumAtoms()
    S=numpy.sum(Distance,axis=1)
    mu=Nbond-Natom+1
    sumk=0.
    for i in range(len(Distance)):
        si=S[i]
        for j in range(i,len(Distance)):
            if adjMat[i,j]==1:
                sumk += 1./numpy.sqrt(si*S[j])
    if mu+1 !=0:
        J=float(Nbond)/float(mu+1)*sumk
    else:
        J=0
    return J


# In[7]:


def CalculateGraphDistance(mol):
    Distance= Chem.GetDistanceMatrix(mol)
    n=int(Distance.max())
    res=0.0
    for i in range(n):
       # print Distance==i+1
        temp=1./2*sum(sum(Distance==i+1))
        #print temp
        res = res+temp**2

    return numpy.log10(res)


# In[8]:


def CalculateDiameter(mol):
    Distance=Chem.GetDistanceMatrix(mol)

    return Distance.max()


# In[9]:


def CalculateRadius(mol):
    Distance=Chem.GetDistanceMatrix(mol)
    temp=[]
    for i in Distance:
        temp.append(max(i))
    return min(temp)


# In[10]:


def CalculatePetitjean(mol):
    diameter=CalculateDiameter(mol)
    radius=CalculateRadius(mol)
    return 1-radius/float(diameter)


# In[11]:


def CalculateXuIndex(mol):
    nAT=mol.GetNumAtoms()
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    Distance= Chem.GetDistanceMatrix(mol)
    sigma=numpy.sum(Distance,axis=1)
    temp1=0.0
    temp2=0.0
    for i in range(nAT):
        temp1=temp1+deltas[i]*((sigma[i])**2)
        temp2=temp2+deltas[i]*(sigma[i])
    Xu=numpy.sqrt(nAT)*numpy.log(temp1/temp2)
    
    return Xu


# In[12]:


def CalculateGutmanTopo(mol):
    nAT=mol.GetNumAtoms()
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    Distance= Chem.GetDistanceMatrix(mol)
    res=0.0
    for i in range(nAT):
        for j in range(i+1,nAT):
            res=res+deltas[i]*deltas[j]*Distance[i,j]

    return numpy.log10(res)


# In[13]:


def CalculatePolarityNumber(mol):
    Distance= Chem.GetDistanceMatrix(mol)
    res=1./2*sum(sum(Distance==3))
    
    return res


# In[14]:


def CalculatePoglianiIndex(mol):
    res=0.0
    for atom in mol.GetAtoms():
        n=atom.GetAtomicNum()
        nV=periodicTable.GetNOuterElecs(n)
        mP=_GetPrincipleQuantumNumber(n)
        res=res+(nV+0.0)/mP
    return res


# In[15]:


def CalculateIpc(mol):
    return numpy.log10(GD.Ipc(mol))


# In[16]:


def CalculateBertzCT(mol):
     return numpy.log10(GD.BertzCT(mol))


# In[17]:


def CalculateHarary(mol):
    Distance=numpy.array(Chem.GetDistanceMatrix(mol),'d')
                
    return 1.0/2*(sum(1.0/Distance[Distance!=0]))


# In[18]:


def CalculateSchiultz(mol):
    Distance=numpy.array(Chem.GetDistanceMatrix(mol),'d')
    Adjacent=numpy.array(Chem.GetAdjacencyMatrix(mol),'d')
    VertexDegree=sum(Adjacent)
    
    return sum(numpy.dot((Distance+Adjacent),VertexDegree))


# In[19]:


def CalculateZagreb1(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    return sum(numpy.array(deltas)**2)


# In[20]:


def CalculateZagreb2(mol):
    ke = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    
    return sum(ke)


# In[21]:


def CalculateMZagreb1(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    res=sum((1./deltas)**2)
    return res


# In[22]:


def CalculateMZagreb2(mol):
    cc = [x.GetBeginAtom().GetDegree()*x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc,'d')
    res = sum((1./cc)**2)
    return res


# In[23]:


def CalculateQuadratic(mol):
    M=CalculateZagreb1(mol)
    N=mol.GetNumAtoms()
    return 3-2*N+M/2.0


# In[24]:


def CalculatePlatt(mol):
    cc = [x.GetBeginAtom().GetDegree()+x.GetEndAtom().GetDegree()-2 for x in mol.GetBonds()]
    return sum(cc)


# In[25]:


def CalculateSimpleTopoIndex(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    res=numpy.prod(deltas)
    
    return numpy.log(res)


# In[26]:


def CalculateHarmonicTopoIndex(mol):
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')  
    nAtoms=mol.GetNumAtoms()
    
    res=nAtoms/sum(1./deltas)
    
    return res


# In[27]:


def CalculateGeometricTopoIndex(mol):
    nAtoms=mol.GetNumAtoms()
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    temp=numpy.prod(deltas)
    res=numpy.power(temp,1./nAtoms)

    return res    


# In[28]:


def CalculateArithmeticTopoIndex(mol):
    nAtoms=mol.GetNumAtoms()
    nBonds=mol.GetNumBonds()
    
    res=2.*nBonds/nAtoms
    return res


# In[29]:


def GetAllTopologicalProperties_pocket(mol):
    try:
        W = CalculateWeiner(mol)
    except:
        print("Could Not Calculate CalculateWeiner")
    try:
        AW = CalculateMeanWeiner(mol)
    except:
        print("Could Not Calculate CalculateMeanWeiner")
    try:
        J = CalculateBalaban(mol)
    except:
        print("Could Not Calculate CalculateBalaban")
    #try:
    #    Tigdi = CalculateGraphDistance(mol)
    #except:
    #    print("Could Not Calculate CalculateGraphDistance")
    try:
        Xu = CalculateXuIndex(mol)
    except:
        print("Could Not Calculate CalculateXuIndex")
    try:
        GMTI = CalculateGutmanTopo(mol)
    except:
        print("Could Not Calculate CalculateGutmanTopo")
    try:
        Pol = CalculatePolarityNumber(mol)
    except:
        print("Could Not Calculate CalculatePolarityNumber")
    try:
        DZ = CalculatePoglianiIndex(mol)
    except:
        print("Could Not Calculate CalculatePoglianiIndex")
    try:
        Ipc = CalculateIpc(mol)
    except:
        print("Could Not Calculate CalculateIpc")
    try:
        BertzCT = CalculateBertzCT(mol)
    except:
        print("Could Not Calculate CalculateBertzCT")
    try:
        Thara = CalculateHarary(mol)
    except:
        print("Could Not Calculate CalculateHarary")
    try:
        Tsch = CalculateSchiultz(mol)
    except:
        print("Could Not Calculate CalculateSchiultz")
    try:
        ZM1 = CalculateZagreb1(mol)
    except:
        print("Could Not Calculate CalculateZagreb1")
    try:
        ZM2 = CalculateZagreb2(mol)
    except:
        print("Could Not Calculate CalculateZagreb2")
    try:
        MZM1 = CalculateMZagreb1(mol)
    except:
        print("Could Not Calculate CalculateMZagreb1")
    try:
        MZM2 = CalculateMZagreb2(mol)
    except:
        print("Could Not Calculate CalculateMZagreb2")
    try:
        Qindex = CalculateQuadratic(mol)
    except:
        print("Could Not Calculate CalculateQuadratic")
    try:
        Platt = CalculatePlatt(mol)
    except:
        print("Could Not Calculate CalculatePlatt")
    #try:
    #    diametert = CalculateDiameter(mol)
    #except:
    #    print("Could Not Calculate CalculateDiameter")
    #try:
    #    radiust = CalculateRadius(mol)
    #except:
    #    print("Could Not Calculate CalculateRadius")
    try:
        petitjeant = CalculatePetitjean(mol)
    except:
        print("Could Not Calculate CalculatePetitjean")
    try:
        Sito = CalculateSimpleTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateSimpleTopoIndex")
    try:
        Hato = CalculateHarmonicTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateHarmonicTopoIndex")
    try:
        Geto = CalculateGeometricTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateGeometricTopoIndex")
    try:
        Arto = CalculateArithmeticTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateArithmeticTopoIndex")
        
    
    return [W, AW, J, Xu, GMTI,
           Pol, DZ, Ipc, BertzCT, Thara,
           Tsch, ZM1, ZM2, MZM1, MZM2,
           Qindex, Platt, 
           petitjeant, Sito,
           Hato, Geto]


# In[30]:


def GetAllTopologicalProperties_ligand(mol):
    try:
        W = CalculateWeiner(mol)
    except:
        print("Could Not Calculate CalculateWeiner")
    try:
        AW = CalculateMeanWeiner(mol)
    except:
        print("Could Not Calculate CalculateMeanWeiner")
    try:
        J = CalculateBalaban(mol)
    except:
        print("Could Not Calculate CalculateBalaban")
    try:
        Tigdi = CalculateGraphDistance(mol)
    except:
        print("Could Not Calculate CalculateGraphDistance")
    try:
        Xu = CalculateXuIndex(mol)
    except:
        print("Could Not Calculate CalculateXuIndex")
    try:
        GMTI = CalculateGutmanTopo(mol)
    except:
        print("Could Not Calculate CalculateGutmanTopo")
    try:
        Pol = CalculatePolarityNumber(mol)
    except:
        print("Could Not Calculate CalculatePolarityNumber")
    try:
        DZ = CalculatePoglianiIndex(mol)
    except:
        print("Could Not Calculate CalculatePoglianiIndex")
    try:
        Ipc = CalculateIpc(mol)
    except:
        print("Could Not Calculate CalculateIpc")
    try:
        BertzCT = CalculateBertzCT(mol)
    except:
        print("Could Not Calculate CalculateBertzCT")
    try:
        Thara = CalculateHarary(mol)
    except:
        print("Could Not Calculate CalculateHarary")
    try:
        Tsch = CalculateSchiultz(mol)
    except:
        print("Could Not Calculate CalculateSchiultz")
    try:
        ZM1 = CalculateZagreb1(mol)
    except:
        print("Could Not Calculate CalculateZagreb1")
    try:
        ZM2 = CalculateZagreb2(mol)
    except:
        print("Could Not Calculate CalculateZagreb2")
    try:
        MZM1 = CalculateMZagreb1(mol)
    except:
        print("Could Not Calculate CalculateMZagreb1")
    try:
        MZM2 = CalculateMZagreb2(mol)
    except:
        print("Could Not Calculate CalculateMZagreb2")
    try:
        Qindex = CalculateQuadratic(mol)
    except:
        print("Could Not Calculate CalculateQuadratic")
    try:
        Platt = CalculatePlatt(mol)
    except:
        print("Could Not Calculate CalculatePlatt")
    try:
        diametert = CalculateDiameter(mol)
    except:
        print("Could Not Calculate CalculateDiameter")
    try:
        radiust = CalculateRadius(mol)
    except:
        print("Could Not Calculate CalculateRadius")
    try:
        petitjeant = CalculatePetitjean(mol)
    except:
        print("Could Not Calculate CalculatePetitjean")
    try:
        Sito = CalculateSimpleTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateSimpleTopoIndex")
    try:
        Hato = CalculateHarmonicTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateHarmonicTopoIndex")
    try:
        Geto = CalculateGeometricTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateGeometricTopoIndex")
    try:
        Arto = CalculateArithmeticTopoIndex(mol)
    except:
        print("Could Not Calculate CalculateArithmeticTopoIndex")
        
    
    return [W, AW, J, Tigdi, Xu, GMTI,
           Pol, DZ, Ipc, BertzCT, Thara,
           Tsch, ZM1, ZM2, MZM1, MZM2,
           Qindex, Platt, diametert,
           radiust, petitjeant, Sito,
           Hato, Geto, Arto]


# In[32]:


def GetTopology(mol, ligand = True):
    
    if ligand:
        _Topology={'W':CalculateWeiner,
               'AW':CalculateMeanWeiner,
               'J':CalculateBalaban,
               'Tigdi':CalculateGraphDistance,
               'Xu':CalculateXuIndex,
               'GMTI':CalculateGutmanTopo,
               'Pol':CalculatePolarityNumber,
               'DZ':CalculatePoglianiIndex,
               'Ipc':CalculateIpc,
               'BertzCT':CalculateBertzCT,
               'Thara':CalculateHarary,
               'Tsch':CalculateSchiultz,
               'ZM1':CalculateZagreb1,
               'ZM2':CalculateZagreb2,
               'MZM1':CalculateMZagreb1,
               'MZM2':CalculateMZagreb2,
               'Qindex':CalculateQuadratic,
               'Platt':CalculatePlatt,
               'diametert':CalculateDiameter,
               'radiust':CalculateRadius,
               'petitjeant':CalculatePetitjean,
               'Sito':CalculateSimpleTopoIndex,
               'Hato':CalculateHarmonicTopoIndex,
               'Geto':CalculateGeometricTopoIndex,
               'Arto':CalculateArithmeticTopoIndex      
        }
    else:
        _Topology={'W':CalculateWeiner,
               'AW':CalculateMeanWeiner,
               'J':CalculateBalaban,
               'Xu':CalculateXuIndex,
               'GMTI':CalculateGutmanTopo,
               'Pol':CalculatePolarityNumber,
               'DZ':CalculatePoglianiIndex,
               'Ipc':CalculateIpc,
               'BertzCT':CalculateBertzCT,
               'Thara':CalculateHarary,
               'Tsch':CalculateSchiultz,
               'ZM1':CalculateZagreb1,
               'ZM2':CalculateZagreb2,
               'MZM1':CalculateMZagreb1,
               'MZM2':CalculateMZagreb2,
               'Qindex':CalculateQuadratic,
               'Platt':CalculatePlatt,
               'petitjeant':CalculatePetitjean,
               'Sito':CalculateSimpleTopoIndex,
               'Hato':CalculateHarmonicTopoIndex,
               'Geto':CalculateGeometricTopoIndex,
               'Arto':CalculateArithmeticTopoIndex      
        }
    result={}
    for DesLabel in _Topology.keys():
        try:
            result[DesLabel]=round(_Topology[DesLabel](mol),3)
        except:
            result[DesLabel]=numpy.nan
    return result


# In[ ]:




