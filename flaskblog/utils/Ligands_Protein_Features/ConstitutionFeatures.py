#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rdkit import Chem
#from rdkit.Chem import rdchem
from rdkit.Chem import Lipinski as LPK


# In[2]:


def CalculateAverageMolWeight(mol, num_atoms=1):
    MolWeight=0
    for atom in mol.GetAtoms():
        MolWeight=MolWeight+atom.GetMass()

    return MolWeight/num_atoms


# In[3]:


def CalculateHydrogenNumber(mol, num_atoms=1):
    i=0
    Hmol=Chem.AddHs(mol)
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==1:
            i=i+1
            
    return i/num_atoms


# In[4]:


def CalculateHalogenNumber(mol, num_atoms=1):
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==9 or atom.GetAtomicNum()==17 or atom.GetAtomicNum()==35 or atom.GetAtomicNum()==53:
            i=i+1
    return i/num_atoms


# In[5]:


def CalculateHeteroNumber(mol, num_atoms=1):
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==6 or atom.GetAtomicNum()==1:
            i=i+1

    return (num_atoms - i)/num_atoms


# In[1]:


def CalculateHeavyAtomNumber(mol,num_atoms=1):
    return mol.GetNumHeavyAtoms()/num_atoms


# In[7]:


def _CalculateElementNumber(mol,AtomicNumber=6):
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNumber:
            i=i+1
            
    return i


# In[8]:


def CalculateFluorinNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=9) / num_atoms


# In[9]:


def CalculateChlorinNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=17) / num_atoms


# In[10]:


def CalculateBromineNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=35) / num_atoms


# In[11]:


def CalculateIodineNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=53) / num_atoms


# In[12]:


def CalculateCarbonNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=6) / num_atoms


# In[13]:


def CalculatePhosphorNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=15) / num_atoms


# In[14]:


def CalculateSulfurNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=16) / num_atoms


# In[15]:


def CalculateOxygenNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=8) / num_atoms


# In[16]:


def CalculateNitrogenNumber(mol, num_atoms=1):
    return _CalculateElementNumber(mol,AtomicNumber=7) / num_atoms


# In[17]:


def CalculateRingNumber(mol):
    return Chem.GetSSSR(mol)


# In[18]:


def CalculateRotationBondNumber(mol, num_bonds=1):
    return LPK.NumRotatableBonds(mol) / num_bonds


# In[19]:


def CalculateHdonorNumber(mol, relative_num_h=1, num_atoms=1):
    num_h_atoms = relative_num_h * num_atoms
    return LPK.NumHDonors(mol) / num_h_atoms


# In[20]:


def CalculateHacceptorNumber(mol, relative_num_h=1, num_atoms=1):
    num_h_atoms = relative_num_h * num_atoms
    return LPK.NumHAcceptors(mol) / num_h_atoms


# In[21]:


def CalculateSingleBondNumber(mol, bonds, num_bonds=1):
    i=0;
    for bond in bonds:

        if bond.GetBondType().name=='SINGLE':
            i=i+1
            
    return i / num_bonds


# In[22]:


def CalculateDoubleBondNumber(mol, bonds, num_bonds=1):
    i=0;
    for bond in bonds:

        if bond.GetBondType().name=='DOUBLE':
            i=i+1
            
    return i / num_bonds


# In[23]:


def CalculateTripleBondNumber(mol, bonds, num_bonds=1):
    i=0;
    for bond in bonds:

        if bond.GetBondType().name=='TRIPLE':
            i=i+1
            
    return i / num_bonds


# In[24]:


def CalculateAromaticBondNumber(mol, bonds, num_bonds=1):
    i=0;
    for bond in bonds:

        if bond.GetBondType().name=='AROMATIC':
            i=i+1
            
    return i


# In[25]:


def _CalculatePathN(mol,PathLength=2):
    return len(Chem.FindAllPathsOfLengthN(mol,PathLength,useBonds=1))


# In[26]:


def CalculatePath1(mol):
    return _CalculatePathN(mol,1)


# In[27]:


def CalculatePath2(mol):
    return _CalculatePathN(mol,2)


# In[28]:


def CalculatePath3(mol):
    return _CalculatePathN(mol,3)


# In[29]:


def CalculatePath4(mol):
    return _CalculatePathN(mol,4)


# In[30]:


def CalculatePath5(mol):
    return _CalculatePathN(mol,5)


# In[31]:


def CalculatePath6(mol):
    return _CalculatePathN(mol,6)


# In[32]:


def GetAllConstitutionProperties_pocket(mol):
    
    num_atoms = mol.GetNumAtoms()
    bonds = mol.GetBonds()
    num_bonds = mol.GetNumBonds()
    AWeight = CalculateAverageMolWeight(mol, num_atoms)
    nhyd = CalculateHydrogenNumber(mol, num_atoms)
    nhal = CalculateHalogenNumber(mol, num_atoms)
    nhet = CalculateHeteroNumber(mol, num_atoms)
    nhev = CalculateHeavyAtomNumber(mol, num_atoms)
    ncof = CalculateFluorinNumber(mol, num_atoms)
    ncocl = CalculateChlorinNumber(mol, num_atoms)
    ncobr = CalculateBromineNumber(mol, num_atoms)
    ncoi = CalculateIodineNumber(mol, num_atoms)
    ncarb = CalculateCarbonNumber(mol, num_atoms)
    nphos = CalculatePhosphorNumber(mol, num_atoms)
    nsulph = CalculateSulfurNumber(mol, num_atoms)
    noxy = CalculateOxygenNumber(mol, num_atoms)
    nnitro = CalculateNitrogenNumber(mol, num_atoms)
    nring = CalculateRingNumber(mol)
    nrot = CalculateRotationBondNumber(mol, num_bonds)
    ndonr = CalculateHdonorNumber(mol, num_atoms, nhyd)
    naccr = CalculateHacceptorNumber(mol, num_atoms, nhyd)
    nsb = CalculateSingleBondNumber(mol, num_bonds, bonds)
    ndb = CalculateDoubleBondNumber(mol, num_bonds, bonds)
    naro = CalculateAromaticBondNumber(mol, num_bonds, bonds)
    ntb = CalculateTripleBondNumber(mol, num_bonds, bonds)
    PC1 = CalculatePath1(mol)
    PC2 = CalculatePath2(mol)
    PC3 = CalculatePath3(mol)
    PC4 = CalculatePath4(mol)
    PC5 = CalculatePath5(mol)
    PC6 = CalculatePath6(mol)
    
    return [
        AWeight,
        nhyd,
        nhal,
        nhet,
        nhev,
        ncof,
        ncocl,
        ncobr,
        ncoi,
        ncarb,
        nphos,
        nsulph,
        noxy,
        nnitro,
        nring,
        nrot,
        ndonr,
        naccr,
        nsb,
        ndb,
        naro,
        ntb,
        PC1,
        PC2,
        PC3,
        PC4,
        PC5,
        PC6
    ]


# In[33]:


def GetAllConstitutionProperties_ligand(mol):
    
    num_atoms = mol.GetNumAtoms()
    bonds = mol.GetBonds()
    num_bonds = mol.GetNumBonds()
    AWeight = CalculateAverageMolWeight(mol, num_atoms)
    nhyd = CalculateHydrogenNumber(mol, num_atoms)
    nhal = CalculateHalogenNumber(mol, num_atoms)
    nhet = CalculateHeteroNumber(mol, num_atoms)
    nhev = CalculateHeavyAtomNumber(mol, num_atoms)
    ncof = CalculateFluorinNumber(mol, num_atoms)
    ncocl = CalculateChlorinNumber(mol, num_atoms)
    ncobr = CalculateBromineNumber(mol, num_atoms)
    ncoi = CalculateIodineNumber(mol, num_atoms)
    ncarb = CalculateCarbonNumber(mol, num_atoms)
    nphos = CalculatePhosphorNumber(mol, num_atoms)
    nsulph = CalculateSulfurNumber(mol, num_atoms)
    noxy = CalculateOxygenNumber(mol, num_atoms)
    nnitro = CalculateNitrogenNumber(mol, num_atoms)
    nring = CalculateRingNumber(mol)
    nrot = CalculateRotationBondNumber(mol, num_bonds)
    ndonr = CalculateHdonorNumber(mol, num_atoms, nhyd)
    naccr = CalculateHacceptorNumber(mol, num_atoms, nhyd)
    nsb = CalculateSingleBondNumber(mol, num_bonds, bonds)
    ndb = CalculateDoubleBondNumber(mol, num_bonds, bonds)
    naro = CalculateAromaticBondNumber(mol, num_bonds, bonds)
    ntb = CalculateTripleBondNumber(mol, num_bonds, bonds)
    PC1 = CalculatePath1(mol)
    PC2 = CalculatePath2(mol)
    PC3 = CalculatePath3(mol)
    PC4 = CalculatePath4(mol)
    PC5 = CalculatePath5(mol)
    PC6 = CalculatePath6(mol)
    
    return [
        AWeight,
        nhyd,
        nhal,
        nhet,
        nhev,
        ncof,
        ncocl,
        ncobr,
        ncoi,
        ncarb,
        nphos,
        nsulph,
        noxy,
        nnitro,
        nring,
        nrot,
        ndonr,
        naccr,
        nsb,
        ndb,
        naro,
        ntb,
        PC1,
        PC2,
        PC3,
        PC4,
        PC5,
        PC6
    ]


# In[34]:


def GetConstitutional_relative(mol):
    
    num_atoms = mol.GetNumAtoms()
    bonds = mol.GetBonds()
    num_bonds = mol.GetNumBonds()
    nhyd = CalculateHydrogenNumber(mol, num_atoms)
    
    result={  
        "AWeight": CalculateAverageMolWeight(mol, num_atoms),
        "nhyd": nhyd,
        "nhal": CalculateHalogenNumber(mol, num_atoms),
        "nhet": CalculateHeteroNumber(mol, num_atoms),
        "nhev": CalculateHeavyAtomNumber(mol, num_atoms),
        "ncof": CalculateFluorinNumber(mol, num_atoms),
        "ncocl": CalculateChlorinNumber(mol, num_atoms),
        "ncobr": CalculateBromineNumber(mol, num_atoms),
        "ncoi": CalculateIodineNumber(mol, num_atoms),
        "ncarb": CalculateCarbonNumber(mol, num_atoms),
        "nphos": CalculatePhosphorNumber(mol, num_atoms),
        "nsulph": CalculateSulfurNumber(mol, num_atoms),
        "noxy": CalculateOxygenNumber(mol, num_atoms),
        "nnitro": CalculateNitrogenNumber(mol, num_atoms),
        "nring": CalculateRingNumber(mol),
        "nrot": CalculateRotationBondNumber(mol, num_bonds),
        "ndonr": CalculateHdonorNumber(mol, num_atoms, nhyd),
        "naccr": CalculateHacceptorNumber(mol, num_atoms, nhyd),
        "nsb": CalculateSingleBondNumber(mol, num_bonds, bonds),
        "ndb": CalculateDoubleBondNumber(mol, num_bonds, bonds),
        "naro": CalculateAromaticBondNumber(mol, num_bonds, bonds),
        "ntb": CalculateTripleBondNumber(mol, num_bonds, bonds),
        "PC1": CalculatePath1(mol),
        "PC2": CalculatePath2(mol),
        "PC3": CalculatePath3(mol),
        "PC4": CalculatePath4(mol),
        "PC5": CalculatePath5(mol),
        "PC6": CalculatePath6(mol)
    }
    
    
    return result


def GetConstitutional(mol):
    
    num_atoms = mol.GetNumAtoms()
    bonds = mol.GetBonds()
    num_bonds = mol.GetNumBonds()
    nhyd = CalculateHydrogenNumber(mol)
    
    result={  
        "AWeight": CalculateAverageMolWeight(mol),
        "nhyd": nhyd,
        "nhal": CalculateHalogenNumber(mol),
        "nhet": CalculateHeteroNumber(mol),
        "nhev": CalculateHeavyAtomNumber(mol),
        "ncof": CalculateFluorinNumber(mol),
        "ncocl": CalculateChlorinNumber(mol),
        "ncobr": CalculateBromineNumber(mol),
        "ncoi": CalculateIodineNumber(mol),
        "ncarb": CalculateCarbonNumber(mol),
        "nphos": CalculatePhosphorNumber(mol),
        "nsulph": CalculateSulfurNumber(mol),
        "noxy": CalculateOxygenNumber(mol),
        "nnitro": CalculateNitrogenNumber(mol),
        "nring": CalculateRingNumber(mol),
        "nrot": CalculateRotationBondNumber(mol),
        "ndonr": CalculateHdonorNumber(mol),
        "naccr": CalculateHacceptorNumber(mol),
        "nsb": CalculateSingleBondNumber(mol, bonds),
        "ndb": CalculateDoubleBondNumber(mol, bonds),
        "naro": CalculateAromaticBondNumber(mol, bonds),
        "ntb": CalculateTripleBondNumber(mol, bonds),
        "PC1": CalculatePath1(mol),
        "PC2": CalculatePath2(mol),
        "PC3": CalculatePath3(mol),
        "PC4": CalculatePath4(mol),
        "PC5": CalculatePath5(mol),
        "PC6": CalculatePath6(mol)
    }
    
    
    return result



# In[ ]:




