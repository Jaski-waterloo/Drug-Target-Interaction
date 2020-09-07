#!/usr/bin/env python
# coding: utf-8

# In[1]:


def compute_charge(mol):
    from rdkit.Chem import AllChem
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception as e:
        logging.exception("Unable to Calculate charge for molcule")
        
    return mol


# In[2]:


def get_partial_charge(atom):
    try:
        value = atom.GetProp(str("_GasteigerCharge"))
        if value == '-nan':
            return 0
        return float(value)
    except KeyError:
        return 0


# In[3]:


def sanitize_mol(mol, ID = ""):
    from rdkit import Chem
    
    try:
        Chem.SanitizeMol(mol)
    except:
        raise ValueError("Can Not Sanitize Mol :", ID)


# In[4]:


def get_coordinates(mol):
    xyz = np.zeros((mol.GetNumAtoms(), 3))
    conf = mol.GetConformer()
    
    for i in range(conf.GetNumAtoms()):
        position = conf.GetAtomPosition(i)
        xyz[i,0] = position.x
        xyz[i,1] = position.y
        xyz[i,2] = position.z
        
    return(xyz)


# In[5]:


def compute_pairwise_distances(protein_xyz, ligand_xyz):
    #Numpy arrays
    from scipy.spatial.distance import cdist
    pairwise_distances = cdist(protein_xyz, ligand_xyz, metric='euclidean')
    return (pairwise_distances)


# In[6]:


import numpy as np


# In[7]:


def compute_morgan_fingerprint(mol, degree=3, bits=1024):
    from rdkit.Chem import AllChem
    bv = AllChem.GetMorganFingerprintAsBitVect(
          mol, degree, nBits=bits)
    return np.array(bv)


# In[8]:


def load_protein(pdb_file):
    #Make sure the pdb file is fixed
    from rdkit import Chem
    my_mol = Chem.MolFromPDBFile(
        str(pdb_file), sanitize=True, removeHs=False)
    
    if my_mol is None:
        raise ValueError("Unable to Read PDB File")
        
    #compute_charge(my_mol)
    
    #xyz = get_coordinates(my_mol)
    
    return xyz, my_mol


# In[ ]:





# In[9]:


def GetAtomResidueId(atom):
    """Return (residue number, residue name, chain id) for a given atom"""
    info = atom.GetPDBResidueInfo()
    res_id = (info.GetResidueNumber(), info.GetResidueName().strip(),
              info.GetChainId())
    return res_id


# In[10]:


def AtomListToSubMol(mol, amap, includeConformer=False):
    from rdkit import Chem
    from itertools import chain, combinations

    """
    Parameters
    ----------
        mol: rdkit.Chem.rdchem.Mol
            Molecule
        amap: array-like
            List of atom indices (zero-based)
        includeConformer: bool (default=True)
            Toogle to include atoms coordinates in submolecule.

    Returns
    -------
        submol: rdkit.Chem.rdchem.RWMol
            Submol determined by specified atom list
    """
    if not isinstance(amap, list):
        amap = list(amap)
    submol = Chem.RWMol(Chem.Mol())
    for aix in amap:
        submol.AddAtom(mol.GetAtomWithIdx(aix))
    for i, j in combinations(amap, 2):
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond:
            submol.AddBond(amap.index(i),
                           amap.index(j),
                           bond.GetBondType())
    if includeConformer:
        for conf in mol.GetConformers():
            new_conf = Chem.Conformer(len(amap))
            for i in range(len(amap)):
                new_conf.SetAtomPosition(i, conf.GetAtomPosition(amap[i]))
                new_conf.SetId(conf.GetId())
                new_conf.Set3D(conf.Is3D())
            submol.AddConformer(new_conf)
    return submol


# In[11]:


def extract_pocket_ligand(pdb_file, cutoff=12., blacklisted_ligands = []):
    """
    It Takes PDB File as input and returns a list
    The list contains the pocket and corresponding ligand
    Takes into consideration all the atoms of the residue present in binding site
    """
    # Get heteroatom residues - connectivity still might be wrong, so GetFrags will fail
    
    import numpy as np
    from scipy.spatial.distance import cdist, chain
    
    mol = Chem.MolFromPDBFile(pdb_file, sanitize = False)
    
    #Molecule needs to be sanitized
    #try:
    #    Chem.SanitizeMol(mol)
    #except:
    #    raise ValueError("Molecule can not be sanitized exiting")
        
    sanitize_mol(mol, pdb_file)
    
    hetatm_residues = {}
    protein_residues = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_id = GetAtomResidueId(atom)
        if info.GetIsHeteroAtom():
            if res_id not in hetatm_residues:
                hetatm_residues[res_id] = []
            hetatm_residues[res_id].append(atom.GetIdx())
        else:
            if res_id not in protein_residues:
                protein_residues[res_id] = []
            protein_residues[res_id].append(atom.GetIdx())


    for res_id in list(hetatm_residues.keys()):  # exhaust keys since we modify
        # Treat single atom residues (waters + metals) as pocket residues
        if len(hetatm_residues[res_id]) == 1:
            protein_residues[res_id] = hetatm_residues[res_id]
            del hetatm_residues[res_id]

    if len(hetatm_residues) == 0:
        raise ValueError('No ligands')

    
    #interaction_dict = {
    #    "pocket": None,
    #    "pocket_residues": None,
    #    "ligand_PDB_ID": None
    #}
    
    interaction_list = []
    
    hetatm_residues_keys = hetatm_residues.keys()
    pprint(hetatm_residues_keys)
    for ligand_key in hetatm_residues_keys:
        #print(i)
        #ligand_key = hetatm_residues_keys[i]
        #print(ligand_key)
        if ligand_key in blacklisted_ligands:
            print(ligand_key, "Is Blacklisted")
            continue
        ligand_amap = hetatm_residues[ligand_key]
        ligand = AtomListToSubMol(mol, ligand_amap, includeConformer=True)
        # should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = ligand.GetConformer()
        ligand_coords = np.array([conf.GetAtomPosition(i)
                                  for i in range(ligand.GetNumAtoms())])

        # Get protein and waters
        blacklist_ids = list(chain(*hetatm_residues.values()))
        protein_amap = np.array([i for i in range(mol.GetNumAtoms())
                                 if i not in blacklist_ids])
        # should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = mol.GetConformer()
        protein_coords = np.array([conf.GetAtomPosition(i)
                                  for i in protein_amap.tolist()])

        # Pocket selection based on cutoff
        mask = (cdist(protein_coords, ligand_coords) <= cutoff).any(axis=1)
        # IDs of atoms within cutoff
        pocket_amap = protein_amap[np.where(mask)[0]].tolist()

        # Expand pocket's residues
        pocket_residues = {}
        for res_id in protein_residues.keys():
            if any(1 for res_aix in protein_residues[res_id]
                    if res_aix in pocket_amap):
                pocket_residues[res_id] = protein_residues[res_id]
        pocket_amap = list(chain(*pocket_residues.values()))

        # Create pocket mol, pocket_amap needs to be mapped to mol Idxs
        pocket = AtomListToSubMol(mol, pocket_amap, includeConformer=True)
        
        #interaction_dict.update({
        #    "pocket": pocket,
        #    "pocket_residues": pocket_residues,
        #    "ligand_PDB_ID": ligand_key[1]
        #})
        
        #try:
        #    Chem.SanitizeMol(pocket)
        #except:
        #    print("Could Not Sanitize Pocket for Ligand: ", ligand_key[1])
        #    continue
        
        #sanitize_mol(pocket, ligand_key[1])
        #compute_charges(pocket)
        
        interaction_list.append([pocket, pocket_residues, ligand_key[1]])

    return interaction_list


# In[12]:


def get_partial_charge(atom):
    try:
        charge = atom.GetProp(str("_GasteigerCharge"))
        if value == '-nan':
            return 0
        return float(charge)
    except:
        return 0


# In[13]:


def get_residue_map(pdb_file, get_residues = True):
    from biopandas.pdb import PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_file)
    files = ppdb.pdb_text.split('\n')[:1000]
    
    sites = {}
    residues = {
        "Site Name":[],
        "Residue Name":{},
        "Residue Chain":{},
        "Residue Number":{}
    }
    

    temp = [i for i in files if "REMARK 800" in i]

    j = 1
    for i in temp:
        tempp = i.split()
        if len(tempp) > 7:
            #print(tempp[-3])
            if len(tempp[-1]) > 4:
                tempp.append(tempp[-1][1:])
                tempp[-2] = tempp[-2][:1]
            site_name = 'AC'+str(j)
            temp_site_key = (int(tempp[-1]), tempp[-3], tempp[-2])
            t = {site_name: temp_site_key}
            #sites.append(['AC'+str(j), tempp[-3]])
            sites.update(t)
            j += 1
            residues["Site Name"].append(temp_site_key)
            residues["Residue Name"].update({temp_site_key:[]})
            residues["Residue Chain"].update({temp_site_key:[]})
            residues["Residue Number"].update({temp_site_key:[]})

    if get_residues: 
        return sites, residues
    else:
        return sites


# In[14]:


def get_residue_dictionary(pdb_file, sites, residues):
    from biopandas.pdb import PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_file)
    files = ppdb.pdb_text.split('\n')[:1000]
    temp = [i for i in files if "SITE" in i and "REMARK" not in i]

    for i in temp:
        tempp = i.split()
        flag = 0
        for i in temp:
            if len(i) > 4:
                flag = 1
                break
        if flag == 1:
            t = []
            for i in tempp:
                if len(i) > 4:
                    t.append(i[0])
                    t.append(i[1:])
                else:
                    t.append(i)
            tempp = t
        for j,k,l in zip(tempp[4::3], tempp[5::3], tempp[6::3]):
            residues["Residue Name"][sites[tempp[2]]].append(j)
            residues["Residue Chain"][sites[tempp[2]]].append(k)
            residues["Residue Number"][sites[tempp[2]]].append(l)
            
    return residues


# In[15]:


def get_residue_keys(residues):
    residues_keys = {}
    from pprint import pprint
    for i in residues["Site Name"]:
        num = len(residues["Residue Number"][i])
        temp = [(int(residues["Residue Number"][i][j]), residues["Residue Name"][i][j],
                 residues["Residue Chain"][i][j]) for j in range(num)]
        temp_dic = {i:temp}
        residues_keys.update(temp_dic)
        del temp
        del temp_dic
    
    return residues_keys


# In[16]:


def extract_pocket_ligand_union(pdb_file, cutoff=12., blacklisted_ligands = []):
    """
    It Takes PDB File as input and returns a list
    The list contains the pocket and corresponding ligand
    Takes into consideration all the atoms of the residue present in binding site
    """
    from rdkit import Chem
    from pprint import pprint
    import numpy as np
    from scipy.spatial.distance import cdist
    from itertools import chain, combinations
    from rdkit.Chem import AllChem
    import pickle

    try:
        mol = Chem.MolFromPDBFile(pdb_file, sanitize = False)
    except:
        try:
            mol = Chem.MolFromPDBBlock(pdb_file)
        except:
            mol = pdb_file
    
    #Molecule needs to be sanitized
    try:
        Chem.SanitizeMol(mol)
    except:
        print("Molecule can not be sanitized exiting")
        return
    
    try:
        sites, residues = get_residue_map(pdb_file)
        residues = get_residue_dictionary(pdb_file, sites, residues)
        ligand_residue_dict = get_residue_keys(residues)
        print("taken union")
    except:
        pass
    
    
        
    # Get heteroatom residues - connectivity still might be wrong, so GetFrags will fail
    
    hetatm_residues = {}
    protein_residues = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_id = GetAtomResidueId(atom)
        if info.GetIsHeteroAtom():
            if res_id not in hetatm_residues:
                hetatm_residues[res_id] = []
            hetatm_residues[res_id].append(atom.GetIdx())
        else:
            if res_id not in protein_residues:
                protein_residues[res_id] = []
            protein_residues[res_id].append(atom.GetIdx())


    for res_id in list(hetatm_residues.keys()):  # exhaust keys since we modify
        # Treat single atom residues (waters + metals) as pocket residues
        if len(hetatm_residues[res_id]) == 1:
            protein_residues[res_id] = hetatm_residues[res_id]
            del hetatm_residues[res_id]

    if len(hetatm_residues) == 0:
        raise ValueError('No ligands')

    
    
    interaction_list = []
    
    hetatm_residues_keys = hetatm_residues.keys()
    pprint(hetatm_residues_keys)
    for ligand_key in hetatm_residues_keys:
        if ligand_key[1] in blacklisted_ligands:
            print(ligand_key, "Is Blacklisted")
            continue
        ligand_amap = hetatm_residues[ligand_key]
        ligand = AtomListToSubMol(mol, ligand_amap, includeConformer=True)
        # we should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = ligand.GetConformer()
        ligand_coords = np.array([conf.GetAtomPosition(i)
                                  for i in range(ligand.GetNumAtoms())])

        # Get protein and waters
        blacklist_ids = list(chain(*hetatm_residues.values()))
        protein_amap = np.array([i for i in range(mol.GetNumAtoms())
                                 if i not in blacklist_ids])
        # we should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = mol.GetConformer()
        protein_coords = np.array([conf.GetAtomPosition(i)
                                  for i in protein_amap.tolist()])

        # Pocket selection based on cutoff
        mask = (cdist(protein_coords, ligand_coords) <= cutoff).any(axis=1)
        # IDs of atoms within cutoff
        
        pocket_amap = []
        try:
            pocket_amap = set(protein_amap[np.where(mask)[0]].tolist())
            ligand_residue_keys = ligand_residue_dict[ligand_key]

            pocket_amap.union(set(list(chain(*[protein_residues[ligand_residue_keys[i]] 
                                   for i in range(len(ligand_residue_keys))]))))
            pocket_amap = list(pocket_amap)
        except:
            pocket_amap = protein_amap[np.where(mask)[0]].tolist()
        
        

        # Expand pocket's residues
        pocket_residues = {}
        for res_id in protein_residues.keys():
            if any(1 for res_aix in protein_residues[res_id]
                    if res_aix in pocket_amap):
                pocket_residues[res_id] = protein_residues[res_id]
        pocket_amap = list(chain(*pocket_residues.values()))

        # Create pocket mol, pocket_amap needs to be mapped to mol Idxs
        pocket = AtomListToSubMol(mol, pocket_amap, includeConformer=True)
        

        try:
            Chem.SanitizeMol(pocket)
        except:
            print("Could Not Sanitize Pocket for Ligand: ", ligand_key[1])
            continue
            
        AllChem.ComputeGasteigerCharges(pocket)
        
        interaction_list.append([pocket, ligand_key[0], ligand_key[1], ligand_key[2], pdb_file])

        #pickle.dump([pocket, ligand_key[1], pdb_file], pickle_file)

    return interaction_list


# In[5]:


def generate_binding_site(ligand_key,mol,hetatm_residues,protein_residues,cutoff=7.0, blacklisted_ligands = []):
    
    from scipy.spatial.distance import cdist
    from itertools import chain, combinations
    
    if ligand_key[1] in blacklisted_ligands:
        print(ligand_key, "Is Blacklisted")
        return
    ligand_amap = hetatm_residues[ligand_key]
    ligand = AtomListToSubMol(mol, ligand_amap, includeConformer=True)
    # we should use GetPositions() here, but it often leads to segfault (RDKit)
    conf = ligand.GetConformer()
    ligand_coords = np.array([conf.GetAtomPosition(i)
                              for i in range(ligand.GetNumAtoms())])
    # Get protein and waters
    blacklist_ids = list(chain(*hetatm_residues.values()))
    protein_amap = np.array([i for i in range(mol.GetNumAtoms())
                             if i not in blacklist_ids])
    # we should use GetPositions() here, but it often leads to segfault (RDKit)
    conf = mol.GetConformer()
    protein_coords = np.array([conf.GetAtomPosition(i)
                              for i in protein_amap.tolist()])
    # Pocket selection based on cutoff
    mask = (cdist(protein_coords, ligand_coords) <= cutoff).any(axis=1)
    # IDs of atoms within cutoff
    
    pocket_amap = []
    try:
        pocket_amap = set(protein_amap[np.where(mask)[0]].tolist())
        ligand_residue_keys = ligand_residue_dict[ligand_key]
        pocket_amap.union(set(list(chain(*[protein_residues[ligand_residue_keys[i]] 
                               for i in range(len(ligand_residue_keys))]))))
        pocket_amap = list(pocket_amap)
    except:
        pocket_amap = protein_amap[np.where(mask)[0]].tolist()
    
    
    # Expand pocket's residues
    pocket_residues = {}
    for res_id in protein_residues.keys():
        if any(1 for res_aix in protein_residues[res_id]
                if res_aix in pocket_amap):
            pocket_residues[res_id] = protein_residues[res_id]
    pocket_amap = list(chain(*pocket_residues.values()))
    # Create pocket mol, pocket_amap needs to be mapped to mol Idxs
    pocket = AtomListToSubMol(mol, pocket_amap, includeConformer=True)
    
    try:
        Chem.SanitizeMol(pocket)
    except:
        print("Could Not Sanitize Pocket for Ligand: ", ligand_key[1])
        return
        
    AllChem.ComputeGasteigerCharges(pocket)
    
    return [pocket, pocket_residues, ligand_key[1]]


# In[2]:


def extract_binding_sites_multiprocessing(pdb_file, cutoff=12., blacklisted_ligands = []):
    
    from rdkit import Chem
    from pprint import pprint
    import numpy as np
    from scipy.spatial.distance import cdist
    from itertools import chain, combinations
    from rdkit.Chem import AllChem
    
    sites, residues = get_residue_map(pdb_file)
    residues = get_residue_dictionary(pdb_file, sites, residues)
    ligand_residue_dict = get_residue_keys(residues)
    
    mol = Chem.MolFromPDBFile(pdb_file, sanitize = False)
    
    #Molecule needs to be sanitized
    try:
        Chem.SanitizeMol(mol)
    except:
        raise ValueError("Molecule can not be sanitized exiting")
        
    # Get heteroatom residues - connectivity still might be wrong, so GetFrags will fail
    
    hetatm_residues = {}
    protein_residues = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_id = GetAtomResidueId(atom)
        if info.GetIsHeteroAtom():
            if res_id not in hetatm_residues:
                hetatm_residues[res_id] = []
            hetatm_residues[res_id].append(atom.GetIdx())
        else:
            if res_id not in protein_residues:
                protein_residues[res_id] = []
            protein_residues[res_id].append(atom.GetIdx())


    for res_id in list(hetatm_residues.keys()):  # exhaust keys since we modify
        # Treat single atom residues (waters + metals) as pocket residues
        if len(hetatm_residues[res_id]) == 1:
            protein_residues[res_id] = hetatm_residues[res_id]
            del hetatm_residues[res_id]

    if len(hetatm_residues) == 0:
        raise ValueError('No ligands')

    
    
    interaction_list = []
    
    hetatm_residues_keys = hetatm_residues.keys()
    pprint(hetatm_residues_keys)
    
    


def get_residue_list(pdb_file, cutoff=4., blacklisted_ligands = ['HOH']):
    """
    It Takes PDB File as input and returns a list
    The list contains the pocket and corresponding ligand
    Takes into consideration all the atoms of the residue present in binding site
    """
    from rdkit import Chem
    from pprint import pprint
    import numpy as np
    from scipy.spatial.distance import cdist
    from itertools import chain, combinations
    from rdkit.Chem import AllChem
    import pickle

    try:
        mol = Chem.MolFromPDBFile(pdb_file, sanitize = False)
    except:
        try:
            mol = Chem.MolFromPDBBlock(pdb_file)
        except:
            mol = pdb_file
    
    #Molecule needs to be sanitized
    try:
        Chem.SanitizeMol(mol)
    except:
        print("Molecule can not be sanitized exiting")
        return
    
    try:
        sites, residues = get_residue_map(pdb_file)
        residues = get_residue_dictionary(pdb_file, sites, residues)
        ligand_residue_dict = get_residue_keys(residues)
        print("taken union")
    except:
        pass
    
    
        
    # Get heteroatom residues - connectivity still might be wrong, so GetFrags will fail
    
    hetatm_residues = {}
    protein_residues = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_id = GetAtomResidueId(atom)
        if info.GetIsHeteroAtom():
            if res_id not in hetatm_residues:
                hetatm_residues[res_id] = []
            hetatm_residues[res_id].append(atom.GetIdx())
        else:
            if res_id not in protein_residues:
                protein_residues[res_id] = []
            protein_residues[res_id].append(atom.GetIdx())


    for res_id in list(hetatm_residues.keys()):  # exhaust keys since we modify
        # Treat single atom residues (waters + metals) as pocket residues
        if len(hetatm_residues[res_id]) == 1:
            protein_residues[res_id] = hetatm_residues[res_id]
            del hetatm_residues[res_id]

    if len(hetatm_residues) == 0:
        raise ValueError('No ligands')

    
    
    interaction_list = []
    
    hetatm_residues_keys = hetatm_residues.keys()
    pprint(hetatm_residues_keys)
    for ligand_key in hetatm_residues_keys:
        if ligand_key[1] in blacklisted_ligands:
            print(ligand_key, "Is Blacklisted")
            continue
        else:
            interaction_list.append(ligand_key)

    return interaction_list



def get_pocket_with_id(pdb_file, key, cutoff=4., blacklisted_ligands = []):
    """
    It Takes PDB File as input and returns a list
    The list contains the pocket and corresponding ligand
    Takes into consideration all the atoms of the residue present in binding site
    """
    from rdkit import Chem
    from pprint import pprint
    import numpy as np
    from scipy.spatial.distance import cdist
    from itertools import chain, combinations
    from rdkit.Chem import AllChem
    import pickle

    try:
        mol = Chem.MolFromPDBFile(pdb_file, sanitize = False)
    except:
        try:
            mol = Chem.MolFromPDBBlock(pdb_file)
        except:
            mol = pdb_file

    #Molecule needs to be sanitized
    try:
        Chem.SanitizeMol(mol)
    except:
        print("Molecule can not be sanitized exiting")
        return

    try:
        sites, residues = get_residue_map(pdb_file)
        residues = get_residue_dictionary(pdb_file, sites, residues)
        ligand_residue_dict = get_residue_keys(residues)
        print("taken union")
    except:
        pass



    # Get heteroatom residues - connectivity still might be wrong, so GetFrags will fail

    hetatm_residues = {}
    protein_residues = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        res_id = GetAtomResidueId(atom)
        if info.GetIsHeteroAtom():
            if res_id not in hetatm_residues:
                hetatm_residues[res_id] = []
            hetatm_residues[res_id].append(atom.GetIdx())
        else:
            if res_id not in protein_residues:
                protein_residues[res_id] = []
            protein_residues[res_id].append(atom.GetIdx())


    for res_id in list(hetatm_residues.keys()):  # exhaust keys since we modify
        # Treat single atom residues (waters + metals) as pocket residues
        if len(hetatm_residues[res_id]) == 1:
            protein_residues[res_id] = hetatm_residues[res_id]
            del hetatm_residues[res_id]

    if len(hetatm_residues) == 0:
        raise ValueError('No ligands')




    hetatm_residues_keys = hetatm_residues.keys()
    pprint(hetatm_residues_keys)
    for ligand_key in hetatm_residues_keys:
        if ligand_key != key:
            print(ligand_key, "Is Not to be taken")
            continue
        ligand_amap = hetatm_residues[ligand_key]
        ligand = AtomListToSubMol(mol, ligand_amap, includeConformer=True)
        # we should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = ligand.GetConformer()
        ligand_coords = np.array([conf.GetAtomPosition(i)
                                  for i in range(ligand.GetNumAtoms())])

        # Get protein and waters
        blacklist_ids = list(chain(*hetatm_residues.values()))
        protein_amap = np.array([i for i in range(mol.GetNumAtoms())
                                 if i not in blacklist_ids])
        # we should use GetPositions() here, but it often leads to segfault (RDKit)
        conf = mol.GetConformer()
        protein_coords = np.array([conf.GetAtomPosition(i)
                                  for i in protein_amap.tolist()])

        # Pocket selection based on cutoff
        mask = (cdist(protein_coords, ligand_coords) <= cutoff).any(axis=1)
        # IDs of atoms within cutoff

        pocket_amap = []
        try:
            pocket_amap = set(protein_amap[np.where(mask)[0]].tolist())
            ligand_residue_keys = ligand_residue_dict[ligand_key]

            pocket_amap.union(set(list(chain(*[protein_residues[ligand_residue_keys[i]]
                                   for i in range(len(ligand_residue_keys))]))))
            pocket_amap = list(pocket_amap)
        except:
            pocket_amap = protein_amap[np.where(mask)[0]].tolist()



        # Expand pocket's residues
        pocket_residues = {}
        for res_id in protein_residues.keys():
            if any(1 for res_aix in protein_residues[res_id]
                    if res_aix in pocket_amap):
                pocket_residues[res_id] = protein_residues[res_id]
        pocket_amap = list(chain(*pocket_residues.values()))

        # Create pocket mol, pocket_amap needs to be mapped to mol Idxs
        pocket = AtomListToSubMol(mol, pocket_amap, includeConformer=True)


        try:
            Chem.SanitizeMol(pocket)
        except:
            print("Could Not Sanitize Pocket for Ligand: ", ligand_key[1])
            continue

        AllChem.ComputeGasteigerCharges(pocket)

        return [pocket, ligand_key[0], ligand_key[1], ligand_key[2], pdb_file]
