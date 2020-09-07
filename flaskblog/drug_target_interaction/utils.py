import os
import secrets
from PIL import Image
from flask import url_for, current_app
from flask_mail import Message
from flaskblog import mail
from flaskblog.utils.SomeHelperFunctions import *
from pebble import concurrent
from pebble import ProcessPool
from concurrent.futures import TimeoutError
from os import listdir

import pymol as pm
import pyvol as pv
from pyvol.pymol_interface import pymol_pocket
import tempfile

from flaskblog.utils.Protein_Features.QuasiSequenceOrderFeatures import GetQuasiSequenceOrderFeatures
from flaskblog.utils.Protein_Features.PseudoAminoAcidCompositionFeatures import GetPseudoAAC, _mean
from flaskblog.utils.Protein_Features.ConjointTriadFeatures import CalculateConjointTriad
from flaskblog.utils.Protein_Features.CompositionTransitionDistributionFeatures import CalculateCompositionTransitionDistribution
from flaskblog.utils.Protein_Features.AutocorrelationFeatures import CalculateAutocorrelation
from flaskblog.utils.Protein_Features.AminoAcidsCompositionFeatures import CalculateAADipeptideComposition

from flaskblog.utils.Ligands_Protein_Features.ChargesFeatures import GetCharge
from flaskblog.utils.Ligands_Protein_Features.ConnectivityFeatures import GetConnectivity
from flaskblog.utils.Ligands_Protein_Features.ConstitutionFeatures import GetConstitutional
from flaskblog.utils.Ligands_Protein_Features.GearyAutocorrelationFeatures import GetGearyAuto
from flaskblog.utils.Ligands_Protein_Features.MOEFeatures import GetMOE
from flaskblog.utils.Ligands_Protein_Features.MoleculePropertyFeatures import GetMolecularProperty
from flaskblog.utils.Ligands_Protein_Features.MoranAutocorrelationFeatures import GetMoranAuto
from flaskblog.utils.Ligands_Protein_Features.MoreauBrotoAutocorrelationFeatures import GetMoreauBrotoAuto
from flaskblog.utils.Ligands_Protein_Features.TopologicalFeatures import GetTopology

from rdkit.Chem import AllChem
from rdkit import Chem

def get_pockets(file, cutoff=4):
	output = get_residue_list(file,4,['HOH'])
	return output

def get_pocket_with_key(file, key):
    output = get_pocket_with_id(file, key)
    return output

def GetAllBindingSiteFeatures(binding_site_list):

    binding_site = binding_site_list[0]
    residue_id = binding_site_list[1]
    ligand_id = binding_site_list[2]
    chain = binding_site_list[3]
    pdb_id = binding_site_list[4]

    result_pocket = {}
    #result_pocket.update({"Name": ligand_file.split('_')[0]})
    result_pocket.update(GetTopology(binding_site, False))
    result_pocket.update(GetMoreauBrotoAuto(binding_site))
    result_pocket.update(GetMoranAuto(binding_site))
    result_pocket.update(GetMolecularProperty(binding_site))
    result_pocket.update(GetMOE(binding_site))
    result_pocket.update(GetGearyAuto(binding_site))
    result_pocket.update(GetConstitutional(binding_site))
    result_pocket.update(GetConnectivity(binding_site))
    result_pocket.update(GetCharge(binding_site))

#     for key in result_pocket.keys():
#         temp = result_pocket[key]
    result_pocket = {"binding_site_"+k: v for k, v in result_pocket.items()}

    print('pdb')
    return [[ligand_id, chain, residue_id], result_pocket]



def get_sequence(file):
    chains = []
    sequence = []
    while True:
        try:
            temp = file.readline()
            chains.append(temp.split('|')[1].split()[1].split(','))
            sequence.append(file.readline()[:-1])
        except:
            break
    return sequence, chains



def GetAllSequenceFeatures(fasta_file):

    output = []

    #for fasta_file in fasta_files:

    ProteinSequences, chainss = get_sequence(fasta_file)

    for ProteinSequence, chains in zip(ProteinSequences, chainss):

        result = {}

        #result.update({"Name": fasta_file[:4]})

        result.update({"Chain": chains})

        result.update(CalculateAADipeptideComposition(ProteinSequence))

        result.update(GetQuasiSequenceOrderFeatures(ProteinSequence))

        result.update(CalculateAutocorrelation(ProteinSequence))

        result.update(CalculateAADipeptideComposition(ProteinSequence))

        result.update(GetPseudoAAC(ProteinSequence))

        result.update(CalculateCompositionTransitionDistribution(ProteinSequence))

        print('sequence')

        result = {"sequence_"+k: v for k, v in result.items()}

        output.append(result)

    return output



def GetAllLigandFeatures(ligand, ligand_name):
    result_ligand = {}
    result_ligand.update({"Name": ligand_name})
    result_ligand.update(GetTopology(ligand))
    result_ligand.update(GetMoreauBrotoAuto(ligand))
    result_ligand.update(GetMoranAuto(ligand))
    result_ligand.update(GetMolecularProperty(ligand))
    result_ligand.update(GetMOE(ligand))
    result_ligand.update(GetGearyAuto(ligand))
    result_ligand.update(GetConstitutional(ligand))
    result_ligand.update(GetConnectivity(ligand))
    result_ligand.update(GetCharge(ligand))

    result_ligand.update({"Volume": AllChem.ComputeMolVolume(ligand)})
    
    result_ligand = {"ligand_"+k: v for k, v in result_ligand.items()}

    print('ligand')
    return result_ligand


def ConverToLigand(text, typ):
    if typ == 'smiles':
        ligand = MolFromSmiles(text)
    elif typ == 'mol2':
        ligand = MolFromMol2Block(text)
    elif typ == 'mol':
        ligand = MolFromMolBlock(text)
    elif typ == 'sdf':
        import tempfile
        fo = tempfile.NamedTemporaryFile(mode='w+')
        fo.write(text)
        file_name = fo.name
        from rdkit import Chem
        from rdkit.Chem import AllChem
        try:
            fo.seek(0)
            ligand = Chem.SDMolSupplier(file_name)[0]
            fo.close()
        except:
            fo.close()
            raise ValueError("SDF File not read")
    else:
        raise ValueError("Unrecognized file type")

    return ligand



def pymol_pocket_cmdline(protein=None, ligand=None, prot_file=None, lig_file=None, min_rad=1.4, max_rad=3.4, constrain_radii=True, mode="largest", coordinates=None, residue=None, resid=None, lig_excl_rad=None, lig_incl_rad=None, min_volume=200, subdivide=False, max_clusters=None, min_subpocket_rad=1.7, max_subpocket_rad=3.4, min_subpocket_surf_rad=1.0, radial_sampling=0.1, inclusion_radius_buffer=1.0, min_cluster_size=50, project_dir=None, output_dir=None, prefix=None, logger_stream_level="INFO", logger_file_level="DEBUG", protein_only=False, display_mode="solid", alpha=1.0, palette=None):

    opts = {
        "protein": protein,
        "ligand": ligand,
        "prot_file": prot_file,
        "lig_file": lig_file,
        "min_rad": min_rad,
        "max_rad": max_rad,
        "constrain_radii": constrain_radii,
        "mode": mode,
        "residue": residue,
        "resid": resid,
        "coordinates": coordinates,
        "lig_excl_rad": lig_excl_rad,
        "lig_incl_rad": lig_incl_rad,
        "min_volume": min_volume,
        "subdivide": subdivide,
        "max_clusters": max_clusters,
        "min_subpocket_rad": min_subpocket_rad,
        "max_subpocket_rad": max_subpocket_rad,
        "min_subpocket_surf_rad": min_subpocket_surf_rad,
        "radial_sampling": radial_sampling,
        "inclusion_radius_buffer": inclusion_radius_buffer,
        "min_cluster_size": min_cluster_size,
        "project_dir": project_dir,
        "output_dir": output_dir,
        "prefix": prefix,
        "logger_stream_level": logger_stream_level,
        "logger_file_level": logger_file_level,
        "protein_only": protein_only,
        "display_mode": display_mode,
        "alpha": alpha,
        "palette": palette
    }

    return opts


def calcVolume(ids):
    volume = 0
    with tempfile.TemporaryDirectory() as tmpdirname:
        try:
            pm.cmd.fetch(ids[0], path=tmpdirname)
            command = pymol_pocket_cmdline(protein=ids[0], mode='specific',
                                           residue=ids[2] + str(ids[3]),
                                           output_dir=str(ids),
                                           project_dir = tmpdirname)
            volume_returned = pymol_pocket(**command)

            volume = volume_returned[0][0].mesh.volume
        except:
            return (ids, None)
    return (ids, volume)
