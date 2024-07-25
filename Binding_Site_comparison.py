import alphaspace2 as al
import mdtraj
import numpy as np
from scipy.spatial.distance import cdist
from scipy.cluster.hierarchy import fcluster, linkage
from alphaspace2.functions import _binCluster, _group
from alphaspace2.Cluster import _DPocket

import features  ## python module with pseudomolecular features
import glob
from collections import defaultdict
import matplotlib.pyplot as plt

protease_list = ['1c70','1hvi','1hvj','1izh','1pro','1siv','2i0a','2psv','2q5k','3lzu']
pka_list = ['1sve','1xh5','2c1a','2erz','2f7e','2f7x','2gfc','2jds','2oh0','2vo7']
er_list = ['2b1z','2p15','2pog','2q70','2yja','3uud','4mg8','4pps','4pxm','4tv1']

lig = mdtraj.load('Binding_Site_comparison/protease/representative_ligand.pdb')
protease_data_ss = {}
protease_data_prot = {}
for pdb_id in protease_list:
    prot = mdtraj.load('Binding_Site_comparison/protease/protein_' + pdb_id + '.pdb')
    protease_data_prot[pdb_id] = prot
    ss_prot = al.Snapshot()
    ss_prot.run(prot, lig)
    protease_data_ss[pdb_id] = ss_prot

lig = mdtraj.load('Binding_Site_comparison/pka/representative_ligand.pdb')
pka_data_ss = {}
pka_data_prot = {}
for pdb_id in pka_list:
    prot = mdtraj.load('Binding_Site_comparison/pka/protein_' + pdb_id + '.pdb')
    pka_data_prot[pdb_id] = prot
    ss_prot = al.Snapshot()
    ss_prot.run(prot, lig)
    pka_data_ss[pdb_id] = ss_prot

lig = mdtraj.load('Binding_Site_comparison/estrogen_receptor/representative_ligand.pdb')
er_data_ss = {}
er_data_prot = {}
for pdb_id in er_list:
    prot = mdtraj.load('Binding_Site_comparison/estrogen_receptor/protein_' + pdb_id + '.pdb')
    er_data_prot[pdb_id] = prot
    ss_prot = al.Snapshot()
    ss_prot.run(prot, lig)
    er_data_ss[pdb_id] = ss_prot

protease_trajectory = al.Trajectory(snapshots=[protease_data_ss[pdb_id] for pdb_id in protease_data_ss.keys()])
protease_trajectory.gen_dpockets(clust_distance=4.7)
dps = sorted([dp for dp in protease_trajectory.dpockets], key=lambda i: sum(i.scores))

protease_contact_pockets = defaultdict(dict)
for dpx, dp in enumerate(dps):
    pockets = list(dp.pockets)
    for px, pdb_id in enumerate(protease_data_ss.keys()):
        if pockets[px].isContact:
            protease_contact_pockets[pdb_id][dpx] = np.array([b.xyz for b in pockets[px].betas])

protease_props_dict = {}
for pdb_id in protease_contact_pockets:
    contact_betas = []
    prot = protease_data_prot[pdb_id]
    for dpx in protease_contact_pockets[pdb_id]:
        contact_betas.extend(protease_contact_pockets[pdb_id][dpx])

    contact_betas = np.array(contact_betas)
    beta_temp_dict = {}
    beta_temp_dict['occluded_asa'] = features._get_pharmacophore_fingerprint(prot, contact_betas)
    beta_temp_dict['usr'] = features._Get_USR_alpha_beta(contact_betas)
    protease_props_dict[pdb_id] = beta_temp_dict

