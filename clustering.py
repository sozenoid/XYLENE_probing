#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 00:07:23 2020

@author: macenrola
"""

def ClusterFps(fps,cutoff=0.2):


    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    lfps = len(fps[0])
    for i in range(1,nfps):
        #sims = DataStructs.BulkKulczynskiSimilarity(fps[i],fps[:i])
        #sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        #sims = DataStructs.BulkBraunBlanquetSimilarity(fps[i],fps[:i])
        sims = DataStructs.BulkDiceSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs


if __name__ == "__main__":
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    import gzip
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina
    
# =============================================================================
#     ms = [x for x in Chem.ForwardSDMolSupplier("/home/macenrola/Documents/ML/ChemTS/ledock_ligand_design_with_sa_scorer_adaincl_good_save/RES/sample2.sdf") if x is not None]
#     fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in ms]
#     #fps = [AllChem.GetMACCSKeysFingerprint(x) for x in ms]
#     clusters=ClusterFps(fps,cutoff=0.4)
#     print(len(clusters))
#     toshow=[]
#     for el in clusters[4]:
#         toshow.append(ms[el])
#      
#     img=Draw.MolsToGridImage(toshow)
#     img.show()
# =============================================================================
    ms = [x for x in Chem.SmilesMolSupplier("/home/macenrola/Documents/ML/ChemTS/ledock_ligand_design_with_sa_scorer_adaincl_good_save/RES/MCTS/1000-BEST-MCTS.smi", delimiter='\t') if x is not None]
    [AllChem.Compute2DCoords(x) for x in ms]
    lgd = [x.GetProp("_Name") for x in ms]
    print(len(lgd))
    img= Draw.MolsToGridImage(ms[500:1000], legends=lgd[500:1000], molsPerRow=10, subImgSize=(500,500), maxMols=500)
    img.show()
