import chimera
import csv
import numpy as np

prot = '1JB0' #1JB0 PDB ID for PSI
p = chimera.openModels.open(prot, type='PDB')[0]

def t_dip(chl):
    nd = np.array(chl.findAtom('ND').coord())
    nb = np.array(chl.findAtom('NB').coord())
    N = np.linalg.norm(nd-nb)
    return (nd-nb)/N

def center(chl):
    mg = chl.findAtom('MG')
    return np.array(mg.coord())
             
with open('chls_1JB0.csv', 'wb') as file:
    writer = csv.writer(file)

    writer.writerow(["Index", "ID", "Chl type", "Chain ID", 
                  "Center0", "Center1", "Center2", "Tdip_vec0", "Tdip_vec1", "Tdip_vec2"])
    i = 0
    for r in p.residues:
        if r.type in ['CLA', 'CHL']:
            r_tdipVec = t_dip(r)
            r_ID = r.type + str(r.id.position) + '.' + r.id.chainId
            r_type = r.type
            r_chainID= r.id.chainId
            r_center = center(r)
            writer.writerow([i, r_ID, r_type, r_chainID, r_center[0], r_center[1],
                             r_center[2], r_tdipVec[0], r_tdipVec[1], r_tdipVec[2]])
            i += 1
