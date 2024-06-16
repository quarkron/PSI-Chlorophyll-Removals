from __future__ import division
import numpy as np
import chimera
import os
import csv
from StructBio.Scenographics.solids import Solids
from StructBio.Utilities.miscFunctions import colorByName


energy_dict = {
    'CLA1011.A': 13201,
    'CLA1012.B': 13787,
    'CLA1013.A': 13744,
    'CLA1021.B': 13964,
    'CLA1022.A': 14105,
    'CLA1023.B': 13863,
    'CLA1101.A': 14133,
    'CLA1102.A': 13936,
    'CLA1103.A': 13603,
    'CLA1104.A': 14006,
    'CLA1105.A': 14136,
    'CLA1106.A': 13970,
    'CLA1107.A': 14129,
    'CLA1108.A': 14247,
    'CLA1109.A': 14644,
    'CLA1110.A': 14257,
    'CLA1111.A': 14050,
    'CLA1112.A': 14170,
    'CLA1113.A': 13847,
    'CLA1114.A': 13565,
    'CLA1115.A': 13976,
    'CLA1116.A': 14315,
    'CLA1117.A': 13745,
    'CLA1118.A': 14028,
    'CLA1119.A': 13963,
    'CLA1120.A': 13866,
    'CLA1121.A': 14315,
    'CLA1122.A': 13911,
    'CLA1123.A': 13803,
    'CLA1124.A': 13692,
    'CLA1125.A': 13981,
    'CLA1126.A': 13832,
    'CLA1127.A': 13978,
    'CLA1128.A': 14513,
    'CLA1129.A': 14695,
    'CLA1130.A': 13909,
    'CLA1131.A': 13640,
    'CLA1132.A': 14229,
    'CLA1133.A': 13847,
    'CLA1134.A': 13973,
    'CLA1135.A': 13687,
    'CLA1136.A': 13835,
    'CLA1137.A': 13999,
    'CLA1138.A': 14125,
    'CLA1139.A': 13644,
    'CLA1140.A': 13911,
    'CLA1201.B': 13879,
    'CLA1202.B': 13658,
    'CLA1203.B': 13875,
    'CLA1204.B': 14328,
    'CLA1205.B': 13721,
    'CLA1206.B': 13801,
    'CLA1207.B': 13811,
    'CLA1208.B': 14237,
    'CLA1209.B': 13953,
    'CLA1210.B': 14271,
    'CLA1211.B': 13786,
    'CLA1212.B': 14186,
    'CLA1213.B': 13982,
    'CLA1214.B': 14216,
    'CLA1215.B': 13688,
    'CLA1216.B': 13871,
    'CLA1217.B': 13990,
    'CLA1218.B': 13872,
    'CLA1219.B': 13828,
    'CLA1220.B': 13774,
    'CLA1221.B': 14411,
    'CLA1222.B': 13796,
    'CLA1223.B': 14017,
    'CLA1224.B': 13648,
    'CLA1225.B': 13880,
    'CLA1226.B': 13302,
    'CLA1227.B': 14826,
    'CLA1228.B': 14324,
    'CLA1229.B': 13783,
    'CLA1230.B': 14214,
    'CLA1231.B': 14047,
    'CLA1232.B': 14123,
    'CLA1233.B': 14008,
    'CLA1234.B': 14139,
    'CLA1235.B': 13523,
    'CLA1236.B': 14017,
    'CLA1237.A': 13791,
    'CLA1238.B': 13812,
    'CLA1239.B': 13800,
    'CLA1501.L': 13605,
    'CLA1502.L': 13832,
    'CLA1503.L': 13787,
    'CLA1601.M': 14201,
    'CLA1401.K': 14026,
    'CLA1402.A': 14131,
    'CLA1301.F': 14100,
    'CLA1302.J': 14190,
    'CLA1303.J': 13990,
    'CLA1701.X': 14008,
    'CLA1801.A': 14025
}

def extract(p):
    cla = []
    for r in p.residues:
        if r.type == 'CLA':
            cla.append(r)
    N = len(cla)
    return cla

def linearScaler(min_xval, max_xval, min_yval, max_yval, xval):
    yval = ( (max_yval - min_yval)*xval + (max_xval*min_yval - min_xval*max_yval) )/(max_xval - min_xval)
    return yval

def center(chl):
    mg = chl.findAtom('MG')
    return np.array(mg.coord())

def t_dip(chl):
    nd = np.array(chl.findAtom('ND').coord())
    nb = np.array(chl.findAtom('NB').coord())
    N = np.linalg.norm(nd-nb)
    return (nd-nb)/N 

def coupling(node1, node2):
    C=116000
    r_ij = center(node1) - center(node2)
    R_ij = np.linalg.norm(r_ij)
    di, dj = t_dip(node1), t_dip(node2)
    return C*(np.inner(di, dj))/((R_ij)**3) - (3*np.inner(r_ij, di)*np.inner(r_ij, dj))/((R_ij)**5)

def overlap_integral(Ei, Ej, v=240, S=160):
    return (1/(2*v*np.sqrt(np.pi)))*np.exp((-S**2 - Ei**2 - Ej**2 +2*S*Ei + 2*Ei*Ej - 2*S*Ej)/(4*v**2))

def transfer_rate(i,j, Ei, Ej):
    return coupling(i,j)**2*overlap_integral(Ei, Ej)*(2*np.pi)**2*299792458*100*1e-12

pdbID = raw_input("Type in four character PDB Identifier: \n")

p = chimera.openModels.open(pdbID, type = "PDB")[0]

cla_array = extract(p)

N = len(cla_array)

centers = []
key_chl = set()

rate_solids = Solids()
r=0.75
cutoff=0.5
for i in [0,1,2,47,48,42,46,86]:
    for j in [0,1,2,47,48,42,46,86]:
        if i > j:
            chl_i = cla_array[i].type + str(cla_array[i].id.position) + '.' + cla_array[i].id.chainId
            chl_j =  cla_array[j].type + str(cla_array[j].id.position) + '.' + cla_array[j].id.chainId
            rate_ij = transfer_rate(cla_array[i],cla_array[j], energy_dict[chl_i], energy_dict[chl_j])
            rate_ji = transfer_rate(cla_array[j],cla_array[i], energy_dict[chl_j], energy_dict[chl_i])
            cen_i = center(cla_array[i])
            cen_j = center(cla_array[j])
            R = np.linalg.norm(cen_i-cen_j)
            centers.append(cen_i)
            cen_jprime = ((R-1.5*1)/R)*(cen_j - cen_i) + cen_i
            cen_iprime = ((R-1.5*1)/R)*(cen_i - cen_j) + cen_j
            if rate_ij > rate_ji:
                if rate_ij > cutoff:
                    print(i,j,chl_i, chl_j, rate_ij)
                    rate_solids.addArrow(cen_i, cen_jprime, rate_ij/25, colorByName("deep sky blue"))
            if rate_ij == rate_ji:
                if rate_ij > cutoff:
                    print(i,j, chl_i, chl_j, rate_ij, rate_ji)
                    rate_solids.addArrow(cen_iprime, cen_jprime, rate_ij/25, colorByName("deep sky blue"))
                    rate_solids.addArrow(cen_jprime, cen_iprime, rate_ji/25, colorByName("deep sky blue"))
            if rate_ij < rate_ji:
                if rate_ji > cutoff:
                    print(i,j,chl_i, chl_j, rate_ij)
                    rate_solids.addArrow(cen_j, cen_iprime, rate_ji/25, colorByName("deep sky blue"))
            if rate_ij > 10:
                key_chl.add(i)
                key_chl.add(j)


for i in range(N):
    nodeName = cla_array[i].type + str(cla_array[i].id.position) + '.' + cla_array[i].id.chainId
    cen_i = center(cla_array[i])
    if nodeName in ['CLA1011.A', 'CLA1021.B']:
        rate_solids.addSphere(cen_i, r, colorByName('black'))
    if nodeName in ['CLA1012.B', 'CLA1013.A', 'CLA1022.A', 'CLA1023.B']:
        rate_solids.addSphere(cen_i, r, colorByName('magenta'))
    if i in [42,86]:
        rate_solids.addSphere(center(cla_array[i]), r, colorByName('green'))






rate_solids.display()
# print(key_chl)
# runCommand("set bg_color white")



            



