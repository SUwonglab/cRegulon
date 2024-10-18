#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 15:26:08 2023

@author: Zaine
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def combine_hex_values(d):
    d_items = sorted(d.items())
    tot_weight = sum(d.values())
    red = int(sum([int(k[:2], 16)*v for k, v in d_items])/tot_weight)
    green = int(sum([int(k[2:4], 16)*v for k, v in d_items])/tot_weight)
    blue = int(sum([int(k[4:6], 16)*v for k, v in d_items])/tot_weight)
    zpad = lambda x: x if len(x)==2 else '0' + x
    return zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])

color = ["8AE0BD","9CABD8","E88DC2","FC9362"]

my_dpi=96;w=700;h=500
plt.figure().set_size_inches(w/72,h/72)
f = open('CL_UMAP.txt')
CL = f.readlines();f.close()
del CL[0]
ct = []
for i in range(len(CL)):
    CL[i] = CL[i].split('\t');
    CL[i][1] = float(CL[i][1]);CL[i][2] = float(CL[i][2]);
    plt.scatter(CL[i][1],CL[i][2],c="#"+color[i])
    plt.text(CL[i][1],CL[i][2],CL[i][0],font="Arial",fontsize=8)
    ct.append(CL[i][0]);del CL[i][0]
CL = np.array(CL)

f = open('CL_A.txt')
L  = f.readlines();f.close()
m = L[0].strip('\n').split('\t');del m[0]
del L[0];ct1 = []
for i in range(len(L)):
    L[i] = L[i].split('\t')
    ct1.append(L[i][0]);del L[i][0]
L = np.array(L).astype('float')
Indel = [ct1.index(ct[i]) for i in range(len(ct))]
L = L[Indel,]

### Plot cRegulon-Cell type association on heatmap
annot_kws={'fontstyle':'Arial'}
hm = sns.heatmap(data = L,
                yticklabels=ct, xticklabels=m,
                cmap="coolwarm",
                linewidths=0, linecolor="white", square=False,
                annot_kws= annot_kws) 
plt.show()

### Plot cRegulon-Cell type association on UMAP
LN = L.T.copy()
for i in range(LN.shape[0]):
    LN[i] = LN[i]/np.sum(LN[i])

CR = np.dot(LN,CL)

for i in range(CR.shape[0]):
    plt.scatter(CR[i][0],CR[i][1],marker="D",c='black')
    plt.text(CR[i][0],CR[i][1],m[i],font="Arial")

for i in range(L.shape[0]):
    for j in range(L.shape[1]):
        if L[i][j] >= 0.02:
            plt.plot([CL[i][0],CR[j][0]],[CL[i][1],CR[j][1]],linestyle="-",
                      linewidth=4*((0.8/0.98)*(L[i][j]-1)+1),
                      c="#"+combine_hex_values({color[i]: 1.0, "000000": 1.0}))
            
plt.xlabel('UMAP-1',font="Arial");plt.xticks(font='Arial')
plt.ylabel('UMAP-2',font="Arial");plt.yticks(font='Arial')
#plt.savefig('Atlas_Landscape.pdf')
plt.show()
