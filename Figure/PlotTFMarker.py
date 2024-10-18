import numpy as np 
import seaborn as sn 
import matplotlib.pyplot as plt

f = open('CL_X.txt')
X = f.readlines();f.close()
Name = X[0].strip('\n').split('\t');del Name[0]
del X[0]
TF = []
for i in range(len(X)):
    X[i] = X[i].split('\t')
    TF.append(X[i][0]);del X[i][0]

X = np.array(X).astype('float')

Marker = ["GATA1","GATA2","TAL1","JUND","JUN",
          "FOSB","FOS","FOSL2",
          "NANOG","POU5F1",
          "FLI1","IRF8",
          "E2F5","E2F2","E2F6","E2F3","ESRRB",
          "KLF9","KLF4","KLF5","KLF7",
          "GLI2","GLI1","GLI3","ZIC3","ZIC2","E2F1"]
XS = []
for i in range(len(Marker)):
    indel = TF.index(Marker[i])
    XS.append(X[indel])
XS = np.array(XS)

annot_kws={'fontstyle':'Arial'}
hm = sn.heatmap(data = XS.T,
                yticklabels=Name, xticklabels=Marker,
                cmap="coolwarm",
                linewidths=0, linecolor="white", square=False,
                annot_kws= annot_kws) 
plt.show()
