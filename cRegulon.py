import os
import gzip
import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import argparse
import sys
from ismember import ismember
import scipy.sparse as sparse
import scipy.io as scio
from numpy.matlib import repmat
import numpy_groupies as npg
from collections import Counter
import mpmath
import subprocess
from multiprocessing import Process
import pybedtools
import shutil
from scipy import stats
from sklearn.decomposition import NMF
from numpy import linalg as LA
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

def Read10X_ATAC(path):
	try:
		with gzip.open(path+'/peaks.bed.gz', 'rt') as f:
			pk = f.readlines();f.close()
			pk = [pk[i].strip('\n').replace('\t','_') for i in range(len(pk))]
	except OSError:
		try:
			with open(path+'/peaks.bed', 'rt') as f:
				pk = f.readlines();f.close()
				pk = [pk[i].strip('\n').replace('\t','_') for i in range(len(pk))]
		except OSError:
			print("No correct peak file!")
	try:
		with gzip.open(path+'/barcodes.tsv.gz', 'rt') as f:
			bc = f.readlines();f.close()
			bc = [bc[i].strip('\n') for i in range(len(bc))]
	except OSError:
		try:
			with open(path+'/barcodes.tsv', 'rt') as f:
				bc = f.readlines();f.close()
				bc = [bc[i].strip('\n') for i in range(len(bc))]
		except OSError:
			print("No correct barcode file!")
	try:
		with gzip.open(path+'/matrix.mtx.gz', 'rt') as f:
			count = np.loadtxt(f, skiprows=2).astype('int')
	except OSError:
		try:
			with open(path+'/matrix.mtx', 'rt') as f:
				count = np.loadtxt(f, skiprows=2).astype('int')
		except OSError:
			print("No correct count file!")

	if count[0][1] != len(bc) or count[0][0] != len(pk):
		print("Error!!!")
	else:
		count = np.delete(count, 0, axis=0).T
		count = csr_matrix((count[2], (count[0]-1, count[1]-1)), shape=(len(pk), len(bc)))
		return pk, bc, count

def Read10X_RNA(path):
	try:
		with gzip.open(path+'/features.tsv.gz', 'rt') as f:
			gs = f.readlines();f.close()
			gs = [gs[i].strip('\n').split('\t')[1] for i in range(len(gs))]
	except OSError:
		try:
			with open(path+'/features.tsv', 'rt') as f:
				gs = f.readlines();f.close()
				gs = [gs[i].strip('\n').split('\t')[1] for i in range(len(gs))]
		except OSError:
			print("No correct peak file!")
	try:
		with gzip.open(path+'/barcodes.tsv.gz', 'rt') as f:
			bc = f.readlines();f.close()
			bc = [bc[i].strip('\n') for i in range(len(bc))]
	except OSError:
		try:
			with open(path+'/barcodes.tsv', 'rt') as f:
				bc = f.readlines();f.close()
				bc = [bc[i].strip('\n') for i in range(len(bc))]
		except OSError:
			print("No correct barcode file!")
	try:
		with gzip.open(path+'/matrix.mtx.gz', 'rt') as f:
			count = np.loadtxt(f, skiprows=2).astype('int')
	except OSError:
		try:
			with open(path+'/matrix.mtx', 'rt') as f:
				count = np.loadtxt(f, skiprows=2).astype('int')
		except OSError:
			print("No correct count file!")

	if count[0][1] != len(bc) or count[0][0] != len(gs):
		print("Error!!!")
	else:
		count = np.delete(count, 0, axis=0).T
		count = csr_matrix((count[2], (count[0]-1, count[1]-1)), shape=(len(gs), len(bc)))
		return gs, bc, count

def ReadTXT(path):
	df = pd.read_csv(path, index_col=0, sep='\t')
	row_names = df.index.tolist()
	column_names = df.columns.tolist()
	matrix = csr_matrix(df.values)
	return row_names, column_names, matrix

def PseudoBulk(name,rna,rna_meta,atac,atac_meta):
	folder = "./PseudoBulk/";os.makedirs("./PseudoBulk/", exist_ok=True)
	f = open(rna_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(tmeta));ct.sort()
	print("We are processing scRNA-seq data with "+str(len(ct))+" cell clusters...")
	if os.path.isfile(rna):
		gene,cell,data = ReadTXT(rna)
	else:
		gene,cell,data = Read10X_RNA(rna)
	gc = open(folder+name+"_"+'CellType.txt','w')
	for i in range(len(ct)):
		gc.write(name+"_"+ct[i]+'\n')
		Indel = [cell.index(cmeta[indel]) for indel,x in enumerate(tmeta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_RNA.txt','w')
		for j in range(len(gene)):
			g.write(gene[j]+'\t'+str(datac[j,0])+'\n')
		g.close()
	gc.close()

	f = open(atac_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]
	ct = list(set(tmeta));ct.sort()
	print("We are processing scATAC-seq data with "+str(len(ct))+" cell clusters...")
	if os.path.isfile(atac):
		peak,cell,data = ReadTXT(atac)
	else:
		peak,cell,data = Read10X_ATAC(atac)
	for i in range(len(ct)):
		Indel = [cell.index(cmeta[indel]) for indel,x in enumerate(tmeta) if x==ct[i]]
		datac = data[:,Indel]
		datac = np.sum(datac,1)/np.sum(datac)*1000000
		g = open(folder+name+"_"+ct[i]+'_PS_ATAC.txt','w')
		for j in range(len(peak)):
			if datac[j,0] >= 2 and "GL" not in peak[j] and "JH" not in peak[j]:
				g.write(peak[j].replace(':','_').replace('-','_')+'\t'+str(datac[j,0])+'\n')

def TFC(x,y):
	eps = 1e-4
	if y > 0:
		fc = x/y-0.5
	else:
		fc = x/eps-0.5
	if fc > 1:
		fc = 1
	if fc < 0:
		fc = 0
	return fc

def TFAS(name,rna,rna_meta,genome):
	f = open('./PseudoBulk/'+name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	f = open("./Data/AllTFs_"+genome+".txt")
	TF = f.readlines();f.close()
	TF = [TF[i].strip('\n') for i in range(len(TF))]

	f = open(rna_meta)
	meta = f.readlines();f.close()
	cmeta = [meta[i].strip('\n').split('\t')[0] for i in range(len(meta))]
	tmeta = [name+"_"+meta[i].strip('\n').split('\t')[1] for i in range(len(meta))]

	if os.path.isfile(rna):
		gene,cell,data = ReadTXT(rna)
	else:
		gene,cell,data = Read10X_RNA(rna)
	cindex = [cell.index(cmeta[i]) for i in range(len(cmeta))];data = data[:,cindex]
	Indel_ct = [[] for c in range(len(ct))]
	for i in range(len(cmeta)):
		Indel_ct[ct.index(tmeta[i])].append(i)

	AllTFAS = [];TFs = []
	for i in range(len(TF)):
		if TF[i] in gene:
			alltfas = [];TFs.append(TF[i])
			indel = gene.index(TF[i])
			for c in range(len(ct)):
				e1 = np.mean(data[indel][:,Indel_ct[c]])
				e2 = np.mean(np.delete(data[indel].toarray(),Indel_ct[c]))
				alltfas.append(TFC(e1,e2))
			AllTFAS.append(alltfas)
	g = open("./PseudoBulk/"+name+"_TFES.txt",'w')
	g.write(name+"\t"+"\t".join(ct)+'\n')
	for i in range(len(TFs)):
		g.write(TFs[i])
		for c in range(len(ct)):
			g.write("\t"+str(AllTFAS[i][c]))
		g.write("\n")
	g.close()

def mode_prep(args):
	PseudoBulk(args.name,args.rna,args.rna_meta,args.atac,args.atac_meta)
	TFAS(args.name,args.rna,args.rna_meta,args.species)


def mfbs(MCFile, TFName, Element_name, motifName, motifWeight, Match2):
	MC = pd.read_csv(MCFile, sep='\t', header=None)
	Mf3 = MC.iloc[:,2]
	Md1, Mf1 = ismember(MC.iloc[:,0], Element_name)
	Md2, Mf2 = ismember(MC.iloc[:,1], pd.DataFrame(motifName))
	Mf1_2 = np.zeros(shape=(MC.shape[0],1))
	Mf2_2 = np.zeros(shape=(MC.shape[0],1))
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mt1 = np.setdiff1d(np.arange(len(motifName)), np.unique(Mf2))
	Mf2 = np.concatenate((np.squeeze(Mf2_2[(Md1 & Md2)]), Mt1 ))
	Mf1 = np.concatenate((np.squeeze(Mf1_2[(Md1 & Md2)]), np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3[np.where(Md1 & Md2)[0]], np.zeros(len(Mt1))))
	Mt1 = np.setdiff1d(np.arange(len(Element_name)), np.unique(Mf1))
	Mf1 = np.concatenate((Mf1, Mt1))
	Mf2 = np.concatenate((Mf2, np.ones(len(Mt1))))
	Mf3 = np.concatenate((Mf3, np.zeros(len(Mt1))))
	MMotif_binding = sparse.csr_matrix((Mf3, (Mf2, Mf1)), shape=(len(motifName), len(Element_name)))
	MMotif_binding = sparse.diags(np.squeeze(np.asarray(1/(motifWeight + 0.1)))) * MMotif_binding
	MMotif_binding = np.log(MMotif_binding.todense()+1)
	MTF_binding = np.zeros((len(TFName), len(Element_name)))
	Mf1_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Mf2_2 = np.zeros(shape=(Match2.shape[0], 1))-1
	Md1, Mf1 = ismember(pd.DataFrame(Match2[:, 0]), motifName)
	Md2, Mf2 = ismember(pd.DataFrame(Match2[:, 1]), TFName)
	Mf1_2[np.where(Md1 == True)[0], 0] = Mf1
	Mf2_2[np.where(Md2 == True)[0], 0] = Mf2
	Mf1 = np.squeeze(Mf1_2)
	Mf2 = np.squeeze(Mf2_2)
	for i in range(len(TFName)):
		Ma = np.where(Mf2 == i)[0]
		if len(Ma) > 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :].max(axis=0)
		elif len(Ma) == 1:
			MTF_binding[i, :] = MMotif_binding[Mf1[Ma].astype(np.int64), :]
		else:
			MTF_binding[i, :] = np.zeros((len(Element_name),))
	MTF_binding = sparse.csr_matrix(MTF_binding)
	return MTF_binding

def mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2):
	TFB = sparse.csr_matrix(([], ([], [])), shape=(len(TFName), len(Element_name)))
	for i in range(N):
		TFB1 = mfbs(".MotifTarget"+str(i+1)+".txt",TFName, Element_name, motifName, motifWeight, Match2)
		TFB = TFB + TFB1
	return TFB

def mf_collect(genome):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	filename = 'knownResults.txt'
	fileID = open(filename)
	C = np.loadtxt(fileID, dtype=str, skiprows=1)
	fileID.close()

	for i in range(C.shape[0]):
		C[i, 6] = float(C[i, 6].strip('%'))
		C[i, 8] = float(C[i, 8].strip('%'))
	Score = np.zeros(shape=(C.shape[0], 3))
	Score[:, 0] = np.array([float(-mpmath.log10(mpmath.mpf(C[i, 2]))) for i in range(C.shape[0])])
	Score[:, 1] = (C[:, 6].astype(float) + 0.1) / (C[:, 8].astype(float) + 0.1)
	Score[Score[:, 0] > 100, 0] = 100
	Score[:, 2] = np.sqrt(Score[:, 0] * Score[:, 1])

	f = np.argsort(Score[:, 2])[::-1]
	Score = Score[f, :]
	Name = C[f, 0]

	filename = 'knownResults_rank.txt'
	fid = open(filename, 'w')
	fid.write('Motif\t-log10(p)\tFoldChange\tScore\n')
	for i in range(Name.shape[0]):
		fid.write('{}\t{}\t{}\t{}\n'.format(Name[i], Score[i, 0], Score[i, 1], Score[i, 2]))
	fid.close()

	Motifmatch = scio.loadmat('../../../Data/MotifMatch_{}_rmdup.mat'.format(species))

	Match2 = np.empty((Motifmatch['Match2'].shape[0], 2), dtype='U100')
	for i in range(Motifmatch['Match2'].shape[0]):
		Match2[i, 0] = Motifmatch['Match2'][i, 0].item()
		Match2[i, 1] = Motifmatch['Match2'][i, 1].item()

	d, f = ismember(Match2[:, 0], Name)
	Score1 = Score[f, :]
	Match2 = Match2[np.where(d == True), :]
	Match2 = np.squeeze(Match2)
	TF, ic = np.unique(Match2[:, 1], return_inverse=True)
	TFScore = npg.aggregate(ic, Score1[:, 2], func='max')
	f = np.argsort(TFScore)[::-1]
	Results = np.column_stack((TF[f], TFScore[f]))
	filename = 'knownResults_TFrank.txt'
	fid = open(filename, 'w')
	for i in range(Results.shape[0]):
		fid.write('{}\t{}\n'.format(Results[i, 0], Results[i, 1]))
	fid.close()
    

def Split(N):
	f = open("region.txt")
	pn = f.readlines();f.close()
	NN = len(pn);NM = int(NN/N)+1
	for i in range(N):
		g = open(".region"+str(i+1)+'.txt','w')
		for j in range(len(pn)):
			if int(j/NM) == i:
				g.write(pn[j])
		g.close()

def Homer_prior(peak,genome,output):
	command="mkdir "+output+";findMotifsGenome.pl "+peak+" "+genome+" ./"+output+" -size given -find ../../Data/all_motif_rmdup -preparsedDir ../../Data/Homer/ > ./"+output+"/"+output+".bed;cat ./"+output+"/"+output+".bed | awk 'NR>1'|cut -f 1,4,6 > "+output+".txt;rm -rf "+output
	logfile=open(output+".log",'w')
	ph=subprocess.Popen(command,shell=True,stderr=logfile)
	return_code=ph.wait()

def Homer(peak,genome,output):
	command="mkdir "+output+";findMotifsGenome.pl "+peak+" "+genome+" ./"+output+" -size given -find ../../Data/all_motif_rmdup > ./" + output + "/" + output + ".bed;cat ./" + output + "/" + output + ".bed | awk 'NR>1'|cut -f 1,4,6 > " + output + ".txt;rm -rf " + output
	logfile=open(output+".log",'w')
	ph = subprocess.Popen(command, shell=True,stderr=logfile)
	return_code = ph.wait()

def MotifFind(genome, num_processes, prior):
	"""
	Args:
		genome (str): 
		num_processes (int): 线程数
		prior (int): 是否有已经扫好的homer motif文件
	"""
	processes = []
	Split(num_processes)
	if prior == 1:
		for p in range(num_processes):
			process = Process(target=Homer_prior, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	else:
		for p in range(num_processes):
			process = Process(target=Homer, args=(".region"+str(p+1)+'.txt',str(genome),".MotifTarget"+str(p+1)))
			processes.append(process)
	for process in processes:
		process.start()
	for process in processes:
		process.join()
        
def GRN(name, celltype, genome,num_processes):
	if "hg" in genome:
		species = "human"
	elif "mm" in genome:
		species = "mouse"
	C = pd.read_csv('openness2.bed', sep='\t', header=None)
	Element_name = C.iloc[:, 0]
	Opn = C.iloc[:, 1]
	Opn_median = C.iloc[:, 2]

	# 读取mat文件
	# Match2, motifName, motifWeight
	MotifMatch_mouse_rmdup = scio.loadmat('../../Data/MotifMatch_{}_rmdup.mat'.format(species))
	# Exp_median, List, R2, TFExp_median, TFName
	TFTG_corr_mouse = scio.loadmat('../../Data/TFTG_corr_{}.mat'.format(species))

	Match2 = np.empty((MotifMatch_mouse_rmdup['Match2'].shape[0], 2), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['Match2'].shape[0]):
		Match2[i, 0] = MotifMatch_mouse_rmdup['Match2'][i, 0].item()
		Match2[i, 1] = MotifMatch_mouse_rmdup['Match2'][i, 1].item()

	motifName = np.empty((MotifMatch_mouse_rmdup['motifName'].shape[0], 1), dtype='U100')
	for i in range(MotifMatch_mouse_rmdup['motifName'].shape[0]):
		motifName[i, 0] = MotifMatch_mouse_rmdup['motifName'][i, 0].item()

	motifWeight = MotifMatch_mouse_rmdup['motifWeight']

	TFName = np.empty((TFTG_corr_mouse['TFName'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['TFName'].shape[0]):
		TFName[i, 0] = TFTG_corr_mouse['TFName'][i, 0].item()

	List = np.empty((TFTG_corr_mouse['List'].shape[0], 1), dtype='U100')
	for i in range(TFTG_corr_mouse['List'].shape[0]):
		List[i, 0] = TFTG_corr_mouse['List'][i, 0].item()

	R2 = TFTG_corr_mouse['R2']
	Exp_median = TFTG_corr_mouse['Exp_median']
	# ---------------------------

	N = num_processes
	TF_binding = mfbs_c(N,TFName, Element_name, motifName, motifWeight, Match2)

	# gene expr
	C = pd.read_csv("../../PseudoBulk/{}_PS_RNA.txt".format(celltype), sep='\t', header=None)
	Symbol = C.iloc[:, 0]
	G = C.iloc[:, 1]

	alhfa = 0.5
	Opn_median = np.log2(1 + Opn_median)
	Opn1 = np.log2(1 + Opn)
	Opn = Opn1 * (Opn1 / (Opn_median + 0.5))

	geneName = np.intersect1d(List, Symbol)
	d, f = ismember(geneName, List)
	R2 = R2[:, f]
	Exp_median = Exp_median[f]

	d, f = ismember(geneName, Symbol)
	G = G[f]
	d1 = sorted(G)
	f1 = np.argsort(np.array(G))
	d2 = sorted(Exp_median)
	f2 = np.argsort(np.squeeze(Exp_median))
	G1 = np.empty(len(G))
	for i in range(len(G)):
		G1[f1[i]] = d2[i][0]
	G = np.multiply((np.power(G1, alhfa)), (G1 / (np.squeeze(Exp_median) + 0.5)))

	d, f = ismember(TFName, geneName)
	TFName = TFName[np.where(d == True)]
	TF_binding = TF_binding[np.where(d == True)[0], :]
	TFExp = G[f]
	R2 = R2[np.where(d == True)[0], :]

	C = pd.read_csv('./Enrichment/knownResults_TFrank.txt', header=None, sep='\t')
	d, f = ismember(TFName, C.iloc[:, 0])
	TF_motif = np.zeros(len(TFName))
	TF_motif[np.where(d == True)[0]] = C.iloc[:, 1][f]
	TFExp = TFExp * TF_motif

	C = pd.read_csv('peak_gene_100k_corr.bed', header=None, sep='\t')
	d, f = ismember(C.iloc[:, 0], Element_name)
	d1, f1 = ismember(C.iloc[:, 1], geneName)
	f_2 = np.zeros(shape=(C.shape[0], 1))
	f1_2 = np.zeros(shape=(C.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1

	f2, ia, ic = np.unique(np.hstack((f_2[(d & d1)], f1_2[(d & d1)])), axis=0, return_index=True, return_inverse=True)
	c3 = npg.aggregate(ic, C.iloc[d & d1, 2], func='min')
	c4 = npg.aggregate(ic, C.iloc[d & d1, 3], func='min')
	c4[np.where(c4 < 0.2)] = 0
	d0 = 500000
	c = np.exp(-1 * c3 / d0) * c4
	Opn[np.isnan(Opn)] = 0
	H1 = sparse.csr_matrix((c, (f2[:, 1], f2[:, 0])), shape=(len(geneName), len(Element_name)))
	TFO = np.multiply(TF_binding.todense(), np.tile(Opn.T, (np.size(TF_binding, 0), 1)))
	H1Tdense = H1.todense().T
	BOH = np.matmul(TFO, H1Tdense)
	Score = np.multiply(np.multiply((np.dot(TFExp.reshape(-1, 1), G.reshape(-1, 1).T)), (2 ** np.abs(R2))), BOH)
	Score[np.isnan(Score)] = 0
	np.savetxt('TFTG_regulationScore.txt', Score, delimiter='\t')
	np.savetxt('TFName.txt', TFName, fmt='%s', delimiter='\n')
	np.savetxt('TGName.txt', geneName, fmt='%s', delimiter='\n')

	TFTGControl = scio.loadmat('../../Data/TFTG_{}_nagetriveControl.mat'.format(species))
	Back_net = TFTGControl['Back_net']
	for i in range(Back_net.shape[0]):
		for j in range(Back_net.shape[1]):
			Back_net[i, j] = Back_net[i, j].item()
	d, f = ismember(Back_net[:, 0], TFName)
	d1, f1 = ismember(Back_net[:, 1], geneName)
	f_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f1_2 = np.zeros(shape=(Back_net.shape[0], 1))
	f_2[np.where(d == True)[0], 0] = f
	f1_2[np.where(d1 == True)[0], 0] = f1
	f2 = np.hstack((f_2[(d & d1)], f1_2[(d & d1)]))
	Score_T_1col = Score.T.reshape(-1, 1)
	aa = (f2[:, 1]) * Score.shape[0] + f2[:, 0]
	Back_score = Score_T_1col[aa.astype(np.int64)].squeeze().T

	Cut = np.percentile(np.asarray(Back_score), 99)
	[b, a] = np.where(Score.T > Cut)
	c = np.where(Score_T_1col > Cut)[0]
	c1 = Score_T_1col[c]
	Net = np.column_stack((TFName[a], geneName[b]))
	a1 = np.sort(H1.T.todense(), axis=0)[::-1]
	a2 = np.argsort(H1.T.todense(), axis=0)[::-1]
	a1 = a1[:10, :]
	a2 = a2[:10, :]
	TFTG_RE = [';'.join(Element_name[np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()]) for i in range(len(a))]
	for i in range(len(a)):
		kk = np.asarray(a2[np.where((TFO[a[i], a2[:, b[i]]] > 0) & (a1[:, b[i]] > 0))[0], b[i]]).squeeze()
		if kk.size == 1:
			TFTG_RE[i] = TFTG_RE[i].replace(';', '')

	d = np.sort(np.asarray(c1).squeeze())[::-1]
	f = np.argsort(np.asarray(c1).squeeze())[::-1]
	Net = np.column_stack((Net[f], np.asarray(d).squeeze(), np.asarray(TFTG_RE)[f]))
	filename = '{}_network.txt'.format(celltype)
	with open(filename, 'wt') as fid:
		fid.write('\t'.join(['TF', 'TG', 'Score', 'FDR', 'REs']) + '\n')
		for i in range(Net.shape[0]):
			fid.write('\t'.join([Net[i, 0], Net[i, 1], str(Net[i, 2]), str((np.sum(Back_score > d[i]) + 1) / len(Back_score)),Net[i, 3]]) + '\n')

def network(name,celltype,self_genome, num_processes=20, prior=0):
	
	folder_name = "./Networks/"+celltype
	os.makedirs(folder_name, exist_ok=True)
	os.chdir('./{}'.format(folder_name))
	
	# openness, region 
	df = pd.read_csv("../../PseudoBulk/{}_PS_ATAC.txt".format(celltype), sep = '\t', header=None, names=['col1', 'col2', 'col3'])
	df[['chr', 'start', 'end']] = df['col1'].str.split('_', expand=True)
	df = df[['chr', 'start', 'end', 'col2']]
	df.to_csv('openness1.bed', sep='\t', header=False, index=False)
	df = pd.read_csv('openness1.bed', sep='\t', header=None)
	sorted_df = df.sort_values(by=[0, 1])
	sorted_df.to_csv('openness1.bed', sep='\t', header=False, index=False)
	sorted_df[3] = sorted_df[0] + '_' + sorted_df[1].astype(str) + '_' + sorted_df[2].astype(str)
	sorted_df[[0, 1, 2, 3]].to_csv('region.txt', sep='\t', header=False, index=False)
	sorted_df[[0, 1, 2]].to_csv('region.bed', sep='\t', header=False, index=False)

	print("step1: motif binding...")
	MotifFind(self_genome, num_processes, prior)
	
	print("step2: Opn files...")
	openness1 = pybedtools.BedTool('openness1.bed')
	opn_median = pybedtools.BedTool('../../Data/Opn_median_{}.bed'.format(self_genome))
	result = openness1.intersect(opn_median, wa=True, wb=True, sorted=True).cut([0, 1, 2, 3, 7]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str) + '_' + result.iloc[:, 3].astype(str)
	result = result.iloc[:, [0, 4]]
	max_values = result.groupby('chrom')['score'].max().reset_index()
	max_values2 = max_values.iloc[:, 0].str.rsplit('_', n=1, expand=True)
	max_values2['score'] = max_values.iloc[:, 1]
	max_values2.to_csv("openness2.bed", sep='\t', header=False, index=False)

	os.mkdir("Enrichment")
	data = pd.read_csv("openness2.bed", sep='\t', header=None)
	data['new_col'] = (data[1] + 0.5) / (data[2] + 0.5)
	sorted_data = data.sort_values('new_col', ascending=False)
	filtered_data = sorted_data[sorted_data.iloc[:, 3] - sorted_data.iloc[:, 2] < 2000].head(10000)
	filtered_data[['chr', 'start', 'end']] = filtered_data.iloc[:, 0].str.split('_', expand=True)
	filtered_data[['chr', 'start', 'end']].to_csv("./Enrichment/region.bed", sep='\t', header=False, index=False)

	print("step3: Motif collect...")
	os.chdir("./Enrichment/")
	if prior==1:
	    command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -preparsedDir ../../../Data/Homer/ -p 30".format(self_genome)
	    logfile=open("homer.log",'w')
	    subprocess.run(command, shell=True, stderr=logfile)
	else:
	    command = "findMotifsGenome.pl region.bed {} ./. -size given -mask -nomotif -mknown ../../../Data/all_motif_rmdup -p 30".format(self_genome)
	    logfile=open("homer.log",'w')
	    subprocess.run(command, shell=True, stderr=logfile)
	mf_collect(self_genome)

	print("step4: Prior...")
	os.chdir('../')
	a = pybedtools.BedTool('region.bed')
	b = pybedtools.BedTool('../../Data/RE_gene_corr_{}.bed'.format(self_genome))
	result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
	result = result.iloc[:, [0, 3, 4, 5]]
	result.to_csv("peak_gene_100k_corr.bed", sep='\t', header=False, index=False)

	b = pybedtools.BedTool('../../Data/Enhancer_RE_gene_corr_{}.bed'.format(self_genome))
	result = a.intersect(b, wa=True, wb=True, sorted=True).cut([0, 1, 2, 6, 7, 8]).to_dataframe()
	result.iloc[:, 0] = result.iloc[:, 0].astype(str) + '_' + result.iloc[:, 1].astype(str) + '_' + result.iloc[:,2].astype(str)
	result = result.iloc[:, [0, 3, 4, 5]]
	result.to_csv("peak_gene_100k_corr.bed", mode='a', sep='\t', header=False, index=False)
	
	print("step5: Network...")
	GRN(name,celltype,self_genome,num_processes)
	
	# 清理中间文件
	
	for filename in os.listdir():
	    if filename.startswith('.'):
	        os.remove(filename)
	shutil.rmtree('Enrichment')
	os.remove("region.bed");os.remove("region.txt");os.remove("peak_gene_100k_corr.bed");os.remove("openness1.bed");os.remove("openness2.bed");
	print('GRN for '+celltype+' done!')

def mode_grn(args):
	network(name=args.name,celltype=args.celltype,self_genome=args.genome,num_processes=args.cores, prior=1)


def connection_specificity_index(R):
	R0 = np.log2(R+1)
	R1 = stats.zscore(R0,1);R1[np.isnan(R1)] = 0.0
	R2 = stats.zscore(R0,0);R2[np.isnan(R2)] = 0.0
	R3 = (R1+R2)/2;R3[R3<0] = 0
	P = np.corrcoef(R3);P[np.isnan(P)] = 0.0
	csi = np.zeros((R3.shape[0],R3.shape[0]))
	for i in range(csi.shape[0]):
		for j in range(i,csi.shape[0]):
			n = 0;cutoff = P[i][j] - 0.05
			csi[i][j] = np.sum((P[i]<cutoff)*(P[j]<cutoff))/csi.shape[0]
			csi[j][i] = csi[i][j]
	return R3, csi

def find_intersection(lists):
	intersection_set = set(lists[0])
	for lst in lists[1:]:
		intersection_set &= set(lst)
	return list(intersection_set)

def CSI(name):
	f = open('./PseudoBulk/'+name+'_CellType.txt')
	ct = f.readlines();f.close()
	AllTF = [];AllTG = [];AllTRS = []
	for c in range(len(ct)):
		ct[c] = ct[c].strip('\n')
		AllTRS.append(np.loadtxt('./Networks/'+ct[c]+'/TFTG_regulationScore.txt'))
		f = open("./Networks/"+ct[c]+"/TFName.txt")
		a = f.readlines();f.close();a = [a[i].strip('\n') for i in range(len(a))]
		AllTF.append(a)
		f = open("./Networks/"+ct[c]+"/TGName.txt")
		a = f.readlines();f.close();a = [a[i].strip('\n') for i in range(len(a))]
		AllTG.append(a)

	TF0 = find_intersection(AllTF);TG0 = find_intersection(AllTG);TF0.sort();TG0.sort()
	TFI = [[AllTF[c].index(TF0[i]) for i in range(len(TF0))] for c in range(len(ct))]
	TGI = [[AllTG[c].index(TG0[i]) for i in range(len(TG0))] for c in range(len(ct))]
	AllTRS = [AllTRS[c][TFI[c],:][:,TGI[c]] for c in range(len(ct))]

	TFSTD = np.array([np.std(AllTRS[c],1) for c in range(len(ct))]).max(0)
	TGSTD = np.array([np.std(AllTRS[c],0) for c in range(len(ct))]).max(0)
	TF1 = [];TFI = []
	for i in range(len(TF0)):
		if TFSTD[i] > 0:
			TF1.append(TF0[i])
			TFI.append(i)
	TG1 = [];TGI = []
	os.makedirs("./CSI/", exist_ok=True)
	g = open("./CSI/TGName.txt","w")
	for i in range(len(TG0)):
		if TGSTD[i] > 0:
			TG1.append(TG0[i]);g.write(TG0[i]+"\n")
			TGI.append(i)
	g.close()
	AllTRS = [AllTRS[c][TFI,:][:,TGI] for c in range(len(ct))]

	STD = np.array([np.std(connection_specificity_index(AllTRS[c])[1],1) for c in range(len(ct))]).max(0)
	TF2 = [];TFI = []
	g = open("./CSI/TFName.txt",'w')
	for i in range(len(TF1)):
		if STD[i] > 0:
			TF2.append(TF1[i]);g.write(TF1[i]+'\n')
			TFI.append(i)
	g.close()
	AllTRS = [AllTRS[c][TFI,:] for c in range(len(ct))]
    
	AllCSI = []
	for c in range(len(ct)):
		TRS, CSI = connection_specificity_index(AllTRS[c])
		np.savetxt('./CSI/'+ct[c]+'_CSI.txt',CSI,delimiter='\t',fmt='%1.8f')
		np.savetxt('./CSI/'+ct[c]+'_TRS.txt',TRS,delimiter='\t',fmt='%1.8f')
		AllCSI.append(CSI)
	return AllCSI, TF2, TG1

def ReadTFAS(Name,TF):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	
	f = open('./PseudoBulk/'+Name+'_TFES.txt')
	TFAS0 = f.readlines();f.close()
	del TFAS0[0]
	TF1 = []
	for i in range(len(TFAS0)):
		TFAS0[i] = TFAS0[i].split('\t')
		TF1.append(TFAS0[i][0]);del TFAS0[i][0]
	TFAS0 = np.array(TFAS0)

	AllTFAS = []
	for i in range(len(TF)):
		indel = TF1.index(TF[i])
		AllTFAS.append(TFAS0[indel])
	return np.array(AllTFAS).T.astype('float')

def CSIS(CS,TS,Name):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	CSTS = []
	for c in range(len(CS)):
		TFES = np.dot(TS[c].reshape((TS[c].shape[0],1)),TS[c].reshape((1,TS[c].shape[0])))
		csist = (CS[c]*TFES)**1/3
		CSTS.append(csist)
		np.savetxt('./CSI/'+ct[c]+'_CSIS.txt',csist,delimiter='\t',fmt='%1.8f')
	return CSTS

eps = 0.001
def Norm(X):
	X = X/np.sqrt(np.sum(X**2,0).reshape(1,X.shape[1]))
	X = X/(np.sum(X,1).reshape(X.shape[0],1)+eps)
	return X
def NormL(L):
	L = L/np.sum(L)
	return L
def Loss(C,x0,l0,mu1):
	t1 = 0;t2 = 0
	for c in range(len(C)):
		t1 += pow(LA.norm(C[c] - x0.dot(l0[c]).dot(x0.T), ord='fro'), 2)
		t2 += mu1*np.sum(C[c] * x0.dot(l0[c]).dot(x0.T))
	return t1 - t2

def RunModuleK(C,K):
	print("Identifying", K, "cRegulons...")
	I = len(C)
	C0 = np.zeros(C[0].shape)
	for c in range(I):
		C0 += C[c]
	C0 = C0/I
	err1 = np.zeros(50)
	for i in range(50):
		model = NMF(n_components=K, init='random', random_state=i, solver='cd', max_iter=50)
		X0 = model.fit_transform(C0);X1 = model.components_
		err1[i] = LA.norm(C0 - np.dot(X0, X1), ord='fro')
	model = NMF(n_components=K, init='random', random_state=np.argmin(err1), solver='cd', max_iter=1000)
	X0 = model.fit_transform(C0);X1 = model.components_
	L0 = []
	for c in range(I):
		Li = np.diag([1/K for i in range(K)])
		L0.append(Li)
	L1 = 0;L2 = 0
	for c in range(I):
		L1 += pow(LA.norm(C[c] - X0.dot(L0[c]).dot(X0.T), ord='fro'), 2)
		L2 += np.sum(C[c] * (X0.dot(L0[c]).dot(X0.T)))
	mu1 = L2/L1
	err = 1000;
	X = X0.copy();L = L0.copy();loss1 = Loss(C,X,L,mu1)
	while err > 1e-4:
		XTX = X.T.dot(X)
		XX1 = np.zeros(X.shape);XX2 = np.zeros(X.shape)
		for c in range(I):
			XX1 += C[c].dot(X).dot(L[c])
			XX2 += X.dot(L[c]).dot(XTX).dot(L[c])
		XNext = X * ((1+0.5*mu1)*XX1)/(eps+XX2)
		XNext = Norm(XNext)
		LNext = []
		for c in range(I):
			LNext.append(L[c]*((1+mu1)*X.T.dot(C[c]).dot(X))/(eps+XTX.dot(L[c]).dot(XTX)))
			LNext[c] = NormL(LNext[c])
                
		loss2 = Loss(C,XNext,LNext,mu1)
		err = np.abs(loss2-loss1)
		loss1 = loss2;
		X = XNext.copy();L = LNext.copy()
	return loss1, X, L

def TFPairs(X, TF, Name):
	X = X.T
	cTF = []
	for kk in range(X.shape[0]):
		XK = X[kk]
		XT = [];TFT = []
		for i in range(80):
			indel = np.argmax(XK)
			XT.append(XK[indel]);TFT.append(TF[indel])
			XK[indel] = -1000
		XT = np.array(XT);XT = XT.reshape(XT.shape[0],1)
		XXT = np.dot(XT,XT.T)
		sample = []
		for j in range(XXT.shape[0]):
			for k in range(j+1,XXT.shape[0]):
				sample.append(XXT[j][k])
		par = stats.gamma.fit(sample)
		g = open("./Results/"+Name+'/cRegulon'+str(kk+1)+'_TFModule.txt','w')
		cTF.append([])
		for j in range(XXT.shape[0]):
			for k in range(j+1,XXT.shape[0]):
				pval = stats.gamma.sf(XXT[j][k],par[0],par[1],par[2])
				if pval <= 0.05:
					g.write(TFT[j]+'\t'+TFT[k]+'\t'+str(XXT[j][k])+'\t'+str(pval)+'\n')
					cTF[kk].append(TFT[j]);cTF[kk].append(TFT[k])
		cTF[kk] = list(set(cTF[kk]))
		if len(cTF) == 0:
			XK = X[kk]
			XT = [];TFT = []
			for i in range(80):
				indel = np.argmax(XK)
				XT.append(XK[indel]);TFT.append(TF[indel])
				XK[indel] = -1000
			XT = np.array(XT);XT = XT.reshape(XT.shape[0],1)
			XXT = np.dot(XT,XT.T)
			sample = []
			for j in range(XXT.shape[0]):
				for k in range(j+1,XXT.shape[0]):
					if XXT[j][k] >= 0.05:
						sample.append(XXT[j][k])
			par = stats.gamma.fit(sample)
			g = open("./Results/"+Name+'/cRegulon'+str(kk+1)+'_TFModule.txt','w')
			for j in range(XXT.shape[0]):
				for k in range(j+1,XXT.shape[0]):
					pval = stats.gamma.sf(XXT[j][k],par[0],par[1],par[2])
					if pval <= 0.05:
						g.write(TFT[j]+'\t'+TFT[k]+'\t'+str(XXT[j][k])+'\t'+str(pval)+'\n')
						cTF[kk].append(TFT[j]);cTF[kk].append(TFT[k])
			g.close();
	return cTF
        
def WriteXL(X, L, TF, TG, Name, cutoff=0.1):
	f = open('./PseudoBulk/'+Name+'_CellType.txt')
	ct = f.readlines();f.close()
	ct = [ct[c].strip('\n') for c in range(len(ct))]
	TRS = [np.loadtxt('./CSI/'+ct[c]+'_TRS.txt') for c in range(len(ct))]
    
	os.makedirs('./Results/'+Name, exist_ok=True)
	g = open('./Results/'+Name+'/X.txt','w')
	g.write('TF\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for j in range(X.shape[1]):
			g.write('\t'+str(round(X[i][j],8)))
		g.write('\n')
	g.close()
	cTFs = TFPairs(X, TF, Name)

	g = open('./Results/'+Name+'/A.txt','w')
	g.write('CellType\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for c in range(len(L)):
		L[c] = np.diag(L[c])
		g.write(ct[c])
		for j in range(X.shape[1]):
			g.write('\t'+str(round(L[c][j],8)))
		g.write('\n')
	g.close()

	os.makedirs('./Results/'+Name+'/Annotation/', exist_ok=True)
	for c in range(len(ct)):
		f = open('./Networks/'+ct[c]+'/'+ct[c]+'_network.txt')
		net = f.readlines();f.close();del net[0]
		for i in range(len(net)):
			net[i] = net[i].strip('\n').split('\t')
			net[i][4] = net[i][4].split(';')
		TGX = np.dot((X.dot(np.diag(L[c]))).T,TRS[c])
		TGX = (TGX-np.min(TGX))/(np.max(TGX)-np.min(TGX))
		for j in range(X.shape[1]):
			if L[c][j] >= cutoff:
				W = TGX[j].copy();me = np.mean(W);sd = np.std(W)
				cTG  = []
				for k in range(len(TG)):
					x = (W[k]-me)/sd
					p = 1-stats.norm.cdf(x)
					if p <= 0.05 or W[k]>=0.95:
						cTG.append(TG[k])
				subnet = []
				for k in range(len(net)):
					if net[k][0] in cTFs[j] and net[k][1] in cTG:
						RE = []
						for l in range(len(net[k][4])):
							if net[k][4][l] != "":
								RE.append(net[k][4][l])
						if len(RE) > 0:
							RE = ';'.join(RE)
							subnet.append(net[k][0]+'\t'+net[k][1]+'\t'+net[k][2]+'\t'+RE+'\n')
				if len(subnet) > 0:
					g = open('./Results/'+Name+'/Annotation/'+ct[c]+'_M'+str(j+1)+"_subnetwork.txt",'w')
					for k in range(len(subnet)):
						g.write(subnet[k])	
					g.close()
	for i in range(X.shape[1]):
		D1 = []
		for j in range(X.shape[0]):
			if X[j,i] > 0.05:
				D1.append(X[j,i])
		me = np.mean(D1);sd = np.std(D1)
		for j in range(X.shape[0]):
			X[j,i] = (X[j,i]-me)/sd
	g = open('./Results/'+Name+'/XN.txt','w')
	g.write('TF\t'+'\t'.join(["M"+str(c+1) for c in range(X.shape[1])])+'\n')
	for i in range(len(TF)):
		g.write(TF[i])
		for j in range(X.shape[1]):
			g.write('\t'+str(X[i][j]))
		g.write('\n')
	g.close()

def fit_line_through_points(points):
    x_coords = np.array([point[0] for point in points])
    y_coords = np.array([point[1] for point in points])
    A = np.vstack([x_coords, np.ones(len(x_coords))]).T
    m, b = np.linalg.lstsq(A, y_coords, rcond=None)[0]
    return m

def SelectK(C,MinK,MaxK):
	if MinK != MaxK:
		print("Selecting best number of cRegulons...")
		FinalLoss = [];FinalK = []
		for kk in range(MinK,MaxK+1):
			K = kk
			Err, X, L = RunModuleK(C,K)
			FinalLoss.append(Err);FinalK.append(K)
		slope = []
		for i in range(MaxK-MinK):
			if i < MaxK-MinK-1:
				points = [[FinalK[i],FinalLoss[i]],[FinalK[i+1],FinalLoss[i+1]],[FinalK[i+2],FinalLoss[i+2]]]
				slope.append(fit_line_through_points(points))
			else:
				points = [[FinalK[i],FinalLoss[i]],[FinalK[i+1],FinalLoss[i+1]]]
				slope.append(fit_line_through_points(points))
		diff = []
		for i in range(MaxK-MinK):
			if i==0:
				diff.append(np.abs(slope[i]))
			else:
				diff.append((slope[i]-slope[i-1])*(float(slope[i]<0)*2-1)/100)
		KF = 0
		for i in range(MaxK-MinK-1):
			if diff[i+1]<1:
				KF = FinalK[i];break;
		print("The best number of cRegulon is:", KF)
		return KF
	else:
		KF = MinK
		return MinK

def mode_model(args):
	CSIF, TFF, TGF = CSI(args.name)
	AllTFAS = ReadTFAS(args.name,TFF)
	CSISF = CSIS(CSIF,AllTFAS,args.name)
	if args.module_number != -1:
		lossF, XF, LF = RunModuleK(CSISF,args.module_number)
		WriteXL(XF, LF, TFF, TGF, args.name)
	elif args.module_max != -1 and args.module_min != -1:
		SK = SelectK(CSISF,args.module_min,args.module_max)
		lossF, XF, LF = RunModuleK(CSISF,SK)
		WriteXL(XF, LF, TFF, TGF, args.name)
	else:
		print("Please provide the TF module number or its range!")

def Annot(name,path_to_rna,K):
    f = open(path_to_rna)
    E = f.readlines();f.close()
    cell = E[0].strip("\n").split("\t");del cell[0];del E[0]
    Gene = []
    for i in range(len(E)):
        E[i] = E[i].strip('\n').split('\t')
        Gene.append(E[i][0]);del E[i][0]
    E = np.array(E).astype('float')

    f = open("./Data/HumanAtlasX.txt")
    C0 = f.readlines();f.close()
    MName = C0[0].strip('\n').split('\t');del MName[0];del C0[0]
    TF = []
    for i in range(len(C0)):
        C0[i] = C0[i].split('\t')
        TF.append(C0[i][0]);del C0[i][0]
    C0 = np.array(C0).astype('float')

    TFI = [[],[],[]]
    for i in range(len(TF)):
        TF[i] = TF[i].strip('\n')
        if TF[i] in Gene:
            TFI[0].append(TF[i])
            TFI[1].append(i)
            TFI[2].append(Gene.index(TF[i]))

    C0 = C0[TFI[1],]

    PS = E.sum(1)/np.sum(E)*1000000
    TG = [];TGI = []
    for i in range(len(Gene)):
        if PS[i]>=1:
            TGI.append(i);TG.append(Gene[i])
    E0 = E[TGI,];E1 = E[TFI[2],]

    E0 = np.log(1 + E0);E0 = E0 / np.max(E0)
    model = NMF(n_components=K, init='random', random_state=6, solver='cd', max_iter=100)
    W0 = model.fit_transform(E0)
    H0 = model.components_
    R0 = np.corrcoef(E1,E0)[0:E1.shape[0],E1.shape[0]:];R0[np.isnan(R0)] = 0.0;R0[R0<0] = 0
    X0 = C0[:,np.random.choice(range(C0.shape[1]),K)]
    A0 = np.ones((X0.shape[1],C0.shape[1]))
    A0 = A0/C0.shape[1]

    mu = pow(LA.norm(C0 - np.dot(X0, A0), ord='fro'), 2) / pow(LA.norm(E0 - R0.T.dot(X0).dot(H0), ord='fro'), 2)
    
    err = 1000;eps = 0.001;
    def Loss(x,r,h,a,e):
        t1 = 1/2 * pow(LA.norm(e - r.T.dot(x).dot(h),ord='fro'), 2)
        t2 = mu/2 * pow(LA.norm(C0 - x.dot(a),ord='fro'), 2)
        return t1 - t2
    def NormX(xx):
        xx = xx/np.sqrt(np.sum(xx**2,0).reshape(1,xx.shape[1]))
        xx = xx/(np.sum(xx,1).reshape(xx.shape[0],1)+eps)
        return xx
    def NormA(aa):
        aa = aa/(np.sum(aa,1).reshape(aa.shape[0],1)+eps)
        return aa


    X = X0.copy();R = R0.copy();H = H0.copy();A = A0.copy()
    loss1 = Loss(X,R,H,A,E0)
    epoch = 0
    while err > 1e-4:
        XNext = X * (R.dot(E0).dot(H.T)+mu*C0.dot(A.T)) / (eps + R.dot(R.T).dot(X).dot(H).dot(H.T)+mu*X.dot(A).dot(A.T))
        XNext = NormX(XNext)
        RNext = R * (X.dot(H).dot(E0.T)) / (eps + X.dot(H).dot(H.T).dot(X.T).dot(R))
        HNext = H * (X.T.dot(R).dot(E0)) / (eps + X.T.dot(R).dot(R.T).dot(X).dot(H))
        ANext = A * (X.T.dot(C0)) / (eps + X.T.dot(X).dot(A))
        ANext = NormA(ANext)

        M = np.zeros((K, K))
        for z in range(K):
            M[z, z] = LA.norm(HNext[z, :])
        HNext = np.dot(LA.inv(M), HNext)

        loss2 = Loss(XNext,RNext,HNext,ANext,E0)
        err = np.abs(loss2-loss1)
        loss1 = loss2
        X = XNext.copy();R = RNext.copy();H = HNext.copy();A = ANext.copy()
    
    os.makedirs('./Results/'+name, exist_ok=True)
    np.savetxt('./Results/'+name+'/X.txt',X,delimiter='\t')
    np.savetxt('./Results/'+name+'/R.txt',R,delimiter='\t')
    np.savetxt('./Results/'+name+'/H.txt',H,delimiter='\t')
    np.savetxt('./Results/'+name+'/A.txt',A,delimiter='\t')

def mode_annot(args):
    Annot(args.name, args.path_rna, args.module_number)

def main():
	parser = argparse.ArgumentParser(description="cRegulon software with four modes: prep, grn, model, and annot")
	subparsers = parser.add_subparsers(dest="mode", help="Select the mode")

	# Preprocessing mode
	parser_prep = subparsers.add_parser('prep', help="Preprocessing mode")
	parser_prep.add_argument('--name','-n',type=str, default = "Run",required=False,help="Task name")
	parser_prep.add_argument('--rna','-r',type=str,default = "",required=True,help='Path to the scRNA-seq data file')
	parser_prep.add_argument('--rna_meta','-rm',type=str, default = "",required=True,help='Path to the scRNA-seq meta file')
	parser_prep.add_argument('--atac','-a',type=str,default = "",required=True,help='Path to the scATAC-seq data file')
	parser_prep.add_argument('--atac_meta','-am',type=str, default = "",required=True,help='Path to the scATAC-seq meta file')
	parser_prep.add_argument('--species','-g',type=str, default = "",required=True,help='The analyzing species')
	parser_prep.set_defaults(func=mode_prep)

	# Mode grn
	parser_grn = subparsers.add_parser('grn', help="GRN mode")
	parser_grn.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
	parser_grn.add_argument('--celltype','-ct',type=str,default = "",required=True,help='Cell type name')
	parser_grn.add_argument('--genome','-g',type=str, default = "",required=True,help='Genome build')
	parser_grn.add_argument('--cores','-p',type=int,default = 4,required=True,help='Path to the scATAC-seq data file')
	parser_grn.set_defaults(func=mode_grn)

	# Mode model
	parser_model = subparsers.add_parser('model', help="Model mode")
	parser_model.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
	parser_model.add_argument('--module_number','-mn',type=int,default = -1,required=False,help='The number of TF modules')
	parser_model.add_argument('--module_max','-mmax',type=int,default = -1,required=False,help='The maximum number of TF modules')
	parser_model.add_argument('--module_min','-mmin',type=int,default = -1,required=False,help='The minimum number of TF modules')
	parser_model.set_defaults(func=mode_model)

	# Mode annot
	parser_annot = subparsers.add_parser('annot', help="Annotation mode")
	parser_annot.add_argument('--name','-n',type=str, default = "Run",required=True,help="Task name")
	parser_annot.add_argument('--path_rna','-rna',type=str,default = '',required=True,help='The number of TF modules')
	parser_annot.add_argument('--module_number','-mn',type=int,default = 25,required=True,help='The maximum number of TF modules')
	parser_annot.set_defaults(func=mode_annot)

	args = parser.parse_args()
	if args.mode:
		args.func(args)
	else:
		parser.print_help()

if __name__ == "__main__":
	main()
