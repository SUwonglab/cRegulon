# cRegulon

## Introduction
This is cRegulon software: an optimization model an optimization model for modeling combinatorial regulation from single-cell multi-omics provides units underpinning cell type landscape.
## Requirements:
1. Homer in your environment PATH <br>
2. Python >=3.0 and packages below: <br>
3. numpy-groupies==0.9.15 <br>
4. numpy==1.26.3 <br>
5. scipy==1.12.0 <br>
6. pandas==1.4.3 <br>
7. ismember==1.0.2 <br>
8. mpmath==1.2.1 <br>
9. pybedtools==0.9.0 <br>
10. scikit-learn==1.4.2 <br>

## Installing cRegulon with the following command:
```bash
wget https://github.com/SUwonglab/cRegulon/archive/master.zip
unzip master.zip
cd cRegulon-master
wget -O cRegulonData.tar.gz https://figshare.com/ndownloader/files/47503895
tar -xzvf cRegulonData.tar.gz
```
***We also provide an example for cRegulon analysis, which can be downloaded from the following command:***
```bash
wget -O cRegulonExample.tar.gz https://figshare.com/ndownloader/files/47503904
tar -xzvf cRegulonExample.tar.gz
```
***After downloading, there will be 5 folders that store raw data, intermediate results, and final results of example: <br>
./example_data/ <br>
./CSI/ <br>
./PseudoBulk/<br>
./Networks/<br>
./Results/***

## Tutorial of cRegulon analysis
### Input single cell data of cRegulon
The typic input file of scRNA-seq data is a gene by cell count matrix: <br>
<table>
  <tr>
    <td>scRNA</td>
    <td>RNACellID1</td>
    <td>RNACellID2</td>
    <td>RNACellID3</td>
  </tr>
  <tr>
    <td>Gene1</td>
    <td>5</td>
    <td>0</td>
    <td>3</td>
  </tr>
  <tr>
    <td>Gene2</td>
    <td>0</td>
    <td>2</td>
    <td>0</td>
  </tr>
  <tr>
    <td>Gene3</td>
    <td>1</td>
    <td>0</td>
    <td>0</td>
  </tr>
</table>
The typic input file of scATAC-seq data is a peak by cell count matrix:
<table>
  <tr>
    <td>scATAC</td>
    <td>ATACellID1</td>
    <td>ATACellID2</td>
    <td>ATACellID3</td>
    <td>ATACellID4</td>
  </tr>
  <tr>
    <td>Peak1</td>
    <td>1</td>
    <td>0</td>
    <td>1</td>
    <td>0</td>
  </tr>
  <tr>
    <td>Peak2</td>
    <td>0</td>
    <td>1</td>
    <td>0</td>
    <td>1</td>
  </tr>
  <tr>
    <td>Peak3</td>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>0</td>
  </tr>
</table>
The peaks are in the format of "chr_start_end". <br>
<br>

***In practice, the file of single cell dataset can be the 10x folder format (example dataset: ./example_data/RA/) or a matrix txt format (example dataset: ./example_data/CL/). In this tutorial, we use RA dataset for demonstration.***

### Input cell type meta data
The typic cell type meta file (./example_data/RA/RA_scRNA_Cluster.txt) of scRNA-seq data is as follows: <br>
<table>
  <tr>
    <td>RNACellID1</td>
    <td>RAC1</td>
  </tr>
  <tr>
    <td>RNACellID2</td>
    <td>RAC2</td>
  </tr>
  <tr>
    <td>RNACellID3</td>
    <td>RAC3</td>
  </tr>
</table>
The typic cell type meta file (./example_data/RA/RA_scATAC_Cluster.txt) of scATAC-seq data is as follows: <br>
<table>
  <tr>
    <td>ATACellID1</td>
    <td>RAC2</td>
  </tr>
  <tr>
    <td>ATACellID2</td>
    <td>RAC3</td>
  </tr>
  <tr>
    <td>ATACellID3</td>
    <td>RAC1</td>
  </tr>
  <tr>
    <td>ATACellID4</td>
    <td>RAC2</td>
  </tr>
</table>

### cRegulon has three main modes, which are three steps of cRegulon model
```bash
python3 cRegulon.py {prep,grn,model} ...
```
 prep: Preprocessing mode <br>
 grn: GRN mode <br>
 model: Model mode <br>

### Step 1: preprocessing and pseudo bulk (***prep mode***)
We run the following script to create pseudo bulk RNA-seq and ATAC-seq data for each cell cluster:

```bash
#cRegulon.py prep [-h] [--name NAME] --rna RNA --rna_meta RNA_META --atac ATAC --atac_meta ATAC_META --species SPECIES (human or mouse)
python3 cRegulon.py prep --name RA --rna ./example_data/RA/scRNA/ --rna_meta ./example_data/RA/RA_scRNA_Cluster.txt --atac ./example_data/RA/scATAC/ --atac_meta ./example_data/RA/RA_scATAC_Cluster.txt -g mouse
```
This process will produce pseudo bulk files (*PS_RNA.txt, *PS_ATAC.txt, *CellType.txt) for each cell cluster in the **PseudoBulk** folder.

### Step 2: GRN construction (***grn mode***)
We run the following script to construct regulatory network for each cell cluster (current we support hg38 and mm10):

```bash
#python3 cRegulon.py grn [-h] --name NAME --celltype CELLTYPE --genome GENOME --cores CORES
for c in `cat ./PseudoBulk/RA_CellType.txt`
do
    python3 cRegulon.py grn -n RA -ct ${c} -g mm10 -p 20
done
```
This process will produce GRN files (*network.txt, TFTG_regulationScore.txt, TFName.txt, TGName.txt) for each cell cluster in the **Networks** folder (The GRN construction is independent for each cell cluster, we can do it **parallelly**).

### Step 3: Running cRegulon model (***model mode***)
We run the following script of cRegulon model:
If we already know or have some expection of the cRegulon number, we can provide this number to cRegulon. For example, we have 9 cRegulons for RA, then we run this script:
```bash
#python3 cRegulon.py model [-h] --name NAME --module_number MODULE_NUMBER
python3 cRegulon.py model -n RA -mn 9
```
If we don't know the cRegulon number, we can provide a range of numbers and cRegulon will use elbow rule to select an optimal number. For example, we guess there may be 4-20 cRegulons for RA, then we run this script:
```bash
#python3 cRegulon.py model [-h] --name NAME --module_max MODULE_MAX --module_min MODULE_MIN
python3 cRegulon.py model -n RA -mmin 4 -mmax 20
```
This will output a folder in "Results" with name you specify: "./RA/" <br>
1. TF combinatorial effects in each cRegulon: ./Results/RA/X.txt <br>
2. Association matrix between cell clusters and cRegulons: ./Results/RA/A.txt <br>
3. TF module of each cRegulon: ./Results/RA/*TFModule.txt <br>
4. Annotation of each cell cluster with cRegulons: ./Results/RA/Annotation/*subnetwork.txt

## Annotation mode of cRegulon (***annot mode***)
If you only have scRNA-seq data, We run the following script to annotate cells with our pre-computed cRegulons from atlas-level dataset:
```bash
#python3 cRegulon.py annot [-h] --name NAME --path_rna PATH_RNA --module_number MODULE_NUMBER
python3 cRegulon.py annot --name PBMC --path_rna ./example_data/PBMC/PBMC_scRNA.txt --module_number 12
```
## Citation:
If you use cRegulon software or cRegulon associated concepts, please cite:

Zhanying Feng, et al. Modeling combinatorial regulation from single-cell multi-omics provides units underpinning cell type landscape. 2024.
