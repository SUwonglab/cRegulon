## An instruction to display cRegulon in our paper
The typical outputs of cRegulon is: <br>
1. X.txt: A TF by module matrix indicating the weights of TFs in each module.
2. A.txt: A cell type by TF module matrix indicating the associations between cell types and cRegulons.
3. Regulon*_TFModule.txt: The TF pairs in TF module of cRegulon*
With these output, we show how to interpret and display the cRegulons.

### Display important TFs in TF modules
We can study the behavior of important TFs in the cRegulons. The important TF can be prior known TF markers or top TFs in TF modules. <br>
Taking cell line experiment in our paper for example, we can run the following sceipt to show important TFs:
```bash
python PlotTFMarker.py
```
### Display important TFs in TF modules
cRegulon will output the TF pairs of the TF modules, such as **CL_cRegulon1_TFModule.txt**. We can display the TF module in the form of network. <br>
Taking cell line experiment in our paper for example, we can use **Cytoscape software (PlotTFModule.cys)** to construct such TF combinatorial network.
### Display cell types and cRegulons associations
Our software will output the associations between cell types and cRegulons. We can display these association by heatmap or the association plot defined in our manuscript. <br>
Taking cell line experiment in our paper for example, we can run the following sceipt to display associations:
```bash
python PlotAssociation.py
```
