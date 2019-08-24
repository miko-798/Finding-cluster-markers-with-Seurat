# Finding-cluster-markers-with-Seurat/Correcting-batch-effects

Lab: Dr. Jill Mesirov, UCSD School of Medicine
Oct. 2018 - June 2019

Hello! 

This README explains the contents in this github repo, to help you better understand how the slides and scripts are organized and assist you to reproduce my results in case you’re interested.

### To get an overview

Please click on Final presentation, if you’d like to see the most recent slides I made that summarize the results - TIL markers found by combining Zheng and Guo datasets (from Dr. Zemin Zhang's lab, cited below), starting from preprocessing the raw counts.

### Folder - Code

This folder contains all the scripts I wrote to do the analysis, including Seurat batch correction, Scanorama batch correction, as well as other helper scripts and functions. 

In Seurat subfolder, please click on “most_recent_code” if you’d like to reproduce the markers found starting from preprocessing the raw counts of Zheng and Guo data. “compare seurat and conos markers” contains the script to find overlapping markers found by both methods, and a .pkl file recording these overlaps.

In Scanorama subfolder, it contains experiments to do batch correction using the method Scanorama. The notebook has more detailed descriptions - adding gaussian noise to one/two batch, adding batch effects to one/two batch, etc.

### Folder - Markers

In this folder, you’ll find the names of the cluster markers saved in .txt format. Please view these markers along with the slides that will tell you which cluster corresponds to which cell type/where the cluster is located on the UMAP.

### Folder - Data

These are direct download from the GEO website of the two papers:

Zheng et al
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638

Guo et al
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99254

---

#### References:

1.	Zheng, C. et al. Landscape of infiltrating T cells in liver cancer revealed by single-cell sequencing. Cell 169, 1342–1356.e16 (2017). DOI: 10.1016/j.cell.2017.05.035
2.	Xinyi Guo et al. Global characterization of T cells in non-small-cell lung cancer by single-cell sequencing, Nature Medicine (2018). DOI: 10.1038/s41591-018-0045-3
3.	Andrew Butler et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species, Nature Biotechnology, Volume36, p411–420 (2018)
4.	Brian Hie et al. Panoramic stitching of heterogeneous single-cell transcriptomic data. bioRxiv. DOI: 10.1101/371179

