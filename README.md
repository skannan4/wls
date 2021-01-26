# Heart fields auto-regulate second heart field cell fate via Wnt secretion

### Introduction
We performed an scRNA-seq study to identify a possible source of Wnt signals that influence second heart field cell fate. To this end, we deleted Wls, a gene required for Wnt secretion, from Mesp1-lineage cells, thereby eliminating Wnt secretion in mesodermal populations. We used this data to identify which mesodermal populations expression Wnt ligands, and determine the downstream effects of Wls knockout. Our results can be found in our preprint (add when ready). Here, we share all of the materials needed to reproduce our analysis.

### Data
All of the relevant files to reproduce the analysis can be downloaded at our [Synapse](https://www.synapse.org/#!Synapse:syn24200678/files/). In particular, the following files may be of relevance:

- `cells_x_genes.mtx`, `cells_x_genes.genes.txt`, and `cells_x_genes.barcodes.txt`: These three files can be found in each of the subfolders starting with "MMiy", each corresponding to a separate 10x experimental run. These are the direct mapped outputs of kallisto|bustools, for those who want to start from the raw mapped counts and pursue their own integration pipelines.

- `combined_data_sparse.mtx`, `combined_data_sparse_genes.txt`, and `combined_data_sparse_cells.txt`: These are files associated with a combined counts table for the three runs (using the kallisto|bustools output). In our pipeline, we input this counts table to Seurat and then integrated using the SCTranform + Integration workflow.

- `combined_pheno.txt` and `seurat_combined_pheno.txt`: Metadata for the cells in the study. The first is a simplified metadata table for all of the cells in the study (including information about run, sample barcode, and experimental condition for each cell). The latter is a more detailed metadata table, and also contains the final Seurat cluster numbers and identified celltype annotations for each cell. Note, however, that not every cell in `combined_pheno.txt` is in `seurat_combined_pheno.txt` - this is not only because poor quality cells were filtered, but also because we filtered RBCs from our analysis. If users are interested in those annotations as well, please contact us directly.

- `wls.RData`: A pre-made workspace that contains a lot of the likely objects of interest for readers. In particular, the final Seurat object can be found as `seurat.integrated.clean`. Additionally, the objects of the trajectory reconstruction in Monocle 2 and subsequent analysis can also be found in this workspace. Note that, for the sake of space, we have removed the Seurat object containing RBCs - however, we are happy to share this with users if desired.

The following relevant files can be found on Github:


- `wls_code.R`: An R file containing the entire analysis workflow for the manuscript. Note that not every file to reproduce this workflow from scratch will be found on Synapse - in particular, the Fastq files associated with the MultiSeq barcodes are not available there. However, we have uploaded that data to [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165300) (GSE165300), so please feel free to grab that data however you choose if you wish to reproduce the entirety of the code. The code has also been largely run and saved in `wls.RData`.


### Dependencies
Most of the libraries used in our codebase can be found from CRAN or Bioconductor. However, we in particular make use of deMULTIplex by Chris McGinnis to work with MultiSeq barcodes. Please see [his github](https://github.com/chris-mcginnis-ucsf/MULTI-seq) for instructions to download. Note that, at one point in the demultiplexing, we needed to use a slightly modified version of deMULTIplex; this code can be found in our [Synapse folder](https://www.synapse.org/#!Synapse:syn24200678/files/). We additionally use the tradeSeq package for trajectory differential gene analysis. Please see [the github](https://github.com/statOmics/tradeSeq) for install instructions. At the time we developed this code, tradeSeq was unable to allow sample name labeling on plots; to this end, I generated my own fork of tradeSeq with a quick modification to allow for plotting. This tweak can be found on [my github](https://github.com/skannan4/tradeSeq); if you wish to fully reproduce the code here, download this version of tradeSeq. However, the new versions of tradeSeq provide more plotting options, so feel free to use their (likely better) version - just keep in mind that some of the plotting code in `wls_code.R` may need to be tweaked appropriately. Lastly, this analysis uses Monocle 2, which technically is "deprecated" but in our hands has been superior to Monocle 3 for certain applications. To get Monocle 2, follow the download instructions [here](http://cole-trapnell-lab.github.io/monocle-release/). You should be able to simultaneously have Monocle 2 and 3 installed.


### How to replicate our workflow
Please note that this workflow is somewhat complicated to reproduce from the very top. This is because the sample demultiplexing steps are somewhat complicated and memory intensive. Likewise, the Seurat integration steps may be somewhat tedious. Thus, we have tried to provide resources to readers so that they can pick up the analysis from a range of different steps. We note those below. To replicate in general, please follow these steps:

1. Download `wls_code.R` and all of the files in the listed Synapse folder.

2. Modify the appropriate line in `wls_code.R`: `setwd("~/Documents/Research/Wls/FinalFiles")` to set the working directory to the same working directory you downloaded the Synapse files into it.

3. If you are interested in replicating the entire workflow, including the demultiplexing steps, download the appropriate Fastq files for the MultiSeq barcode data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165300). Note that for some of the demultiplexing, you will need our custom modification to deMULTIplex, which can be found in the Synpase folder (`classifyCells.R`).

4. If you are interested in replicating the workflow starting from immediately after demultiplexing, simply load in `combined_data.txt` and `combined_pheno.txt` into R and proceed.

Alternatively, note that you can start by loading `wls.RData` into R or an IDE of your choice. This should get you the final Seurat object (`seurat.integrated.clean`) as well as a range of other objects generated from the analysis. This object is compatible with the code in `wls_code.R`, and thus you should be able to immediately re-run any of the downstream analysis or figure plotting steps straight from this workspace.

Please feel free to email or raise an issue if any of the code doesn't work as claimed!
