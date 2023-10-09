# README for Atlantic Salmon vs Coho Salmon snRNAseq

**By Sarah Salisbury**

&#x1F4D8; This repository contains the scripts used to analyze snRNAseq libraries generated from skin and fin samples from Atlantic Salmon (_Salmo salar_) and Coho Salmon (_Oncorhynchus kisutch_) for the manuscript titled "Keratinocytes Drive the Epithelial Hyperplasia key to Sea Lice Resistance in Coho Salmon".

## Repository Architecture and Scripts


&#128194; - **directory**

&#128195; - **script/code**

&#128203; - **data table**

<pre>
├──&#128194;<b>STEP_1_STAR</b> <- contains scripts for analysis using STAR
│   ├──&#128195;<b>1_STAR_AS.md</b> <- script for indexing Atlantic Salmon genome and mapping Atlantic Salmon snRNAseq libraries with STAR
│   ├──&#128195;<b>2_STAR_CO.md</b> <- script for indexing Coho Salmon genome and mapping Coho Salmon snRNAseq libraries with STAR
├──&#128194;<b>STEP_2_SEURAT</b> <- contains scripts for analysis using Seurat
│   ├──&#128194;<b>ATLANTIC_SALMON</b> <- contains scripts used to analyze Atlantic Salmon samples and generate cell clustering for this species
│   │   ├──&#128195;<b>1_Initial_QC_and_Clustering.Rmd</b> <- R Markdown script to do initial QC and cell clustering for Atlantic Salmon samples
│   │   ├──&#128195;<b>2_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type from initial clustering
│   │   ├──&#128194;<b>3_Filter1</b> <- contains scripts used to remove dubious clusters after initial clustering and then recluster remaining cells
│   │   │   ├──&#128195;<b>1_Filter1_Recluster.Rmd</b> <- R Markdown script to remove dubious clusters after initial clustering and then recluster remaining cells
│   │   │   ├──&#128195;<b>2_Filter1_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after filtering and reclustering cells
│   │   │   ├──&#128194;<b>3_Subcluster_Immune_Cells_Only</b> <- contains scripts to subcluster those clusters identified as immune cells
│   │   │   │   ├──&#128195;<b>1_Filter1_Immune_Cells_Only_Recluster.Rmd</b> <- R Markdown script to take only immune cells and recluster them
│   │   │   │   ├──&#128195;<b>2_Filter1_Immune_Cells_Only_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after reclustering immune cells
│   │   │   │   ├──&#128195;<b>3_Filter1_Immune_Cells_Only_ID_Immune_Cell_Types.Rmd</b> <- R Markdown script to generate list of cell barcodes for each cell type after reclustering immune cells
│   │   │   ├──&#128194;<b>4_Filter1_Incorporating_Immune_Cell_Subclustering</b> <- contains scripts to incorporate immune cell subclusters into larger Seurat object with all other cell types
│   │   │   │   ├──&#128195;<b>1_Filter1_Incorporate_Immune_Substructuring_Cell_Types.Rmd</b> <- R Markdown script to assign immune cells to appropriate subcluster ID within the larger Seurat object containing all cell types
│   │   │   │   ├──&#128195;<b>2_Filter1_Post_Immune_Cell_Incorporation_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after incorporating immune cell subclusters
│   │   │   │   ├──&#128195;<b>3_Filter1_Post_Immune_Cell_Incorporation_Differential_Expression.Rmd</b> <- R Markdown script to generate lists of differentially expressed genes between the control and each infection time point for each cell type
│   │   │   │   ├──&#128195;<b>4_Filter1_Filter_Differential_Expression_Results.Rmd</b> <- R Markdown script to filter genes identified as differentially expressed
│   │   ├──&#128203;<b>FEATUREEnsembltosymbolAS.csv</b> Table listing ENSEMBL ID (V1) and Seurat-assigned gene symbol (V2)
│   ├──&#128194;<b>CELL_CYCLE_GENES</b> <- contains tables listing genes used for Cell Cycle Scoring analyses
│   │   ├──&#128203;<b>CellCycleGenesCOVRAI.csv</b> <- table listing cell cycle genes used for Seurat analyses of Coho Salmon samples only
│   │   ├──&#128203;<b>CellCycleGenesCOVRAI_ASCOconversionforCO_proof_final.csv</b> <- table listing cell cycle genes used for Seurat analyses of Coho Salmon samples in preparation for integrating Atlantic Salmon and Coho Salmon samples together
│   │   ├──&#128203;<b>CellCycleGenesV3ENVRAI.csv</b> <- table listing cell cycle genes used for Seurat analyses of Atlantic Salmon samples only
│   ├──&#128194;<b>COHO_SALMON</b> <- contains scripts used to analyze Coho Salmon samples and generate cell clustering for this species
│   │   ├──&#128195;<b>1_Initial_QC_and_Clustering.Rmd</b> <- R Markdown script to do initial QC and cell clustering for Atlantic Salmon samples
│   │   ├──&#128195;<b>2_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type from initial clustering
│   │   ├──&#128194;<b>3_Filter1</b> <- contains scripts used to remove dubious clusters after initial clustering and then recluster remaining cells
│   │   │   ├──&#128195;<b>1_Filter1_Recluster.Rmd</b> <- R Markdown script to remove dubious clusters after initial clustering and then recluster remaining cells
│   │   │   ├──&#128195;<b>2_Filter1_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after filtering and reclustering cells
│   │   ├──&#128194;<b>4_Filter2</b> <- contains scripts used to remove dubious clusters after Filter1 and then recluster remaining cells
│   │   │   ├──&#128195;<b>1_Filter2_Recluster.Rmd</b> <- R Markdown script to remove dubious clusters after first filtration and then recluster remaining cells
│   │   │   ├──&#128195;<b>2_Filter2_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after filtering and reclustering cells
│   │   │   ├──&#128194;<b>3_Subcluster_Immune_Cells_Only</b> <- contains scripts to subcluster those clusters identified as immune cells
│   │   │   │   ├──&#128195;<b>1_Filter2_Immune_Cells_Only_Recluster.Rmd</b> <- R Markdown script to take only immune cells and recluster them
│   │   │   │   ├──&#128195;<b>2_Filter2_Immune_Cells_Only_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after reclustering immune cells
│   │   │   │   ├──&#128195;<b>3_Filter1_Immune_Cells_Only_ID_Immune_Cell_Types.Rmd</b> <- R Markdown script to generate list of cell barcodes for each cell type after reclustering immune cells
│   │   │   ├──&#128194;<b>4_Cluster_12_Only</b> <- contains scripts to incorporate subcluster cluster 12
│   │   │   │   ├──&#128195;<b>1_Filter2_Cluster12_Only_Recluster.Rmd</b> <- R Markdown script to take only cluster 12 cells and recluster them
│   │   │   │   ├──&#128195;<b>2_Filter2_Cluster12_Only_ID_Cell_Types.Rmd</b> <- R Markdown script to identify marker genes for each cell type after reclustering cluster 12 cells
│   │   │   ├──&#128194;<b>5_Filter2_Incorporating_Immune_Cell_Subclustering</b> <- contains scripts to incorporate immune cell subclusters into larger Seurat object with all other cell types
│   │   │   │   ├──&#128195;<b>1_Filter2_Incorporate_Immune_and_Cluster12_Substructuring_Cell_Types.Rmd</b> <- R Markdown script to assign immune cells and cluster 12 cells to subcluster ID within the larger Seurat object containing all cell types
│   │   │   │   ├──&#128195;<b>2_Filter2_Post_Immune_Cell_and_Cluster12_Incorporation_Marker_Gene_Detection.Rmd</b> <- R Markdown script to identify marker genes for each cell type after incorporating immune cell and cluster 12 cell subclusters
│   │   │   │   ├──&#128195;<b>3_Filter2_Post_Immune_Cell_and_Cluster12_Incorporation_Differential_Expression.Rmd</b> <- R Markdown script to generate lists of differentially expressed genes between the control and each infection time point for each cell type
│   │   │   │   ├──&#128195;<b>4_Filter1_Filter_Differential_Expression_Results.Rmd</b> <- R Markdown script to filter genes identified as differentially expressed
│   │   ├──&#128203;<b>FEATUREEnsembltosymbolCO.csv</b> Table listing ENSEMBL ID (V1) and Seurat-assigned gene symbol (V2)
│   ├──&#128194;<b>INTEGRATE_ATLANTIC_AND_COHO_SALMON</b> <- contains scripts used to analyze Atlantic Salmon and Coho Salmon samples together and generate cell clustering combining cells from both species
│   │   ├──&#128195;<b>1_Atlantic_Salmon_Preparation.Rmd</b> <- R Markdown script to do initial QC and cell clustering for Atlantic Salmon samples
│   │   ├──&#128195;<b>2_Coho_Salmon_Preparation.Rmd</b> <- R Markdown script to identify marker genes for each cell type from initial clustering
│   │   ├──&#128195;<b>3_Integrating_Atlantic_and_Coho_Salmon_Data.Rmd</b> <- contains scripts used to remove dubious clusters after initial clustering and then recluster remaining cells
│   │   ├──&#128195;<b>4_Marker_Gene_Detection_for_Integrated_Data.Rmd</b> <- contains scripts used to remove dubious clusters after initial clustering and then recluster remaining cells
│   ├──&#128194;<b>ORTHOGROUPS</b> <- contains table listing 1:1 orthologs between Atlantic Salmon and Coho Salmon
│   │   ├──&#128203;<b>Orthogroups_1geneperspecies.tsv</b> <- table listing 1:1 orthologs between Atlantic Salmon and Coho Salmon
│   ├──&#128194;<b>PLOTS</b> <- contains examples of scripts used to generate various types of plots to illustrate data
│   │   ├──&#128194;<b>Dotplots</b> <- contains scripts used generate dotplot figures
│   │   │   ├──&#128195;<b>Genes_vs_cell_types_dotplot_Atlantic_Salmon.Rmd</b> <- R Markdown script to generate dotplot of gene expression of a list of genes in all cell types in Atlantic Salmon
│   │   │   ├──&#128195;<b>Genes_vs_cell_types_dotplot_Coho_Salmon.Rmd</b> <- R Markdown script to generate dotplot of gene expression of a list of genes in all cell types in Coho Salmon
│   │   │   ├──&#128195;<b>Timepoints_vs_cell_types_dotplot_multiple_cell_types.Rmd</b> <- R Markdown script to generate dotplot comparing gene expression of a given gene at each timepoint and in multiple cell types
│   │   │   ├──&#128195;<b>Timepoints_vs_cell_types_dotplot_single_cell_type.Rmd</b> <- R Markdown script to generate dotplot comparing gene expression of a given gene at each timepoint and in a single cell type
│   │   ├──&#128194;<b>Stacked_bar_plots_cell_counts</b> <- contains scripts used generate bar plots illustrating the number of cells from each sample in a particular cell type
│   │   │   ├──&#128195;<b>Atlantic_Salmon_stacked_bar_plot_cell_count.Rmd</b> <- R Markdown script to generate bar plot illustrating the number of cells from each Atlantic Salmon sample in a particular cell type
│   │   │   ├──&#128195;<b>Coho_Salmon_stacked_bar_plot_cell_count.Rmd</b> <- R Markdown script to generate bar plot illustrating the number of cells from each Coho Salmon sample in a particular cell type
│   │   ├──&#128194;<b>UMAPS</b> <- contains scripts used illustrate UMAPs of cell types
│   │   │   ├──&#128195;<b>Atlantic_Salmon_UMAP.Rmd</b> <- R Markdown script to generate UMAP plot for Atlantic Salmon samples
│   │   │   ├──&#128195;<b>Coho_Salmon_UMAP.Rmd</b> <- R Markdown script to generate UMAP plot for Coho Salmon samples
│   │   ├──&#128194;<b>Violin_Plots</b> <- contains scripts used generate various violin plots (of gene expression, umi/feature counts, etc.)
│   │   │   ├──&#128195;<b>Cell_type_and_sample_violin_plots_Atlantic_Salmon.Rmd</b> <- R Markdown script to generate violin plots of umi/feature counts per cell type/sample and the expression per cell type of the top 20 marker genes for each cell type using Atlantic Salmon samples
│   │   │   ├──&#128195;<b>Cell_type_and_sample_violin_plots_Coho_Salmon.Rmd</b> <- R Markdown script to generate violin plots of umi/feature counts per cell type/sample and the expression per cell type of the top 20 marker genes for each cell type using Coho Salmon samples
│   │   │   ├──&#128195;<b>Stacked_Violin_Plot_Atlantic_Salmon.Rmd</b> <- R Markdown script to generate multiple violin plots of the expression of a list of genes in all Atlantic Salmon cell types
│   │   │   ├──&#128195;<b>Stacked_Violin_Plot_Coho_Salmon.Rmd</b> <- R Markdown script to generate multiple violin plots of the expression of a list of genes in all Coho Salmon cell types
</pre>

## How to Navigate this Repository

STAR analyses (located within STEP_1_STAR) must be completed before Seurat analyses (located within STEP_2_SEURAT).

Within each directory, the order in which scripts (or the scripts within a subdirectory) should be completed is indicated by the number at the beginning of the script/subdirectory name.

For example, the directory ```STEP_2_SEURAT/ATLANTIC_SALMON``` has within it the following scripts and subdirectories:

&#128195; 1_Initial_QC_and_Clustering.Rmd

&#128195; 2_Marker_Gene_Detection.Rmd

&#128194;3_Filter1

&#128194;4_Filter1_Incorporate_Immune_Cell_Substructuring

Therefore, for this directory, script ```1_Initial_QC_and_Clustering.Rmd``` should be completed first, followed by the script ```2_Marker_Gene_Detection.Rmd```, followed by _all_ scripts within the subdirectory ```3_Filter1```, followed by _all_ scripts within the subdirectory ```4_Filter1_Incorporate_Immune_Cell_Substructuring```. Please note that within subdirectories the scripts/subsubdirectories will also be numbered in this manner. Therefore, scripts may be conducted in the order they are mentioned within the "Repository Architecture and Scripts" section above.
