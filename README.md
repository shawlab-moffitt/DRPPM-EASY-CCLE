# DRPPM-EASY-CCLE

This is an extention of the [DRPPM Expression Analysis ShinY (EASY) App](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY) which integrates data from the [Cancer Cell Line Encyclopedia (CCLE)](https://sites.broadinstitute.org/ccle/) for sample selection perform differential gene expression and gene set enrichment analyses. This R Shiny app is very similar in features to the main app, expect for the addition of the first tab which allows for CCLE sample selection. Based on this selection, the expression and meta data will be subset and imported in the back end to the app for further analysis and visualization. 

# Installation

* Download ZIP file from https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE
* Unzip and load into directory as a project in R Studio
* Open the ‘App.R’ script and write in user input files and options as directed at the top of the script
  * ‘App.R’ script begins with example files loaded in from the ExampleData folder
* Press ‘Run App’ button in R Studio to run in application or browser window and enjoy!
  * The app script will install any missing packages that the user may not have locally

# Requirements

* `R` - https://cran.r-project.org/src/base/R-4/
* `R Studio` - https://www.rstudio.com/products/rstudio/download/

# R Dependencies

|  |  |  |  |  |
| --- | --- | --- | --- | --- |
| shiny_1.6.0 | shinythemes_1.2.0 | shinyjqui_0.4.0 | shinycssloaders_1.0.0 | tools_4.1.0 |
| dplyr_1.0.7 | tidyr_1.1.3 | readr_2.0.1 | tibble_3.1.3 | DT_0.18 |
| ggplot2_3.3.5 | plotly_4.9.4.1 | enrichplot_1.12.2 | pheatmap_1.0.12 | ggrepel_0.9.1 |
| enrichR_3.0 | limma_3.48.3 | clusterProfiler_4.0.5 | limma_3.48.3 | GSVA_1.40.1 |
| BiocManager_1.30.16 | reshape2_1.4.4 | ggpubr_0.4.0 |  |  |


# Required Files

## Required By User

* **GMT file or Gene Set Data (optional):**
  * If the user chooses to user their own gene set file it must be formatted correctly.
    * If using a .gmt file you can find example formatting by the Broad Institute as seen [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).
    * The other option, usefull if generating your own gene sets, is making a two column tab delimited file with the first column being the gene set name repeating for every gene symbol that would be placed in the second column. Examples of this format can be seen in some of the gene set files we have provided [here](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/tree/main/GeneSets).
  * If the user chooses to use their own gene set file, it is recommended that they use the Getting Started Script for the main EASY app, [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/GettingStartedScripts/GeneSetRDataListGen.R), to generate an R data list which is needed to perform ssGSEA analysis.
  * To simplify this optional input there is a tab within the app's GSEA section for the user to upload their own gene set file instead of hard coding it in.
 * **CCLE Data:**
  * This can be found in the [CCLE_data folder](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/tree/main/CCLE_data)
  * The CCLE data was obtained from the [DepMap Portal](https://depmap.org/portal/download/)
    * We selected the expression data, Mutation data, and Sample Info
    * This data was parse to subset only similar sample name and the mutation data was processed for each cell line annotating if they have the mutation and if so, is it damaging or not.
  * This data is used for the sample selection and allows users to select the type or lineage of cancer and a certain meta variable to group the samples by.

## Required and Provided

* **MSigDB Files:** 
  * These gene set files were gathered from the [Molecular Signatures Database (MSigDB)](http://www.gsea-msigdb.org/gsea/msigdb/index.jsp) as separate collections and processed through R to generate a master gene set file to use for GSEA and ssGSEA analysis.
  * These files begin with "msigdb_" and can be found [here](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/tree/main/GeneSets).
    * Please note that some gene sets are available for *Homo sapiens* and *Mus musculus* which are designated by HS or MM respectively.
* **Tab 2 Gene Set files:**
  * The tab 2 gene set initially written to show gene sets from the Cell Marker Database but can be adjusted by the user
  * We also provide LINCS L1000 gene sets derived from small molecule perturbations following drug treatment.
  * If the user choses to adjust this gene set they must also ensure there is an RData list file provided for it as well.
    * This file can be generated with one of the [getting started scripts](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/GettingStartedScripts/GeneSetRDataListGen.R) that is described on the main EASY app's page.

### Generating a Gene Set RData List for ssGSEA

If you choose to use your own Gene Sets either in .gmt or tab delimited format as described above, in order to perform ssGSEA analysis the Gene Set must be converted into an RData list object when loaded into the app. This list can take several minutes to generate depending on the size of the gene set file, so there is a separate script from the main EASY app to perform this with [GeneSetRDataListGen.R](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY/blob/main/GettingStartedScripts/GeneSetRDataListGen.R). The only user input is the Gene Set file path and name, whether or not it has a header, and the desired outfile path and name. Once they are input, the code can be run as a whole and it will produce and save an RData list which can be input to the R Shiny app. These RData lists have already been generated for provided Gene Sets and can be found in the [Gene Sets folder](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/tree/main/GeneSets) of this repsitory.

# Prepping the R Shiny App

With the main purpose of this app to pre-load the CCLE data for comparison, there does not need to be any user adjustment to the scrip unless you would like to perform GSEA with a gene set file of your choice. Even so, there is an option to upload your own gene set file when running the app so you do not need to hard code it in.

# App Features

With this being an extention of the DRPPM-EASY app, many of the features are similar other than the CCLE sample selection which we will feature below. For more description on the features please see the main [GitHub Page](https://github.com/shawlab-moffitt/DRPPM-EASY-ExprAnalysisShinY).

## CCLE Sample Selection

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/blob/main/Example_App_Pictures/EASE_CCLE_Samples.png?raw=true)

1. Users may select samples based on cancer lineage or disease type.
   * Based on this selection different drop down boxes will appear
2. Once a lineage or disease is selected, users must select a codition to group the samples by.
   * If the user chooses to select a gene of interest to view the variantion type a drop down box of the genes available will be displayed where the user may select one.
3. Once the choices are made the user must update the data selection by pressing the button shown.
4. When the data is updated the meta data and expresison matirx that will be used will be show below.

### CCLE Sample Selection and Ouput example

![alt text](https://github.com/shawlab-moffitt/DRPPM-EASY-CCLE/blob/main/Example_App_Pictures/CCLE_example_output.png?raw=true)

A. In this use case of the App we have selected all sub-liniages of the Ovary which are catagorized based on their TP53 mutation status.

B. These samples were used in single sample GSEA (ssGSEA) with a gene set defining DNA damage response by Amundson *et al.*

C. Next we selected non-small cell lung cancer samples to catagorize based on their KRAS mutation status.

D. We then performed GSEA with the MSigDB Hallmark gene set for unfolded protein response on the subset samples.

E. Additionally, we could visualize the ssGSEA score related to GOBP Protein LLocalization to Cytoplasmic Stress Granule formation.


# Future Enhancments

* Allows for download of the selcted CCLE data
* Adjust sample names. 
  * They currently represent the DepMap ID but we would like to add the cell line name for easier identification
* Possibly gether more sample information to allow for a more extensive codition selection criteria
* Incorporate other data set from large cohort studies
* Add further survival and splice junction analysis



# Quesions and Comments

Please email Alyssa Obermayer at alyssa.obermayer@moffitt.org if you have any further comments or questions in regards to the R Shiny Application.
