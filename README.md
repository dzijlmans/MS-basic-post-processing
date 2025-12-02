# MS-basic-post-processing

So you analyzed your mass spectrometry raw files and want to some post-processing steps. Look no further. This script will generate QC plots and volcano plots highlighting your proteins of interest, as well as a browsable Excel file.

## Getting started

Getting started is easy. Download the R script and analysis parameters file to your working directory. To do this:
- Click on the big green button in the top right that says ```<> Code``` and then click on ```Download ZIP```.
- Place the files in your working directory
- Fill in the analysis parameter file accordingly.
- Open the R script with RStudio, select all and [run it](https://www.youtube.com/watch?v=w6QGe-pXgdI).

<ins>**Note:**</ins> - Installing the DEP package for the first time may take a while. Should be smooth sailing after that

## Parameters explained

Explanations are also present in the Excel file, when hovering over the parameter boxes.

```analysis_method``` - One of "DiaNN", "MaxQuant" or "Proteome discoverer".

```conditions``` - The amount of experimental variables. This can be either 2 (Condition A, Replicates) or 3 (Condition B, Condition A, Replicates). These variables should be present in the .raw file names.

```filtering_type``` - Whether a protein should have a value in all samples ("complete") or all samples in at least one condition/group ("condition").

```mass_spec``` - Required when selecting DiaNN as analysis method. This should be present in the names of the .raw files. 

```comparison``` - The type of comparison that will be tested. This can be all possible pairwise comparisons ("all"), a manual selection ("manual") or limited to the comparisons versus the control ("control").

```contrasts``` - In case you choose "manual" as your comparison, add the comparisons you would like to make, separated by a semicolon. E.g. C1_vs_DMSO;DMSO_vs_IgG

```control``` - In case you choose "control" as your comparison, add the name of the control sample, e.g. DMSO, IgG, untreated, KO.

```unique_peptides``` - The minimum amount of unique peptides a protein should be quantified with.

```volcano``` - Specifiy which proteins to highlight in the volcano plot. Can be either a list of supplied protein names ("protein list"), all proteins above a supplied p-value and fold change cutoff ("specify significance") or the top N most differential proteins below a set p-value ("TopN").

```proteins_to_highlight``` - list of protein names to be highlighted, separated by a semicolon (;). E.g. EZH2;MBD3

```p_value``` - p-value cutoff for significance. **REQUIRED**

```log2_FC``` - log2 FC value cutoff for significance. **REQUIRED**

```TopN``` - TopN number of proteins to highlight.

```highlight_imputed``` - If set to TRUE, a secondary volcano plot will be generated, highlighting which proteins were (partially) imputed.

## Script input requirements

```DiaNN``` - Output files should end on ```*pg_matrix.tsv``` and ```*pr_matrix.tsv```.

```MaxQuant``` - Output file should end on ```*proteinGroups.txt```.

```Proteome Discoverer``` - Output file should end on ```*proteins.txt```.

### Abundance/Intensity column requirements

Please ensure that the abundance/intensity column names contain the following information, in order and separated by an underscore:

- <ins>Optional:</ins> Condition B (e.g. Cell line)
- Condition A (e.g. WT and KO or treated and untreated)
- Replicate
             
Example: ```HCT116_untreated_1```.

For DiaNN, this information should also be included in the ```.raw``` files. Example: ```20250205_Astral_AUR_DZ114_Wildtype_untreated_1.raw```

## Per contrast comparisons
From version v0.2.0 onwards, the analysis will run per contrast (= comparison), as opposed to on the dataset as a whole. This was done to avoid unnecessary overfiltering and over- or under-imputing.

If you have multiple contrasts in your experiment, this means less proteins will be filtered out if you filtered for at least one complete condition. Additionally, imputation is more accurate since the imputed values are only based on the conditions involved and missingness patterns reflect only the biology of that comparison.

## Credit

This script uses a lot of the DEP package, so credit to them for developing the package [(Zhang X, Smits A, van Tilburg G, Ovaa H, Huber W, Vermeulen M (2018). “Proteome-wide identification of ubiquitin interactions using UbIA-MS.” Nature Protocols, 13, 530–550)](https://www.nature.com/articles/nprot.2017.147).
