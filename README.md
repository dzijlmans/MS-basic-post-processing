So you analyzed your mass spectrometry raw files and want to some post-processing steps. Look no further. This script will generate QC plots and volcano plots highlighting your proteins of interest.

## Getting started

Getting started is easy. Download the R script and analysis parameters file to your working directory. To do this:
- Click on the big green button in the top right that says ```<> Code``` and then click on ```Download ZIP```.
- Place the files in your working directory
- Fill in the analysis parameter file accordingly.
- Open the R script with RStudio, select all and [run it](https://www.youtube.com/watch?v=w6QGe-pXgdI).

**Note:** - Installing the DEP package fr the first time may take a while. Should be smooth sailing after that

## Parameters explained

Explanations are also present in the Excel file, when hovering over the parameter boxes.

```analysis_method``` - One of "DiaNN", "MaxQuant" or "Proteome discoverer".

```conditions``` - The amount of experimental variables. This can be either 2 (Condition A, Replicates) or 3 (Condition B, Condition A, Replicates). These variables should be present in the .raw file names.

```filtering_type``` - Whether a protein should have a value in all samples ("complete") or all samples in at least one condition/group ("condition").

```mass_spec``` - This should be present in the names of the .raw files. Required when selecting DiaNN as analysis method.

```comparison``` - The type of comparison that will be tested. This can be all possible pairwise comparisons ("all") or limited to the comparisons versus the control ("control").

```control``` - In case you choose "control" as your comparison, add the name of the control sample, e.g. DMSO, IgG, untreated, KO.

```unique_peptides``` - The minimum amount of unique peptides a protein should be quantified with.

```volcano``` - Specifiy which proteins to highlight in the volcano plot. Can be either a list of supplied protein names ("protein list"), all proteins above a supplied p-value and fold change cutoff ("specify significance") or the top N most differential proteins below a set p-value ("TopN").

```proteins_to_highlight``` - list of protein names to be highlighted, separated by a semicolon (;).

```p_value``` - p-value cutoff for significance. **REQUIRED**

```log2_FC``` - log2 FC value cutoff for significance. **REQUIRED**

```TopN``` - TopN number of proteins to highlight.



### Credit

This script uses a lot of the DEP package, so credit to them for developing the package [(Zhang X, Smits A, van Tilburg G, Ovaa H, Huber W, Vermeulen M (2018). “Proteome-wide identification of ubiquitin interactions using UbIA-MS.” Nature Protocols, 13, 530–550)](https://www.nature.com/articles/nprot.2017.147).
