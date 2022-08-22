# CTULIP - Canine TUmor CLassIfication Predictor <img src = "images/tulip.svg" alt = "tulip" width = "50" height = "50">

CTULIP (the Canine TUmor cLassIfication Predictor) classifies Canine RNA-Seq samples into tumor types.

## Description

CTULIP is a modified framework of [TULIP](https://github.com/CBIIT/TULIP) which is a 1D convolutional neural network for classifying Canine RNA-Seq data with [14,761 genes](https://github.com/CBIIT/CTULIP/blob/main/gene_lists/caninegenes.txt). CTULIP can classify either list into [17](https://github.com/CBIIT/CTULIP/blob/main/labels/17_tumors.csv) or [18](https://github.com/CBIIT/CTULIP/blob/main/labels/18_tumors.csv) (that includes Glioblastoma multiforme) tumor types. 

The resource transfer team trained and validated the models in CTULIP on over 9,000 TCGA RNA-Seq files from the [Genomic Data Commons](https://portal.gdc.cancer.gov/) (GDC) in February 2022. To use CTULIP, the user must provide a file of RNA-Seq data expressed as FPKM-UQ (fragments per kilobase of transcript per million mapped reads upper quartile) for one or more samples. CTULIP then converts FPKM-UQ values to TPM (transcripts per million), performs log10 normalization, and reformats the data into the correct dimensions before applying the selected model. CTULIP generates two files, as described in the following table:

*The naming convention of the output files reflect the model used and the file contents.*

| Example File for 17 Tumor Types| Example File for 18 Tumor Types| File Contents |
| ------------- | ------------- | ------------- |
| [predictions_17_ctoh_1-1.csv](https://github.com/CBIIT/CTULIP/blob/main/example_results/predictions_17_ctoh_1-1.csv)  | [predictions_18_ctoh_1-1.csv](https://github.com/CBIIT/CTULIP/blob/main/example_results/predictions_18_ctoh_1-1.csv)   | Only the predicted primary tumor types and their probability scores. |
| [predictions_full_17_ctoh_1-1.csv](https://github.com/CBIIT/CTULIP/blob/main/example_results/predictions_full_17_ctoh_1-1.csv)  | [predictions_full_18_ctoh_1-1.csv](https://github.com/CBIIT/CTULIP/blob/main/example_results/predictions_full_17_ctoh_1-1.csv)    | The probabilitity scores for each tumor type for a full reference. |

## Software Setup

To set up the Python environment for TULIP:
1. Install [conda](https://docs.conda.io/en/latest/) package manager. 
2. Clone this repository. 

   ```bash
   git clone https://github.com/CBIIT/CTULIP.git
   ```

3. Create the environment as shown below.

   ```bash
   conda env create -f environment.yml -n ctulip
   conda activate ctulip
   ```

These commands generate an initial folder structure without a *models* folder. 

## Downloading Model Weights

To download the model weights needed for running TULIP:
1. Create an account on the Model and Data Clearinghouse [MoDaC](https://modac.cancer.gov). 
2. Run the following command.

   ```bash
   python modac_file_download.py
   ```
   
3. When prompted, enter your MoDaC credentials. The script creates a *models* folder and downloads files from the [CTULIP](https://modac.cancer.gov/assetDetails?dme_data_id=NCI-DME-MS01-21771349) asset in MoDaC to this folder. 



## Folder Structure

After performing the Software Setup steps and downloading the model weights, the following folder structure is available:

```
.
├── example_data/                 # folder containing example input files
├── example_results/              # folder containing example output files
├── gene_lists/                   # folder containing lists genes that is mapped 1-on-1 between humans and Canine
├── labels/                       # folder containing lists for 17 and 18 tumor types
├── models/                       # folder containing model weights, that will be downloaded running modac_file_download.py script
├── utils                         # Python helper scripts
├── environment.yml               # Python and libraries to run TULIP
├── modac_file_download.py        # Python script to download model weights from MoDaC
├── tulip.py                      # Python script of CTULIP
└── ...

```

## Data Setup

CTULIP accepts the RNA-seq data expressed as FPKM-UQ, in CSV, XLSX, and TSV file formats. 

Arrange the data with the Ensembl IDs in the first column and the expression values starting from the third column, as shown in the following example:

| gene_id | gene_name | CPT0019990006 | CPT0017440009 | CPT0077290006 |
| --------| ----------|----------|----------|----------|
| ENSG00000000003.13 | TSPAN6 | 157778.5731 | 76515.8557 | 205326.5947 |
| ENSG00000000005.5 | TNMD | 2828.4868 | 3321.4867 | 5517.4428 |
| ENSG00000000419.11 | DPM1 | 508866.9116 | 332778.5383 | 468852.2266 |

The [example_data](https://github.com/CBIIT/CTULIP/tree/main/example_data) folder provides example file. 

If the data file contains any duplicate Ensembl IDs, CTULIP removes them. Additionally, if the data file is missing any Ensembl IDs, CTULIP adds them and sets the expression values to 0 for each sample. 

## Running CTULIP

To run CTULIP, run a command in the following format:

   ```bash
   python ctulip.py -i <path/to/file> [options]
   ```

Consider the following example commands:

 * 17 tumor types
   ```bash
   python tulip.py -i example_data/test_data_all.csv -t 17 -o example_results/
   ```
 * 18 tumor types
   ```bash
   python tulip.py -i example_data/test_data_all.csv -t 18 -o example_results/
   ```
The results of the above example commands can be found in [example_results](https://github.com/CBIIT/CTULIP/blob/main/example_results/).

CTULIP parameters:

| Required? | Parameter | Description | Default |
| ------------- | ------------- | ------------- | ------------- |
| Yes  | -i, --input | The full path of the gene expression matrix file (FPKM-UQ) in the required format.  | (None) |
| No  | -t, --types | The number of tumor types, 17 or 18, to use for classification.  | 17 |
| No  | -o, --output_dir | The full path to the output directory. | (Current directory) |
| No  | -m, --min_score | The minimum probability score (0.0 to 1.0) for keeping the predicted primary tumor type. | (None) |

## Acknowledgments
This work has been supported in part by the Joint Design of Advanced Computing Solutions for Cancer (JDACS4C) program established by the U.S. Department of Energy (DOE) and the National Cancer Institute (NCI) of the National Institutes of Health.
