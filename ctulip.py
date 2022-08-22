"""
CTULIP: Canine TUmor CLassIfication Predictor
This tool predicts the primary tumor type (17 or 32 types) based on RNA-seq Canine data.
This tool takes expression data of 14,761 genes, that is common to Canine-to-Human 1 to 1 mapped and Human protein coding genes  
Authors: Satish RG, Sara Jones 
"""
from __future__ import print_function
import os, warnings 
warnings.simplefilter(action="ignore", category=FutureWarning)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from keras.layers import Input, Dense, Dropout, Activation, Conv1D, MaxPooling1D, Flatten
from keras.models import Sequential, Model

#############
# Functions #
#############

def get_arg():
    """
    Define command-line arguments.
    """
    parse = argparse.ArgumentParser(description=__doc__, prog = "tulip.py", usage = "python %(prog)s -i <path/to/file>", add_help = True)
    reqargs = parse.add_argument_group("required")
    reqwdef = parse.add_argument_group("required with default")
    optargs = parse.add_argument_group("optional")
    reqargs.add_argument("-i", "--input", type=str, required=True, help="Provide the full path of the gene expression matrix file (FPKM-UQ) in the required format.")
    reqwdef.add_argument("-t", "--types", type=int, default=17, help="Indicate the number of tumor types, 17 or 18(inlcudes GBM class), you want to use for classification.")
    ##reqwdef.add_argument("-g", "--genes", type=str, default="pc", help="Indicate 'all' if you want to use all 60K genes or 'pc' for protein coding genes only.")
    optargs.add_argument("-o", "--output_dir", type=str, help="Provide the full path to the output directory.")
    optargs.add_argument("-m", "--min_score", type=float, help="Give the minimum probability score (0.0 to 1.0) for keeping the predicted primary tumor type.")
    all_args = parse.parse_args()
    return all_args

def check_path(filepath):
    """
    Check if filepath exists.
    Args: filepath (string) - path to file
    Returns: filepath
    """
    if os.path.exists(filepath):
        return filepath
    else:
        print("{} does not exist.".format(filepath))
        sys.exit()

def check_extension(filepath):
    """
    Check if input file has the following extensions (csv, xlsx, tsv) and read in file as dataframe.
    Args: filepath (string) - path to file
    Returns: df (dataframe)
    """
    if filepath.lower().endswith(".csv"):
        df = pd.read_csv(filepath)
        return df
    elif filepath.lower().endswith(".tsv"):
        df = pd.read_csv(filepath, header=0, sep="\t")
        return df
    elif filepath.lower().endswith(".xlsx"):
        df = pd.read_excel(filepath)
        return df
    else:
        print("File extension not accepted.")
        print("Only files with the following extensions are accepted: tsv, csv, and xlsx.")
        sys.exit()

def check_format(df):
    """
    Check if the file is in the correct format.
    Args: df (dataframe)
    """
    assert df[df.columns[0]][0].startswith("ENSG") == True, "File should contain Ensembl IDs in the first column."
    assert (df[df.columns[2:]].dtypes == "float64").all() == True, "Columns of expression values contain non-numeric data."

def check_genes(df, geneList):
    """
    Check that the list of genes matches the list of genes in the dataframe.
    If any missing genes, add rows with 0 values. 
    Args: df (dataframe)
          geneList (list)
    Returns: df (dataframe)
    """
    # check that gene ids are the same
    if df[df.columns[0]].tolist() == geneList:
        return df
    # sort gene ids and check that gene ids are the same
    elif sorted(df[df.columns[0]].tolist()) == geneList:
        # sort dataframe by gene ids
        return df.sort_values(by=df[df.columns[0]].name)
    else:
        # list of missing genes
        missing_list = list(set(geneList).difference(sorted(df[df.columns[0]].tolist())))
        print("Number of missing genes: {}".format(len(missing_list)))
        print("Adding missing genes with expression values set at 0.")
        for gene_id in missing_list:
            # make list with gene id
            row_list = [gene_id]
            # extend list with 0's
            row_list.extend([0] * (df.shape[1]-1))
            # add row to dataframe
            df.loc[len(df)] = row_list
        # drop any rows with duplicate gene ids
        df = df.drop_duplicates(subset=df[df.columns[0]].name, keep="first")
        # sort dataframe by gene ids
        return df.sort_values(by=df[df.columns[0]].name)

def prep_data(df):
    """
    Convert values to TPM and perform log10 normalization. Prepare dataset for classification.
    Args: df (dataframe)
    Returns: data (numpy array)
    """
    # convert to TPM
    tpm_df = df.div(df.sum(axis=1), axis=0)
    tpm_df = tpm_df * 1000000
    # change any values equal to or less than 0 to a very small number
    tpm_df[tpm_df <= 0] = 0.000001
    # log10 normalization
    tpm_df = tpm_df.astype(np.float64).apply(np.log10)
    # change any negative values to 0
    tpm_df[tpm_df < 0] = 0
    # expand to right dimensions
    data = np.expand_dims(tpm_df, axis=2)
    return data

def create_model(numClasses):
    """
    Creates model using number of genes as input and the number of tumor types as output.
    
    Parameters:
    numClasses (int): number of cancer types (17 or 18) as output
     (int): number of genes as input 
    
    Returns: model
    """
    filters = 128
    filter_len = 20
    stride = 1
    numGenes=14761

    model = Sequential()
    model.add(Conv1D(filters = filters,
                     kernel_size = filter_len,
                     strides = stride,
                     padding="valid",
                     input_shape=(numGenes, 1)))
    model.add(Activation("relu"))
    model.add(MaxPooling1D(pool_size = 1))
    model.add(Conv1D(filters = filters,
                     kernel_size = filter_len,
                     strides = stride,
                     padding = "valid"))
    model.add(Activation("relu"))
    model.add(MaxPooling1D(pool_size = 10))
    model.add(Flatten())
    model.add(Dense(200))
    model.add(Activation("relu"))
    model.add(Dropout(0.1))
    model.add(Dense(20))
    model.add(Activation("relu"))
    model.add(Dropout(0.1))
    model.add(Dense(numClasses))
    model.add(Activation("softmax"))
    
    return model

def load_trained_weights(numClasses, model, weightsDict):
    """
    Load trained weights to the model based on the number of genes and number of cancer types.
    
    Parameters:
    numClasses (int): number of cancer types (17 or 18) as output
    model: model
    weightsDict (dict): model type and the corresponding path to weights file
    
    Returns: model
    """
    if numClasses == 17:
        weights_path = weightsDict["pc_17_weights"]
    elif numClasses == 18:
        weights_path = weightsDict["pc_18_weights"]
    model.load_weights(weights_path)
    
    print("Loading weights from: {}".format(weights_path))
    
    return model

################
# Main program #
################  

def main():

    # paths to model weights
    cmg_17_weights_path = "models/17CT.model.h5"
    cmg_18_weights_path = "models/18CT.model.h5"

    weights_dict = {"pc_17_weights": cmg_17_weights_path,
                    "pc_18_weights": cmg_18_weights_path}

    # paths to gene lists
    all_gene_list_path = "gene_lists/caninegenes.txt"
    #pc_gene_list_path = "gene_lists/protein_coding_genes.txt"

    # paths to labels
    labels_17_path = "labels/17_tumors.csv"
    labels_18_path = "labels/18_tumors.csv"

    # read in arguments
    arg = get_arg()

    # check path to input file
    input = check_path(arg.input)
    # read in input file as dataframe
    input_df = check_extension(input)
    # check format of input
    check_format(input_df)

    gene_list = pd.read_csv(all_gene_list_path, header=None, sep="\t")[[0]][0].tolist()
    num_genes = len(gene_list)

    
    # if the number of genes in the input file is greater than the number of genes selected for classification, filter to genes
    if input_df.shape[0] > num_genes:
        input_df = input_df[input_df[input_df.columns[0]].isin(gene_list)]

    # check genes and transpose dataframe
    t_input_df = check_genes(input_df, gene_list).T
    # drop rows of gene ids and names and reset index
    t_input_df = (t_input_df.drop(t_input_df.index[0:2], axis=0)).reset_index()
    # rename index to samples
    t_input_df = t_input_df.rename(columns={t_input_df.columns[0]: "samples"})
    # sample names
    samples = t_input_df[["samples"]]
    # expression columns only
    features = t_input_df[t_input_df.columns[1:]]
    # prepare data for model
    input_data = prep_data(features)

    # check if minimum probability score is between 0 and 1
    min_score = None
    if arg.min_score:
        min_score = arg.min_score
        if (min_score > 1.0) or (min_score < 0):
            print("Minimum probability score must be between 0 and 1.")
            min_score = None

    # read in labels based on selection
    if arg.types == 17:
        labels_list = pd.read_csv(labels_17_path, header=None)[[1]][1].tolist()
        num_classes = 17
    elif arg.types == 18:
        labels_list = pd.read_csv(labels_18_path, header=None)[[1]][1].tolist()
        num_classes = 18
    else:
        print("{} is not an option.".format(arg.types))
        print("Indicate 17 or 18 for the number of tumor types.")
        sys.exit()

    if (num_classes == 17):
        print("Using model trained on intersection of  human protein coding and 1-1 mapped canine to human genes only and 17 tumor types.")
    elif (num_classes == 18):
        print("Using model trained on intersection of  human protein coding and 1-1 mapped canine to human genes only and 18 tumor types (inlcudes just GBM as extract group).")

    # create model
    selected_model = create_model(num_classes)
    selected_model.summary()

    # load weights
    selected_model = load_trained_weights(num_classes, selected_model, weights_dict)

    print("Getting predictions.")

    # get predictions
    predictions = selected_model.predict(input_data)

    # create dataframe of predictions and add tumor type labels
    predictions_df = pd.DataFrame(predictions, columns = labels_list)
    # add samples
    final_predictions_df = pd.concat([samples, predictions_df], axis=1)

    # get top predicted class and its probability
    final_predictions_df["predicted_class"] = final_predictions_df.iloc[:,1:num_classes].idxmax(axis = 1)
    # get full names for each predicted class
    final_predictions_df["probability"] = final_predictions_df.iloc[:,1:num_classes].max(axis = 1)

    if min_score:
        # if probability is less than minimum score, changed predicted class to None
        final_predictions_df["predicted_class"] = np.where(final_predictions_df["probability"] < min_score, "None", final_predictions_df["predicted_class"])

    print("Saving predictions.")

    gs="ctoh_1-1"

    if arg.output_dir:
        # if outdir directory doesn't exist, create new one
        if not os.path.exists(arg.output_dir):
            os.makedirs(arg.output_dir)
            print("A new directory has been created.")
        # includes probabilities for each tumor type
        final_predictions_df.to_csv(os.path.join(arg.output_dir, "predictions_full_{}_{}.csv".format(arg.types,gs)), header = True, index = False)
        # includes just predicted class and its probability
        final_predictions_df[["samples","predicted_class","probability"]].to_csv(os.path.join(arg.output_dir, "predictions_{}_{}.csv".format(arg.types, gs)), header = True, index = False)
    else:
        final_predictions_df.to_csv("predictions_full_{}_{}.csv".format(arg.types,gs), header = True, index = False)
        final_predictions_df[["samples","predicted_class","probability"]].to_csv("predictions_{}_{}.csv".format(arg.types,gs), header = True, index = False)

    print("Done.")

if __name__ == "__main__":
    main()
