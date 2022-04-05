# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import os
from datetime import date, datetime
from os import getcwd, listdir, mkdir, path
from shutil import copy2, move
from Bio import SeqIO

import pandas as pd
from ruamel import yaml  # conda install -c conda-forge ruamel.yaml

# define location of sample sheet and workflow config

def update_sample_sheet():
    """
    This function updated the samples.csv file in your config with all .fasta files
    included in data/query and data/reference. Please make sure all query files are
    single fastas per sample. Reference fasta files can be multiple sequenence fasta 
    files.
    """

    config = snakemake.config

    QUERY_PATH = str(config["data-handling"]["data-query"])
    REFERENCE_PATH = str(config["data-handling"]["data-reference"])

    incoming_files = [f for f in listdir(QUERY_PATH)]

    if not incoming_files:
        print("No files in data/query")
        new_files_reference = pd.DataFrame(columns=["sample_name", "file", "type"])
    else:
        print("Updating sample sheet")
        # create dataframe
        new_files_query = pd.DataFrame(incoming_files, columns=["file"])

        # get id of sample by splitting the file handle
        new_files_query["sample_name"] = new_files_query["file"].apply(
            lambda x: (x.rsplit(".", 1)[0])
        )

        new_files_query["type"] = "query"
        # add path of file
        new_files_query["file"] = QUERY_PATH + "/" + new_files_query["file"]

        print(new_files_query)
        print("\t{} query samples added".format(len(new_files_query)))

    incoming_files = [f for f in listdir(REFERENCE_PATH)]

    if not incoming_files:
        print("No files in data/reference")
        new_files_reference = pd.DataFrame(columns=["sample_name", "file", "type"])
    else:
        print("Updating sample sheet")
        # create dataframe
        new_files_reference = pd.DataFrame(incoming_files, columns=["file"])

        # get id of sample by splitting the file handle
        new_files_reference["sample_name"] = new_files_reference["file"].apply(
            lambda x: (x.rsplit(".", 1)[0])
        )

        new_files_reference["type"] = "reference"
        # add path of file
        new_files_reference["file"] = REFERENCE_PATH + "/" + new_files_reference["file"]

        print(new_files_reference)
        print("\t{} reference file added".format(len(new_files_reference)))

        sample_sheet = pd.concat([new_files_query, new_files_reference])

        sample_sheet.to_csv(snakemake.input[0], index=False, columns=["sample_name", "file", "type"])


update_sample_sheet()
