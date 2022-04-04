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

# define location of sample sheet and workflow configx
SAMPLE_SHEET = snakemake.input[0]  # "config/pep/samples.csv"


def update_sample_sheet(SAMPLE_SHEET, verbose=True, dry_run=False):
    """
    This function
        - copies files from the incoming data directory to the snakemake data directory and
        - updates the sample sheet with the files copied to the snakemake data directory.

    The paths of these directory must be defined in the config yaml of the workflow as follows:

        data-handling:
            # path of incoming data
            incoming: [YOUR PATH]]
            # path to store data in the workflow
            data: [YOUR PATH]

    Args:
        SAMPLE_SHEET ([type]): Path to the location of the sample sheet
        CONFIG_YAML ([type]): Path to the location of the config yaml
        verbose (bool, optional): Provide additional details. Defaults to True.
    """

    config = snakemake.config

    DATA_PATH = str(config["data-handling"]["data"])
    IN_PATH = str(config["data-handling"]["incoming"])

    ##################################
    ### Check directories and data ###
    ##################################

    if verbose:
        print("Checking directories")

    # check if directory exist
    for given_path in [IN_PATH, DATA_PATH]:
        if not path.exists(given_path):
            raise Exception("Data directory (%s) not found" % given_path)

    # check if there is new data in the incoming data directory:
    # get files that are in incoming and do not contain 'ndetermined' and '.fastq.gz' in their name and are not under a specific filesize
    incoming_files = []
    for f in listdir(IN_PATH):
        if (
            path.isfile(path.join(IN_PATH, f))
            and ".fasta" in f
        ):
            incoming_files.append(f)

            # add subfolder for multi fasta in data path
            DATA_PATH += f.split(".fasta")
            if not path.isdir(DATA_PATH):
                mkdir(DATA_PATH)
            with open(f) as input_handle:
                for record in SeqIO.parse(input_handle, "fasta"):
                    print(record)


        else:
            print(f, "not used. Please check format!")



    # get files that are in outgoing directory
    data_files = [f for f in listdir(DATA_PATH) if path.isfile(path.join(DATA_PATH, f))]

    # print prompt, which data is in incoming and in outgoing and thus is not moved
    files_not_to_copy = [f for f in data_files if f in incoming_files]

    if files_not_to_copy:
        if verbose:
            print(
                "Following files are already located in %s and are not moved:"
                % DATA_PATH
            )
            i = 0
            for f in files_not_to_copy:
                print("\t%s" % DATA_PATH + f)
                i += 1
            print("\tIn total: {}".format(i))

    files_to_copy = [f for f in incoming_files if f not in data_files]

    ##################################
    ######### update the csv #########
    ##################################
    # check if there are no files to copy, thus list is empty
    if not files_to_copy:
        print("No (new) files to copy")
    else:

        if verbose:
            print("Updating sample sheet")
        # create dataframe
        new_files_df = pd.DataFrame(files_to_copy, columns=["file"])

        # get only files, that contain .fastq.gz
        new_files_df = new_files_df[new_files_df["file"].str.contains(".fastq.gz")]
        new_files_df = new_files_df[~new_files_df["file"].str.contains("Undetermined")]

        # get id of sample, thus split at first '_'
        new_files_df["sample_name"] = new_files_df["file"].apply(
            lambda x: (x.split("_", 1)[0])
        )

        # add path of file
        new_files_df["path"] = DATA_PATH + "/" + new_files_df["file"]

        # identify R1 or R2
        new_files_df["read"] = new_files_df["file"].apply(
            lambda x: "R1" if "R1" in x else "R2"
        )

        # set multiindex
        new_files_df.set_index(
            [new_files_df["sample_name"], new_files_df["read"]], inplace=True
        )

        # drop not need columns
        new_files_df.drop(columns=["file", "sample_name", "read"], inplace=True)

        # unstack multiindex
        new_files_df = new_files_df.unstack(1)
        new_files_df.sort_index(inplace=True)
        new_files_df.columns = ["fq1", "fq2"]
        new_files_df["date"] = today
        new_files_df["is_amplicon_data"] = 1
        new_files_df.loc[
            new_files_df.index.str.contains("No-RKI", case=False),
            ["include_in_high_genome_summary"],
        ] = "0"
        new_files_df.loc[
            ~new_files_df.index.str.contains("No-RKI", case=False),
            ["include_in_high_genome_summary"],
        ] = "1"
        print(new_files_df)

        new_sample_sheet = (
            pd.read_csv(SAMPLE_SHEET, index_col="sample_name")
            .append(new_files_df)
            .sort_values(by=["date", "sample_name"])
        )
        new_sample_sheet.index = new_sample_sheet.index.astype("str")

        # remove last line of sample.csv
        new_sample_sheet.drop("NAME", inplace=True, errors="ignore")

        # check for duplicates
        # TODO Generalize for more than two samples
        new_sample_sheet.index = new_sample_sheet.index.where(
            ~new_sample_sheet.index.duplicated(),
            new_sample_sheet.index.astype("str") + "_2",
        )
        # save to csv
        if verbose:
            print("\t{} samples added".format(len(new_files_df)))

        if not dry_run:
            new_sample_sheet.to_csv(snakemake.input[0])

        ##################################
        ## copying and moving the files ##
        ##################################

        if verbose:
            print("Copying files to " + DATA_PATH)

        # move all data in incoming path to data folder in snakemake
        # if ends with .fastq.gz and does not contain Undetermined
        i = 0
        for file in files_to_copy:
            if file.endswith(".fastq.gz") and not "ndetermined" in file:
                # if verbose:
                #     print("\t%s" % IN_PATH + file)
                if not dry_run:
                    copy2(IN_PATH + file, DATA_PATH)
                i += 1
        if verbose:
            print("\t{} files copied".format(i))

    # archiving incoming data
    all_incoming_files = [
        f for f in listdir(IN_PATH) if path.isfile(path.join(IN_PATH, f))
    ]
    if not all_incoming_files:
        print("No files to move")

    else:
        if verbose:
            print("Moving files to " + ARCHIVE_PATH + today)

        if not path.isdir(ARCHIVE_PATH + today):
            mkdir(ARCHIVE_PATH + today)

        archive_files = [
            f
            for f in listdir(ARCHIVE_PATH + today)
            if path.isfile(path.join(ARCHIVE_PATH + today, f))
        ]
        timestamp = datetime.now().strftime("%H:%M:%S")

        # move all files from incoming to archive
        i = 0
        for file in all_incoming_files:
            if not dry_run:
                if file in archive_files:
                    # if file is already in the archive add a timestemp when moving
                    move(
                        IN_PATH + file,
                        ARCHIVE_PATH + today + "/" + file + "-" + timestamp,
                    )
                else:
                    move(IN_PATH + file, ARCHIVE_PATH + today)
            i += 1
            pass
        if verbose:
            print("\t{} files moved".format(i))

        # save sample sheet in archive folder as backup
        if not dry_run:
            if files_to_copy:
                new_sample_sheet.to_csv(ARCHIVE_PATH + today + "/samples.csv")


update_sample_sheet(SAMPLE_SHEET)