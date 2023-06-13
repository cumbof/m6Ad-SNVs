#!/usr/bin/env python3
"""
Search for SNPs occurring outside of DRACH sites that destruction the DRACH RNA structure
from paired to unpaired or the other way around
"""

# TODO We are currently considering RMVar IDs for which SNPs do not occur on DRACH sites
#      We should also check whether a SNP occurring on a DRACH site could affect the structure of another DRACH site

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"
__version__ = "0.1.0"
__date__ = "Jun 13, 2023"

import argparse as ap
import copy
import errno
import itertools
import os
import re
import requests
import subprocess
import sys
import tempfile
import wget
import zipfile
from io import StringIO

from Bio import SeqIO
from tqdm import tqdm

# Tool name
TOOL_ID = "m6A"

# List of external software dependencies
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/
# https://github.com/ViennaRNA/ViennaRNA
DEPENDENCIES = [
    "datasets",
    "RNAfold"
]

# URL to the Human m6A table
RMVAT_HUMAN_M6A_TABLE = "https://rmvar.renlab.org/download/RMVar_Human_basic_info_m6A.txt"

# Permalinks to the forna v1.0.1 assets CSS and JS
FORNA_ASSETS = [
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.css",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.css.map",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.js",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.js.map"
]


def read_params(argv):
    """
    Read and test input arguments

    :param argv:    List of arguments
    :return:        The ArgumentParser object
    """

    p = ap.ArgumentParser(
        prog=TOOL_ID,
        description=(
            "Search for SNPs occurring outside of DRACH sites that destruction the DRACH RNA structure "
            "from paired to unpaired or the other way around"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--table",
        type=os.path.abspath,
        required = "--download" not in argv,
        help="Path to the RMVar table with m6A modifications",
    )
    p.add_argument(
        "--download",
        action="store_true",
        default=False,
        help="Download the RMVar table into the working directory",
    )
    p.add_argument(
        "--co-occurrence",
        action="store_true",
        default=False,
        dest="co_occurrence",
        help="Merge m6A modifications on the same reference sequence",
    )
    p.add_argument(
        "--expand-left",
        type=int,
        default=0,
        dest="expand_left",
        help="Expand the reference sequence of N bases to the left",
    )
    p.add_argument(
        "--expand-right",
        type=int,
        default=0,
        dest="expand_right",
        help="Expand the reference sequence of N bases to the right",
    )
    p.add_argument(
        "--paired-unpaired",
        action="store_true",
        default=False,
        dest="paired_unpaired",
        help=(
            "Report a RMVar ID if at least one paired DRACH site in the reference structure "
            "become unpaired in the alterative structure"
        ),
    )
    p.add_argument(
        "--unpaired-paired",
        action="store_true",
        default=False,
        dest="unpaired_paired",
        help=(
            "Report a RMVar ID if at least one unpaired DRACH site in the reference structure "
            "become paired in the alterative structure"
        ),
    )
    p.add_argument(
        "--out-table",
        type=os.path.abspath,
        dest="out_table",
        help="Path to the output file (optional, it overwrites the output file if it exists)",
    )
    p.add_argument(
        "--out-structures",
        type=os.path.abspath,
        dest="out_structures",
        help=(
            "Path to the output folder with RNAfold structures (optional, it overwrites the output files if they exist). "
            "In order to visualize the secondary structures of both the reference and alterative sequences, "
            "go to http://rna.tbi.univie.ac.at/forna/, click on \"Add Sequence and Secondary Structure\" and "
            "paste the content of a .fasta file, then click on \"Colors -> Custom\" and paste the content of "
            "the corresponding .color file. PNG figures are also produced."
        ),
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version="{} version {} ({})".format(TOOL_ID, __version__, __date__),
        help="Print the {} tool version and exit".format(TOOL_ID),
    )
    return p.parse_args(argv)


def format_m6A(rm_id=None, data=None, header=None):
    """
    Format m6A info

    :param rm_id:   RMVar ID
    :param data:    Dictionary with m6A info
    :param header:  List with header fields
    :return:        m6A info
    """

    if header and not data:
        return "rm_id\t{}".format("\t".join(header))

    if rm_id and header and data:
        rm_info = list()

        for info in header:
            if isinstance(data[info], list):
                rm_info.append(",".join(data[info]))

            elif isinstance(data[info], dict):
                rm_info.append(",".join(["{}:{}".format(key, data[info][key]) for key in data[info]]))

            else:
                rm_info.append(str(data[info]))

        return "{}\t{}".format(rm_id, "\t".join(rm_info))


def forna_index(init=False, header=None, row=None, close=False, filepath=None):
    """
    Create a table with target m6A info and links to the HTML pages with structures

    :param init:        Initialize the table
    :param header:      List with header info
    :param row:         Dictionary with target m6A information
    :param close:       Close the table
    :param filepath:    Path to the index HTML file
    """

    if header:
        header = copy.deepcopy(header)

        # Remove info
        for info in ("reference_sequence", "reference_base", "alterative_sequence", "alterative_base"):
            header.remove(info)

    if init:
        with open(filepath, "w+") as index_file:
            index_file.write(
                """
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>Target m6A sites</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
        <link href="https://cdn.datatables.net/v/dt/dt-1.13.4/datatables.min.css" rel="stylesheet"/>
        <script src="https://cdn.datatables.net/v/dt/dt-1.13.4/datatables.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
    </head>
    <body>
        <div class="container">
            <h2>Target m6A index</h2>
            <p>Click on a RMVar ID to inspect reference and alterative structures</p>
            <table class="table table-condensed table-bordered table-striped" id="myTable">
                """
            )

            if header:
                index_file.write("<thead><tr><th>{}</th></tr></thead>\n".format("</th><th>".join(header)))

            index_file.write("<tbody>\n")

    if row and header:
        row = copy.deepcopy(row)

        # Remove info
        for info in ("reference_sequence", "reference_base", "alterative_sequence", "alterative_base"):
            del row[info]

        if "rm_id" in row:
            row["rm_id"] = '<a href="data/{rm_id}.html" target="_blank">{rm_id}</a>'.format(rm_id=row["rm_id"])

        with open(filepath, "a+") as index_file:
            index_file.write("<tr>{}</tr>".format("".join(["<td>{}</td>".format(row[info]) for info in header])))

    if close:
        with open(filepath, "a+") as index_file:
            index_file.write(
                """
                </tbody>
            </table>
        </div>

        <script>
            $(document).ready(function(){
                $("#myTable").dataTable({
                    scrollX: true
                });
            });
        </script>
    </body>
</html>
                """
            )


def forna_template(data=None) -> str:
    """
    Return the forna HTML template to render structures

    :param data:    Dictionary with m6A info
    :return:        The HTML template
    """

    return """
<!DOCTYPE html>
<html lang="en">
    <head>
        <meta charset="utf-8">
        <title>RNAPlot - {rm_id}</title>
        <link rel="stylesheet" href="../assets/fornac.css" />
        <style>
            svg {{
                width: 100%;
                height: 100%;
                max-width: 400px;
                border: 1px solid gray;
            }}
        </style>
    </head>
    <body>
        <div style="max-width: 800px; margin: auto;">
            <h2>{rm_id}</h2>
            <div id="reference">
                <h3>Reference</h3>
                <p>Sequence: {reference_sequence}</p>
                <p>Structure: {reference_structure}</p>
                <p>Base: {reference_base}</p>
            </div>
            <div id="alterative">
                <h3>Alterative</h3>
                <p>Sequence: {alterative_sequence}</p>
                <p>Structure: {alterative_structure}</p>
                <p>Base: {alterative_base}</p>
            </div>

            <script src="https://unpkg.com/jquery"></script>
            <script src="https://unpkg.com/d3@3.5"></script>
            <script src="https://unpkg.com/d3-grid"></script>
            <script src="../assets/fornac.js"></script>

            <script type="text/javascript">
                var reference_structure = "{reference_structure}";
                var reference_sequence = "{reference_sequence}";
                var reference_colors = "{reference_colors}";

                var reference_data = {{
                    "structure": reference_structure,
                    "sequence": reference_sequence,
                    "name": "Reference"
                }};

                let reference_container = new fornac.FornaContainer(
                    "#reference",
                    {{"animation": false, "zoomable": true, "initialSize": [300,300]}}
                );

                reference_container.addRNA(reference_data.structure, reference_data);
                reference_container.addCustomColorsText(reference_colors);

                var alterative_structure = "{alterative_structure}";
                var alterative_sequence = "{alterative_sequence}";
                var alterative_colors = "{alterative_colors}";

                var alterative_data = {{
                    "structure": alterative_structure,
                    "sequence": alterative_sequence,
                    "name": "Alterative"
                }};

                let alterative_container = new fornac.FornaContainer(
                    "#alterative",
                    {{"animation": false, "zoomable": true, "initialSize": [300,300]}}
                );

                alterative_container.addRNA(alterative_data.structure, alterative_data);
                alterative_container.addCustomColorsText(alterative_colors);
            </script>
        </div>
    </body>
</html>
    """


def main():
    # Load command line parameters
    args = read_params(sys.argv[1:])

    if args.download:
        print("Downloading RMVat m6A table: {}".format(RMVAT_HUMAN_M6A_TABLE))

        try:
            m6A_table_filepath = os.path.join(os.getcwd(), os.path.basename(RMVAT_HUMAN_M6A_TABLE))

            if not os.path.isfile(m6A_table_filepath):
                wget.download(RMVAT_HUMAN_M6A_TABLE, os.getcwd())

            args.table = m6A_table_filepath

        except:
            raise Exception("Unable to retrieve data from the RMVat database")

    else:
        if args.table and not os.path.isfile(args.table):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.table)

    if args.out_structures:
        # Path to the folder with .fasta and .color files
        out_data_folder = os.path.join(args.out_structures, "data")
        os.makedirs(out_data_folder, exist_ok=True)

        # Path to the folder with "forna" assets
        forna_folder = os.path.join(args.out_structures, "assets")
        os.makedirs(forna_folder, exist_ok=True)

        # Retrieve "forna" CSS and JS
        for asset in FORNA_ASSETS:
            try:
                asset_filepath = os.path.join(args.out_structures, os.path.basename(asset))

                if not os.path.isfile(asset_filepath):
                    wget.download(asset, out=os.path.abspath(forna_folder))

            except:
                # Do not block the execution here since the host could not have access to internet
                # We can still produce the .fasta and .color files
                print("Unable to download {}".format(asset))

    m6A_data = dict()

    modifications_count = 0

    print("\nLoading MRVat m6A table")

    with open(args.table) as m6A_table:
        header = m6A_table.readline().strip().split("\t")

        for line in m6A_table:
            if line.strip():
                line_split = [value.strip() for value in line.split("\t")]

                if len(line_split) < len(header):
                    missing_info = ["-"] * (len(line_split) - len(header))
                    line_split.extend(missing_info)

                try:
                    if line_split[header.index("modification_type")] != "m6A" or line_split[header.index("species")] != "Human":
                        raise Exception("Is this the right table?")

                    process_m6A = True

                    # Get m6A modification occurring on exons
                    if "exon" not in line_split[header.index("gene_region")]:
                        process_m6A = False

                    # Get m6A modification for which we have a disease association
                    #if not line_split[header.index("tumor")] or line_split[header.index("tumor")] == "-":
                    #    process_m6A = False

                    if process_m6A:
                        reference_sequence = line_split[header.index("reference_sequence")]
                        alterative_sequence = line_split[header.index("alterative_sequence")]

                        if len(reference_sequence) + args.expand_left + args.expand_right > 1023:
                            raise Exception(
                                "Sequence length cannot exceed 1023 bases due to ViennaRNA:RNAfold limitation. "
                                "Reference and alterative sequence length in RMVar table is {}. "
                                "Please reduce --expand-left and --expand-right accordingly".format(len(reference_sequence))
                            )

                        # Search for DRACH sites
                        # DRACH regex involving 5 positions only: [AGU][AG]AC[ACU]
                        regex = re.compile("[AGT][AG]AC[ACT]")

                        # Focus on SNPs occurring outside the DRACH sites
                        snp_in_DRACH = False

                        for match in regex.finditer(reference_sequence):
                            # Check whether the alterative base occurs in a DRACH site
                            if reference_sequence[match.start():match.start() + 5] != alterative_sequence[match.start():match.start() + 5]:
                                snp_in_DRACH = True

                        if not snp_in_DRACH:
                            if reference_sequence not in m6A_data:
                                m6A_data[reference_sequence] = dict()

                            m6A_data[reference_sequence][line_split[header.index("rm_id")]] = {
                                    "rm_id": line_split[header.index("rm_id")],
                                    "chromosome": line_split[header.index("chromosome")],
                                    "strand": line_split[header.index("strand")],
                                    "position": line_split[header.index("position")],
                                    "snp_start": line_split[header.index("snp_start")],
                                    "snp_end": line_split[header.index("snp_end")],
                                    "reference_sequence": reference_sequence,
                                    "alterative_sequence": alterative_sequence,
                                    "reference_base": line_split[header.index("reference_base")],
                                    "alterative_base": line_split[header.index("alterative_base")],
                                    "gene": line_split[header.index("gene")].split(","),
                                    "gene_database_id": {
                                        value.split(":")[0]: value.split(":")[1] for value in line_split[header.index("gene_database_id")].split(",")
                                    },
                                    "gene_type": line_split[header.index("gene_type")].split(","),
                                    "gene_region": line_split[header.index("gene_region")].split(","),
                                    "tumor": line_split[header.index("tumor")].split(","),
                                    "structure_state": "",
                                    "drach_sites": 0
                                }

                            modifications_count += 1

                except:
                    raise Exception("Malformed table")

    print("Considering {} m6A modifications over {} sequences".format(modifications_count, len(m6A_data)))

    if m6A_data:
        # Take track of RNAfold results
        rnafold_structures = dict()

        first_ref = list(m6A_data.keys())[0]
        first_mod = list(m6A_data[first_ref].keys())[0]
        header_info = list(m6A_data[first_ref][first_mod].keys())

        if args.out_table:
            with open(args.out_table, "w+") as outfile:
                outfile.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        format_m6A(header=header_info),
                        "rnafold_reference_mfe_structure",
                        "rnafold_reference_minimum_free_energy",
                        "rnafold_alterative_mfe_structure",
                        "rnafold_alterative_minimum_free_energy",
                        "relative_drach_positions"
                    )
                )

        if args.out_structures:
            out_data_index = os.path.join(args.out_structures, "index.html")

            # Init target index page
            forna_index(init=True, header=header_info, filepath=out_data_index)

        for reference_sequence in tqdm(m6A_data, desc="Sequences"):
            has_been_expanded = False
            expanded_left = ""
            expanded_right = ""

            if args.expand_left == 0 and args.expand_right == 0:
                gene_reference_sequence = reference_sequence

            else:
                # Get the gene symbols
                gene_symbols = list(set([gene for rm_id in m6A_data[reference_sequence] for gene in m6A_data[reference_sequence][rm_id]["gene"] if gene.strip()]))

                if not gene_symbols:
                    # Cannot retrieve the gene sequence without the gene symbol
                    gene_reference_sequence = reference_sequence

                else:
                    # Retrieve the gene sequence from NCBI
                    # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/gene/datasets_download_gene_symbol/
                    # Retrieve the sequences for all the involved genes with a single query
                    # Enlarge the reference and alterative sequences to N bases
                    with tempfile.NamedTemporaryFile() as genes_list, tempfile.NamedTemporaryFile() as genes_info:
                        with open(genes_list.name, "wt") as genes_list_file:
                            for gene in gene_symbols:
                                genes_list_file.write("{}\n".format(gene))

                        try:
                            # Retrieve the gene sequence from NCBI Datasets
                            subprocess.check_call(
                                [
                                    "datasets",
                                    "download",
                                    "gene",
                                    "symbol",
                                    "--taxon",
                                    "human",
                                    "--include",
                                    "gene",
                                    "--no-progressbar",
                                    "--inputfile",
                                    genes_list.name,
                                    "--filename",
                                    genes_info.name
                                ],
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL
                            )

                            retrieved = True

                        except:
                            # Unable to retrieve data from NCBI Datasets
                            print("Warning: genes \"{}\" not found in NCBI Datasets".format(gene_symbols))

                            retrieved = False

                        if retrieved:
                            try:
                                with zipfile.ZipFile(genes_info.name, "r") as archive:
                                    gene_stringio = StringIO(archive.read("ncbi_dataset/data/gene.fna").decode("utf-8"))

                                    # Load the gene sequence
                                    gene_content = SeqIO.parse(gene_stringio, "fasta")

                                    # Search for the reference sequence into the whole gene sequence
                                    regex = re.compile(reference_sequence.replace("_", ""))

                                    match_found = False

                                    for record in gene_content:
                                        for match in regex.finditer(str(record.seq)):
                                            # Take track of the expanded sub sequences
                                            expanded_left = str(record.seq)[match.start() - args.expand_left : match.start()]
                                            expanded_right = str(record.seq)[match.end() : match.end() + args.expand_right]

                                            # Expand the reference sequence
                                            gene_reference_sequence = "{}{}{}".format(expanded_left, reference_sequence, expanded_right)

                                            match_found = True

                                            has_been_expanded = True

                                            # There should be just one match
                                            # If it gets here, the reference sequence exists
                                            break

                                        if match_found:
                                            break

                                    if not match_found:
                                        # This should never happen
                                        # There is something wrong with the sequence retrieved from NCBI Datasets
                                        # Do not expand the sequence and use the reference sequence in RMVar
                                        gene_reference_sequence = reference_sequence

                                    gene_stringio.close()

                                not_found = False
                            
                            except:
                                # Unable to read the zip package
                                print("Warning: sequence not found in archive for genes \"{}\"".format(gene_symbols))

                                not_found = True

                        if not retrieved or not_found:
                            # It was unable to retrieve the gene sequence
                            # Usually because of "Error: No genes found that match selection"
                            gene_reference_sequence = reference_sequence

            m6A_entries = m6A_data[reference_sequence]

            if args.co_occurrence and len(m6A_entries) > 1:
                # There could be multiple alterative sequences for a single reference sequence
                # Combine all possible combination of modifications occurring on the same reference sequence
                # This is going to be combinatorial but the maximum number of modifications per gene is very limited
                rm_ids = list(m6A_data[reference_sequence].keys())

                for how_many_entries in range(2, len(rm_ids) + 1):
                    for subset_entries in list(itertools.combinations(rm_ids, how_many_entries)):
                        combined_entry = dict()

                        for rm_id in subset_entries:
                            m6A_sub_entry = m6A_data[reference_sequence][rm_id]

                            if not combined_entry:
                                combined_entry = m6A_sub_entry

                            else:
                                combined_entry["rm_id"] = "{},{}".format(combined_entry["rm_id"], m6A_sub_entry["rm_id"])
                                combined_entry["snp_start"] = "{},{}".format(combined_entry["snp_start"], m6A_sub_entry["snp_start"])
                                combined_entry["snp_end"] = "{},{}".format(combined_entry["snp_end"], m6A_sub_entry["snp_end"])

                                for alt_pos in range(len(combined_entry["alterative_sequence"])):
                                    if m6A_sub_entry["alterative_sequence"][alt_pos] != combined_entry["alterative_sequence"][alt_pos]:
                                        string_list = list(combined_entry["alterative_sequence"])

                                        # Update the alterative sequence
                                        string_list[alt_pos] = m6A_sub_entry["alterative_sequence"][alt_pos]

                                        combined_entry["alterative_sequence"] = "".join(string_list)

                                combined_entry["reference_base"] = "{},{}".format(combined_entry["reference_base"], m6A_sub_entry["reference_base"])
                                combined_entry["alterative_base"] = "{},{}".format(combined_entry["alterative_base"], m6A_sub_entry["alterative_base"])

                                combined_entry["tumor"].extend(m6A_sub_entry["tumor"])
                                combined_entry["tumor"] = list(set(combined_entry["tumor"]))

                        m6A_entries[combined_entry["rm_id"]] = combined_entry

            for rm_id in m6A_entries:
                m6A_entry = m6A_entries[rm_id]

                alterative_sequence = m6A_entry["alterative_sequence"]

                if len(alterative_sequence) != len(gene_reference_sequence):
                    # Expand the alterative sequence
                    if args.expand_left > 0:
                        alterative_sequence = "{}{}".format(expanded_left, alterative_sequence)

                    if args.expand_right > 0:
                        alterative_sequence = "{}{}".format(alterative_sequence, expanded_right)

                # Override original reference and alterative sequences with the expanded ones
                m6A_entry["reference_sequence"] = gene_reference_sequence
                m6A_entry["alterative_sequence"] = alterative_sequence

                rnafold_structures[rm_id] = dict()

                with tempfile.NamedTemporaryFile() as rnafold_in, tempfile.NamedTemporaryFile() as rnafold_out:
                    with open(rnafold_in.name, "wt") as rnafold_in_file:
                        # Write reference sequence
                        rnafold_in_file.write(">reference_sequence\n{}\n".format(gene_reference_sequence))
                        # Write alterative sequence
                        rnafold_in_file.write(">alterative_sequence\n{}\n".format(alterative_sequence))

                    try:
                        # Run ViennaRNA:RNAfold on reference and alterative sequence
                        subprocess.check_call(
                            [
                                "RNAfold",
                                "-i",
                                rnafold_in.name
                            ],
                            stdout=rnafold_out,
                            stderr=subprocess.DEVNULL
                        )

                    except subprocess.CalledProcessError as e:
                        raise Exception("An error has occurred while processing {} with RNAfold".format(rm_id)).with_traceback(e.__traceback__)

                    # Get mfe structures in bracket notation
                    with open(rnafold_out.name, "rt") as rnafold_out_file:
                        last_seq_id = None

                        for line in rnafold_out_file:
                            line = line.strip()

                            if line:
                                if line.startswith(">"):
                                    last_seq_id = line[1:]

                                if line.startswith(".") or line.startswith("("):
                                    line_split = line.split(" ")

                                    rnafold_structures[rm_id][last_seq_id] = {
                                        "mfe_structure": line_split[0],
                                        "minimum_free_energy": line_split[-1][1:-1]  # kcal/mol
                                    }

                if rnafold_structures[rm_id]:
                    reference_structure = rnafold_structures[rm_id]["reference_sequence"]["mfe_structure"]
                    alterative_structure = rnafold_structures[rm_id]["alterative_sequence"]["mfe_structure"]

                    if reference_structure != alterative_structure:
                        # Folded structures of reference and alterative sequences are different here for sure
                        # Check whether the alterative structure is unpaired on the DRACH site
                        regex = re.compile("[AGT][AG]AC[ACT]")

                        # Focus on SNPs occurring outside DRACH sites
                        at_least_one_DRACH_unpaired = False
                        at_least_one_DRACH_paired = False

                        # Take track of the positions of the DRACH sites in the reference sequence
                        # Since the SNPs do not occur on DRACH sites, reference and alterative sequences
                        # have the same DRACH sites in the same positions
                        DRACH_positions = list()

                        for match in regex.finditer(gene_reference_sequence):
                            # We just need the start position since it is always 5 bases long
                            DRACH_positions.append(match.start() + 1)

                            # Check whether the structure at the DRACH site on changed from paired to unpaired or the other way around
                            # A dot in the dot-bracket notation represents an unpaired position
                            how_many_dots_reference = reference_structure[match.start():match.start() + 5].count(".")
                            how_many_dots_alterative = alterative_structure[match.start():match.start() + 5].count(".")

                            if args.paired_unpaired:
                                if how_many_dots_reference == 0 and how_many_dots_alterative > 0:
                                    at_least_one_DRACH_unpaired = True
                            
                            if args.unpaired_paired:
                                if how_many_dots_reference > 0 and how_many_dots_alterative == 0:
                                    at_least_one_DRACH_paired = True

                        if at_least_one_DRACH_unpaired or at_least_one_DRACH_paired:
                            m6A_entry["structure_state"] = "paired-unpaired" if at_least_one_DRACH_unpaired else "unpaired-paired"
                            m6A_entry["drach_sites"] = len(DRACH_positions)

                            if args.out_table:
                                with open(args.out_table, "a+") as outfile:
                                    outfile.write(
                                        "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                            format_m6A(rm_id=rm_id, header=header_info, data=m6A_entry),
                                            reference_structure,
                                            rnafold_structures[rm_id]["reference_sequence"]["minimum_free_energy"],
                                            alterative_structure,
                                            rnafold_structures[rm_id]["alterative_sequence"]["minimum_free_energy"],
                                            ",".join([str(drach_pos) for drach_pos in DRACH_positions])
                                        )
                                    )

                            if args.out_structures:
                                # Also take track of different bases positions
                                snp_positions = list()

                                for pos in range(len(gene_reference_sequence)):
                                    if gene_reference_sequence[pos].lower() != alterative_sequence[pos].lower():
                                        snp_positions.append(pos + 1)

                                fasta_filepath = os.path.join(out_data_folder, "{}.fasta".format(rm_id))
                                color_filepath = os.path.join(out_data_folder, "{}.color".format(rm_id))
                                forna_filepath = os.path.join(out_data_folder, "{}.html".format(rm_id))

                                with open(fasta_filepath, "w+") as fasta_file, open(color_filepath, "w+") as color_file, open(forna_filepath, "w+") as forna_file:
                                    # Define DRACH sites colors
                                    drach_colors = " ".join(["{}-{}:green".format(drach_pos, drach_pos + 4) for drach_pos in DRACH_positions])
                                    # Define SNPs colors
                                    snp_colors = " ".join(["{}:orange".format(snp_pos) for snp_pos in snp_positions])

                                    # Highlight DRACH sites and SNPs in both reference and alterative structures
                                    color_file.write("{} {}\n".format(drach_colors, snp_colors))

                                    # Write reference sequence and structure
                                    fasta_file.write(
                                        ">{}-reference\n{}\n{}\n".format(
                                            rm_id,
                                            gene_reference_sequence,
                                            reference_structure
                                        )
                                    )

                                    # Write alterative sequence and structure
                                    fasta_file.write(
                                        ">{}-alterative\n{}\n{}\n".format(
                                            rm_id,
                                            alterative_sequence,
                                            alterative_structure
                                        )
                                    )

                                    # Also define the html page
                                    forna_file.write(
                                        forna_template(data=m6A_entry).format(
                                            rm_id=rm_id,
                                            reference_structure=reference_structure,
                                            reference_sequence=gene_reference_sequence,
                                            reference_base=m6A_entry["reference_base"],
                                            reference_colors="{} {}".format(drach_colors, snp_colors),
                                            alterative_structure=alterative_structure,
                                            alterative_sequence=alterative_sequence,
                                            alterative_base=m6A_entry["alterative_base"],
                                            alterative_colors="{} {}".format(drach_colors, snp_colors)
                                        )
                                    )

                                # Add row to forna index page
                                forna_index(header=header_info, row=m6A_entry, filepath=out_data_index)

        if args.out_structures:
            # Finalize forna index page
            forna_index(close=True, filepath=out_data_index)


if __name__ == "__main__":
    main()
