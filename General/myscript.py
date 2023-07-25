#!/usr/bin/env python3
"""
Describe your script here
"""

__author__ = "Fabio Cumbo (fabio.cumbo@gmail.com)"

__version__ = "0.1.0"
__date__ = "Jul 7, 2023"

import argparse as ap
import copy
import errno
import multiprocessing as mp
import os
import re
import subprocess
import sys
import tempfile
import tqdm
import wget
from functools import partial

import pandas as pd
import pysam
from Bio import Seq

# Tool name
TOOL_ID = "myscript"

# List of external software dependencies
# https://github.com/ViennaRNA/ViennaRNA
DEPENDENCIES = [
    "RNAfold"
]

# Permalinks to the forna v1.0.1 assets CSS and JS
FORNA_ASSETS = [
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.css",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.css.map",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.js",
    "https://raw.githubusercontent.com/ViennaRNA/fornac/36df3c5d73d2f651c3c3b5266e7d705e5bb1d3d1/dist/fornac.js.map"
]

# Table columns
TABLE_COLS = [
    "entry_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "gene_symbol",
    "transcript_id",
    "reference_sequence",
    "reference_structure",  # MFE structure in dot-bracket notation as reported by RNAfold
    "reference_minimum_free_energy",  # In kcal/mol as reported by RNAfold
    "reference_pattern_sites",  # Number of sites in the reference (and alternate) sequence (e.g., DRACH)
    "reference_pattern_positions",
    "alternate_position",
    "alternate_sequence",
    "alternate_structure",  # MFE structure in dot-bracket notation as reported by RNAfold
    "alternate_minimum_free_energy",  # In kcal/mol as reported by RNAfold
    "alternate_pattern_sites",  # Number of sites in the reference (and alternate) sequence (e.g., DRACH)
    "alternate_pattern_positions",
    "structure_state",  # Paired-to-unpaired or unpaired-to-paired
    "delta_g",
    "synonymous",  # In case of single-nucleotide variations
    "frameshift",  # In case of indels
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
            "Describe your script here"
        ),
        formatter_class=ap.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--bed",
        type=os.path.abspath,
        required=True,
        help="Path to the BED file with the coordinates of the genomic regions of interest",
    )
    p.add_argument(
        "--bed-skip-lines",
        type=int,
        default=1,
        dest="bed_skip_lines",
        help="Skip the first number of lines in the input BED file",
    )
    p.add_argument(
        "--vcf",
        type=os.path.abspath,
        required=True,
        help="Path to the VCF file with the SNPs coordinates",
    )
    p.add_argument(
        "--genome",
        type=os.path.abspath,
        required=True,
        help="Path to the fasta genome file",
    )
    p.add_argument(
        "--search",
        type=str,
        default="[AGT][AG]AC[ACT]",
        help=(
            "Search for a specific pattern in the sequences. "
            "It accepts regular expressions. It searches for DRACH sites by default"
        ),
    )
    p.add_argument(
        "--out-table",
        type=os.path.abspath,
        dest="out_table",
        required = "--out-html" not in argv,
        help="Path to the output file (optional, it overwrites the output file if it exists)",
    )
    p.add_argument(
        "--out-html",
        type=os.path.abspath,
        dest="out_html",
        required = "--out-table" not in argv,
        help="Path to the output folder with RNAfold structures (optional, it overwrites the output files if they exist)",
    )
    p.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Make it parallel",
    )
    p.add_argument(
        "-v",
        "--version",
        action="version",
        version="{} version {} ({})".format(TOOL_ID, __version__, __date__),
        help="Print the {} tool version and exit".format(TOOL_ID),
    )
    return p.parse_args(argv)


def regex_is_valid(regex: str) -> bool:
    """
    Check if the input string is a valid regular expression
    It must be a fixed-length regex

    :param regex:   The input regular expression
    :return:        True if regex is a valid regular expression 
    """

    try:
        re.compile(regex)

        # Check whether it is a fixed-length regex
        if "+" in regex or "*" in regex:
            return False

        return True

    except re.error:
        return False


def forna_index(init=False, header=None, row=None, close=False, filepath=None):
    """
    Create a table with target info and links to the HTML pages with structures

    :param init:        Initialize the table
    :param header:      List with header info
    :param row:         List with target info
    :param close:       Close the table
    :param filepath:    Path to the index HTML file
    """

    text = ""

    remove_positions = list()

    if header:
        header = copy.deepcopy(header)

        remove_info = [
            "reference_sequence",
            "reference_structure",
            "reference_pattern_positions",
            "alternate_sequence",
            "alternate_structure",
            "alternate_pattern_positions",
            "frameshift"
        ]

        remove_positions = [header.index(info) for info in remove_info]

        # Remove info
        for info in remove_info:
            header.remove(info)

    if init:
        text += """
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>Target sites</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js"></script>
        <link href="https://cdn.datatables.net/v/dt/dt-1.13.4/datatables.min.css" rel="stylesheet"/>
        <script src="https://cdn.datatables.net/v/dt/dt-1.13.4/datatables.min.js"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
        <style>
            tfoot input {
                width: 100%;
                padding: 3px;
                box-sizing: border-box;
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h2>Target index</h2>
            <p>Click on the IDs to inspect reference and alternate structures</p>
            <table class="table table-condensed table-bordered table-striped display" id="myTable" style="width: 100%">
                    """

        if header:
            text += "<thead><tr><th>{}</th></tr></thead>\n".format("</th><th>".join(header))
        
        text += "<tbody>\n"

        if filepath:
            with open(filepath, "w+") as index_file:
                index_file.write(text)

    if row and header:
        row = copy.deepcopy(row)

        # Entry ID is always in the first position
        row[0] = '<a href="data/{entry_id}.html" target="_blank">{entry_id}</a>'.format(entry_id=row[0])

        text += "<tr>{}</tr>".format("".join(["<td>{}</td>".format(row[pos]) for pos in range(0, len(row)) if pos not in remove_positions]))

        if filepath:
            with open(filepath, "a+") as index_file:
                index_file.write(text)

    if close:
        text += "</tbody>"
        
        if header:
            text += "<tfoot><tr><th>{}</th></tr></tfoot>\n".format("</th><th>".join(header))

        text += """
            </table>
        </div>

        <script>
            $(document).ready(function(){
                $("#myTable").DataTable({
                    scrollX: true,
                    initComplete: function () {
                        this.api()
                            .columns()
                            .every(function () {
                                var column = this;
                                var title = column.footer().textContent;
                                // Create input element and add event listener
                                $('<input type="text" placeholder="Search ' + title + '" />')
                                    .appendTo($(column.footer()).empty())
                                    .on('keyup change clear', function () {
                                        if (column.search() !== this.value) {
                                            column.search(this.value).draw();
                                        }
                                    });
                            });
                    },
                });
            });
        </script>
    </body>
</html>
                """

        if filepath:
            with open(filepath, "a+") as index_file:
                index_file.write(text)
    
    if not filepath:
        return text

    return None


def forna_template() -> str:
    """
    Return the forna HTML template to render structures

    :return:    The HTML template
    """

    return """
<!DOCTYPE html>
<html lang="en">
    <head>
        <title>RNAPlot - {entry_id}</title>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
        <link rel="stylesheet" href="../assets/fornac.css" />
        <script src="https://unpkg.com/jquery"></script>
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
        <script src="https://unpkg.com/d3@3.5"></script>
        <script src="https://unpkg.com/d3-grid"></script>
        <script src="../assets/fornac.js"></script>
        <style>
            svg {{
                width: 100%;
                height: 100%;
                max-width: 900px;
                border: 1px solid gray;
            }}

            td {{
                padding: 0 15px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h2>Target ID: {entry_id}</h2>

            <div>
                <a href="https://www.ncbi.nlm.nih.gov/clinvar/variation/{entry_id}/" target="_blank">
                    https://www.ncbi.nlm.nih.gov/clinvar/variation/{entry_id}/
                </a>
            </div>

            <div id="reference">
                <h4>Reference</h4>
                <table>
                    <tr>
                        <td>Sequence</td>
                        <td><textarea disabled rows="4" cols="80" style="font-family: monospace">{original_reference_sequence}</textarea></td>
                    </tr>
                    <tr>
                        <td>Base</td>
                        <td><textarea disabled rows="1" cols="80" style="font-family: monospace">{reference_base}</textarea></td>
                    </tr>
                    <tr>
                        <td>Structure</td>
                        <td><textarea disabled rows="4" cols="80" style="font-family: monospace">{reference_structure}</textarea></td>
                    </tr>
                    <tr>
                        <td>Colors</td>
                        <td><textarea disabled rows="2" cols="80" style="font-family: monospace">{reference_colors}</textarea></td>
                    </tr>
                </table>

                <button type="button" class="btn btn-link" onclick="showStructure('reference_structure')">Show reference structure</button>
            </div>

            <div id="reference_structure" style="display: none"></div>

            <div id="alternate">
                <h4>Alternate</h4>
                <table>
                    <tr>
                        <td>Sequence</td>
                        <td><textarea disabled rows="4" cols="80" style="font-family: monospace">{original_alternate_sequence}</textarea></td>
                    </tr>
                    <tr>
                        <td>Base</td>
                        <td><textarea disabled rows="1" cols="80" style="font-family: monospace">{alternate_base}</textarea></td>
                    </tr>
                    <tr>
                        <td>Structure</td>
                        <td><textarea disabled rows="4" cols="80" style="font-family: monospace">{alternate_structure}</textarea></td>
                    </tr>
                    <tr>
                        <td>Colors</td>
                        <td><textarea disabled rows="2" cols="80" style="font-family: monospace">{alternate_colors}</textarea></td>
                    </tr>
                </table>
                
                <button type="button" class="btn btn-link" onclick="showStructure('alternate_structure')">Show alternate structure</button>
            </div>

            <div id="alternate_structure" style="display: none"></div>

            <script type="text/javascript">
                var ref_is_built = false;
                var alt_is_built = false;

                function buildReference() {{
                    var reference_structure = "{reference_structure}";
                    var reference_sequence = "{reference_sequence}";
                    var reference_colors = "{reference_colors}";

                    var reference_data = {{
                        "structure": reference_structure,
                        "sequence": reference_sequence,
                        "name": "Reference"
                    }};

                    let reference_container = new fornac.FornaContainer(
                        "#reference_structure",
                        {{"animation": false, "zoomable": true, "initialSize": [900,700]}}
                    );

                    reference_container.addRNA(reference_data.structure, reference_data);
                    reference_container.addCustomColorsText(reference_colors);
                }}

                function buildAlternate() {{
                    var alternate_structure = "{alternate_structure}";
                    var alternate_sequence = "{alternate_sequence}";
                    var alternate_colors = "{alternate_colors}";

                    var alternate_data = {{
                        "structure": alternate_structure,
                        "sequence": alternate_sequence,
                        "name": "Alternate"
                    }};

                    let alternate_container = new fornac.FornaContainer(
                        "#alternate_structure",
                        {{"animation": false, "zoomable": true, "initialSize": [900,700]}}
                    );

                    alternate_container.addRNA(alternate_data.structure, alternate_data);
                    alternate_container.addCustomColorsText(alternate_colors);
                }}

                function showStructure(sectionId) {{
                    var section = document.getElementById(sectionId);
                    if (section.style.display === "none") {{
                        if (sectionId === "reference_structure" && !ref_is_built) {{
                            buildReference();
                            ref_is_built = true;
                        }} else if (sectionId === "alternate_structure" && !alt_is_built) {{
                            buildAlternate();
                            alt_is_built = true;
                        }}
                        section.style.display = "block";
                    }} else {{
                        section.style.display = "none";
                    }}
                }}
            </script>
        </div>
    </body>
</html>
    """


def process_vcf_entry(
    vcf_entry_id,
    reference_sequence,
    alternate_sequence,
    bed_entry,
    vcf_entry,
    search_for="[AGT][AG]AC[ACT]",
    out_table=None,
    out_data_folder=None
):
    """
    Process a VCF entry

    :param vcf_entry_id         ID of the VCF entry
    :param reference_sequence:  Reference sequence
    :param alternate_sequence:  Alternate sequence (reference with modification)
    :param bed_entry:           Genomic coordinates of a specific region
    :param vcf_entry:           VCF entry
    :param search_for:          Regex (must be of a fixed length)
    :param out_table:           Path to the output table
    :param out_data_folder:     Path to the output folder with the HTML index
    """

    out_table_row = list()
    out_index_row = None

    # Check whether the reference sequence actually contains at least one pattern site
    regex = re.compile(search_for)

    if not regex.search(reference_sequence):
        return vcf_entry_id, list(), None

    # Search for patterns
    # Focus on SNPs occurring outside the patterns
    snp_in_pattern = False

    # Define the length of the fixed-length regex
    pattern_len = 0

    open_bracket = False

    for i in range(len(search_for)):
        if not open_bracket:
            pattern_len += 1

        if search_for[i] == "[":
            open_bracket = True
        
        if search_for[i] == "]":
            open_bracket = False

    for pos in range(len(reference_sequence)):
        if pos + pattern_len == len(reference_sequence):
            break

        ref_window = reference_sequence[pos: pos + pattern_len]
        alt_window = alternate_sequence[pos: pos + pattern_len]

        if re.match(search_for, ref_window):
            if "_" in alt_window:
                snp_in_pattern = True

            elif ref_window != alt_window:
                snp_in_pattern = True

        elif "_" in ref_window:
            # Expand window
            count_missing_bases = reference_sequence.count("_")

            ref_window = reference_sequence[pos: pos + pattern_len + count_missing_bases].replace("_", "")[:pattern_len]

            if re.match(search_for, ref_window):
                snp_in_pattern = True

        if snp_in_pattern:
            break

    if not snp_in_pattern:
        reference_no_missing_bases = reference_sequence.replace("_", "")
        reference_structure = None
        reference_min_free_energy = None

        alternate_no_missing_bases = alternate_sequence.replace("_", "")
        alternate_structure = None
        alternate_min_free_energy = None

        if reference_no_missing_bases.count("TGA") != alternate_no_missing_bases.count("TGA") \
            or reference_no_missing_bases.count("TAA") != alternate_no_missing_bases.count("TAA") \
            or reference_no_missing_bases.count("TAG") != alternate_no_missing_bases.count("TAG"):
            # Do not report the entry if the SNP forms a termination codon (TGA/TAA/TAG) or changes termination codon to a sense codon
            return vcf_entry_id, list(), None

        # Does the modification cause a frameshift in case of indel?
        frameshift = len(reference_no_missing_bases) != len(alternate_no_missing_bases)

        if frameshift:
            # Do not report the entry in case of frameshifts
            return vcf_entry_id, list(), None

        with tempfile.NamedTemporaryFile() as rnafold_in, tempfile.NamedTemporaryFile() as rnafold_out:
            with open(rnafold_in.name, "wt") as rnafold_in_file:
                # Write reference sequence
                rnafold_in_file.write(">reference_sequence\n{}\n".format(reference_no_missing_bases))

                # Write alternate sequence
                rnafold_in_file.write(">alternate_sequence\n{}\n".format(alternate_no_missing_bases))

            try:
                # Run ViennaRNA:RNAfold on reference and alternate sequence
                subprocess.check_call(
                    [
                        "RNAfold",
                        "-i",
                        rnafold_in.name,
                        "--noLP",
                    ],
                    stdout=rnafold_out,
                    stderr=subprocess.DEVNULL
                )

            except subprocess.CalledProcessError as e:
                raise Exception("An error has occurred while running RNAfold").with_traceback(e.__traceback__)

            # Get mfe structures and minimum free energies in bracket notation
            with open(rnafold_out.name, "rt") as rnafold_out_file:
                last_seq_id = None

                for line in rnafold_out_file:
                    line = line.strip()

                    if line:
                        if line.startswith(">"):
                            last_seq_id = line[1:]

                        if line.startswith(".") or line.startswith("("):
                            line_split = line.split(" ")

                            if last_seq_id == "reference_sequence":
                                reference_structure = line_split[0]
                                reference_min_free_energy = float(line_split[-1][1:-1])  # kcal/mol

                            elif last_seq_id == "alternate_sequence":
                                alternate_structure = line_split[0]
                                alternate_min_free_energy = float(line_split[-1][1:-1])  # kcal/mol

        if reference_structure != alternate_structure:
            # Folded structures of reference and alternate sequences are different here for sure
            # Check whether the alternate structure is paired/unpaired on the target patterns
            # Focus on SNPs occurring outside patterns
            at_least_one_pattern_unpaired = False
            at_least_one_pattern_paired = False

            # Take track of the positions of the patterns in the reference and alternate sequences
            reference_pattern_positions = list()
            alternate_pattern_positions = list()

            for match in regex.finditer(reference_no_missing_bases):
                reference_pattern_positions.append(match.start() + 1)

                ref_match_pos = match.start()

                if len(reference_no_missing_bases) == len(alternate_no_missing_bases):
                    alt_match_pos = match.start()

                elif len(reference_no_missing_bases) < len(alternate_no_missing_bases):
                    delta = len(alternate_no_missing_bases) - len(reference_no_missing_bases)
                    missing_base_pos = reference_sequence.index("_")

                    if match.start() > missing_base_pos:
                        alt_match_pos = match.start() + reference_sequence.count("_")

                    else:
                        alt_match_pos = match.start()

                elif len(reference_no_missing_bases) > len(alternate_no_missing_bases):
                    delta = len(reference_no_missing_bases) - len(alternate_no_missing_bases)
                    missing_base_pos = alternate_sequence.index("_")

                    if match.start() > missing_base_pos:
                        alt_match_pos = match.start() - alternate_sequence.count("_")

                    else:
                        alt_match_pos = match.start()

                # Check whether the structure at the pattern site on changed from paired to unpaired or the other way around
                # A dot in the dot-bracket notation represents an unpaired position (free)
                how_many_dots_reference = reference_structure[ref_match_pos:ref_match_pos + pattern_len].count(".")
                how_many_dots_alternate = alternate_structure[alt_match_pos:alt_match_pos + pattern_len].count(".")

                if how_many_dots_reference == 0 and how_many_dots_alternate > 0:
                    at_least_one_pattern_unpaired = True

                if how_many_dots_reference > 0 and how_many_dots_alternate == 0:
                    at_least_one_pattern_paired = True

            for match in regex.finditer(alternate_no_missing_bases):
                # Also count the number of pattern sites in the alternate sequence
                alternate_pattern_positions.append(match.start() + 1)

            if at_least_one_pattern_unpaired or at_least_one_pattern_paired:
                structure_state = "paired-unpaired" if at_least_one_pattern_unpaired else "unpaired-paired"
                reference_pattern_sites = len(reference_pattern_positions)
                alternate_pattern_sites = len(alternate_pattern_positions)

                # Take track of different bases positions
                reference_snp_positions = list()
                alternate_snp_positions = list()

                # Check if synonymous or non-synonymous
                synonymous = "no"

                if len(reference_no_missing_bases) == len(alternate_no_missing_bases):
                    for pos in range(len(reference_no_missing_bases)):
                        if reference_no_missing_bases[pos].lower() != alternate_no_missing_bases[pos].lower():
                            reference_snp_positions.append(pos + 1)
                            alternate_snp_positions.append(pos + 1)

                    # Check in case of single-nucleotide variation only
                    synonymous = "yes" if str(Seq.Seq(reference_no_missing_bases)) == str(Seq.Seq(alternate_no_missing_bases)) else "no"

                elif len(reference_no_missing_bases) > len(alternate_no_missing_bases):
                    missing_bases = alternate_sequence.count("_")
                    first_missing_base_pos = alternate_sequence.index("_")
                    reference_snp_positions.extend([pos for pos in range(first_missing_base_pos, first_missing_base_pos + missing_bases + 1)])

                    alternate_snp_positions.append(first_missing_base_pos)

                elif len(reference_no_missing_bases) < len(alternate_no_missing_bases):
                    missing_bases = reference_sequence.count("_")
                    first_missing_base_pos = reference_sequence.index("_")
                    alternate_snp_positions.extend([pos for pos in range(first_missing_base_pos, first_missing_base_pos + missing_bases + 1)])

                    reference_snp_positions.append(first_missing_base_pos)

                # Define the row for the output table
                out_table_row = [
                    vcf_entry_id,
                    bed_entry["chromosome"],
                    bed_entry["start"],
                    bed_entry["end"],
                    bed_entry["strand"],
                    bed_entry["gene_symbol"],
                    bed_entry["transcript_id"],
                    reference_no_missing_bases,
                    reference_structure,
                    round(reference_min_free_energy, 4),
                    reference_pattern_sites,
                    ",".join([str(pos) for pos in reference_pattern_positions]),
                    vcf_entry["pos"],
                    alternate_no_missing_bases,
                    alternate_structure,
                    round(alternate_min_free_energy, 4),
                    alternate_pattern_sites,
                    ",".join([str(pos) for pos in alternate_pattern_positions]),
                    structure_state,
                    round(abs(reference_min_free_energy) - abs(alternate_min_free_energy), 4),
                    synonymous,
                    "yes" if frameshift else "no"
                ]

                if out_data_folder:
                    # Build the current VCF entry html page
                    # Also define the row for the output index table
                    forna_filepath = os.path.join(out_data_folder, "{}.html".format(vcf_entry_id))

                    with open(forna_filepath, "w+") as forna_file:
                        # Define pattern sites colors in the reference sequence
                        reference_pattern_colors = " ".join(
                            ["{}-{}:green".format(pattern_pos, pattern_pos + pattern_len - 1) for pattern_pos in reference_pattern_positions]
                        )

                        # Define SNPs color in the reference sequence
                        reference_snp_colors = " ".join(["{}:orange".format(snp_pos) for snp_pos in reference_snp_positions])

                        # Define pattern sites colors in the alternate sequence
                        alternate_pattern_colors = " ".join(
                            ["{}-{}:green".format(pattern_pos, pattern_pos + pattern_len - 1) for pattern_pos in alternate_pattern_positions]
                        )

                        # Define SNPs color in the alternate sequence
                        alternate_snp_colors = " ".join(["{}:orange".format(snp_pos) for snp_pos in alternate_snp_positions])

                        # Also define the html page
                        forna_file.write(
                            forna_template().format(
                                entry_id=vcf_entry_id,
                                original_reference_sequence=reference_sequence,
                                reference_sequence=reference_no_missing_bases,
                                reference_base="".join([reference_no_missing_bases[pos - 1] for pos in reference_snp_positions]),
                                reference_structure=reference_structure,
                                reference_colors="{} {}".format(reference_pattern_colors, reference_snp_colors),
                                original_alternate_sequence=alternate_sequence,
                                alternate_sequence=alternate_no_missing_bases,
                                alternate_base="".join([alternate_no_missing_bases[pos - 1] for pos in alternate_snp_positions]),
                                alternate_structure=alternate_structure,
                                alternate_colors="{} {}".format(alternate_pattern_colors, alternate_snp_colors)
                            )
                        )

                    # Define the index row
                    out_index_row = forna_index(header=TABLE_COLS, row=out_table_row, filepath=None)

    return vcf_entry_id, out_table_row, out_index_row


def main():
    # Load command line parameters
    args = read_params(sys.argv[1:])

    for input_file in [args.bed, args.vcf, args.genome]:
        if not os.path.isfile(input_file):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), input_file)

    if not regex_is_valid(args.search):
        raise ValueError("The --search argument must be a fixed-length regex")

    if args.out_table:
        with open(args.out_table, "w+") as table_file:
            table_file.write("{}\n".format("\t".join(TABLE_COLS)))

    out_data_folder = None

    if args.out_html:
        out_data_index = os.path.join(args.out_html, "index.html")

        # Path to the folder with .html files
        out_data_folder = os.path.join(args.out_html, "data")
        os.makedirs(out_data_folder, exist_ok=True)

        # Path to the folder with "forna" assets
        forna_folder = os.path.join(args.out_html, "assets")
        os.makedirs(forna_folder, exist_ok=True)

        # Init target index page
        forna_index(init=True, header=TABLE_COLS, filepath=out_data_index)

        # Retrieve "forna" CSS and JS
        for asset in FORNA_ASSETS:
            try:
                asset_filepath = os.path.join(args.out_html, os.path.basename(asset))

                if not os.path.isfile(asset_filepath):
                    wget.download(asset, out=os.path.abspath(forna_folder))

            except:
                # Do not block the execution here since the host could not have access to internet
                print("Unable to download {}".format(asset))

    # Consider the SNPs occurring on the genomic regions defined in the input BED file
    # Load the BED file into a pandas DataFrame
    # Assume the first three columns are Chromosome, Start, End, and Strand. Does not need any other info
    bed = pd.read_csv(
        args.bed,
        header=None,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["chromosome", "start", "end", "name", "score", "strand"],
        sep="\t",
        skiprows=args.bed_skip_lines
    )

    # Open vcf file
    vcf = pysam.VariantFile(args.vcf)

    # Open fasta file
    genome = pysam.FastaFile(args.genome)

    process_vcf_entry_partial = partial(
        process_vcf_entry,
        search_for=args.search,
        out_table=args.out_table,
        out_data_folder=out_data_folder
    )

    out_table_rows = list()

    out_index_rows = list()

    processed = set()

    # Run prediction on folds in parallel
    with mp.Pool(processes=args.nproc) as pool, tqdm.tqdm() as pbar:
        # Wrapper around the update function of tqdm
        def progress(*args):
            pbar.update()

        jobs = list()

        for entry in vcf:
            if not entry.alts:
                # In case of deletions, alts is None
                # Define alts as a tuple with a "." character
                entry.alts = (".")

            for alt in entry.alts:
                # Check whether there exists at least one BED region for the current VCF entry before spawning a new job
                # Retrieve the genomic region in which the VCF entry occurs in
                # VCF lines are always assumed to be '+' strand, as VCF doesn't specify that attribute
                bed_regions = bed.loc[
                    (bed["chromosome"] == "chr{}".format(entry.chrom)) & 
                    (bed["start"] <= entry.pos) & (entry.pos <= bed["end"])
                ]

                alt = alt.upper()
                ref = entry.ref.upper()

                if not bed_regions.empty:
                    bed_count = 1

                    # Iterate over the target genomic regions
                    # Process the VCF entry in case of an intersection with at least one genomic region in the BED file
                    for _, row in bed_regions.iterrows():
                        # Use the start and end positions of the genomic region by default
                        region_start = row["start"]
                        region_end = row["end"]

                        # In case of indels only
                        indel_len = max(len(ref), len(alt))

                        if entry.pos + indel_len > region_end:
                            new_region_end = region_end + indel_len - (region_end - entry.pos)

                            if region_end != new_region_end:
                                region_end = new_region_end

                        # Retrieve the sequence from the input genome
                        # Within pysam, coordinates are 0-based half-open intervals
                        # i.e., the first base of the reference sequence is numbered zero;
                        # and the base at position start is part of the interval, but the base at end is not.
                        reference_sequence = genome.fetch("chr{}".format(entry.chrom), region_start - 1, region_end).upper()

                        # Check whether the REF base in the VCF entry actually matches with the base in the reference sequence under that specific position
                        if ref != "." and reference_sequence[entry.pos - region_start:entry.pos - region_start + len(ref)] != ref:
                            print(entry)
                            print(row)
                            print(region_start)
                            print(region_end)
                            print(reference_sequence)
                            raise Exception("References mismatch! Are you using the right reference genome sequence?")

                        # Initialize the alternate sequence
                        alternate_sequence = reference_sequence

                        if ref == ".":
                            # Insertion
                            reference_sequence = "{}{}{}".format(
                                reference_sequence[:entry.pos - region_start],
                                "_" * len(alt),
                                reference_sequence[entry.pos - region_start:]
                            )

                            alternate_sequence = "{}{}{}".format(
                                alternate_sequence[:entry.pos - region_start],
                                alt,
                                alternate_sequence[entry.pos - region_start:]
                            )

                        elif alt == ".":
                            # Deletion
                            alternate_sequence = "{}{}{}".format(
                                alternate_sequence[:entry.pos - region_start],
                                "_" * len(ref),
                                alternate_sequence[entry.pos - region_start + len(ref):]
                            )

                        else:
                            # Replacement
                            if len(ref) > len(alt):
                                delta = len(ref) - len(alt)

                                alternate_sequence = "{}{}{}{}".format(
                                    alternate_sequence[:entry.pos - region_start],
                                    alt,
                                    "_" * delta,
                                    alternate_sequence[entry.pos - region_start + len(ref):]
                                )

                            elif len(ref) < len(alt):
                                delta = len(alt) - len(ref)

                                reference_sequence = "{}{}{}".format(
                                    reference_sequence[:entry.pos - region_start + len(ref)],
                                    "_" * delta,
                                    reference_sequence[entry.pos - region_start + len(ref):]
                                )

                                alternate_sequence = "{}{}{}".format(
                                    alternate_sequence[:entry.pos - region_start],
                                    alt,
                                    alternate_sequence[entry.pos - region_start + len(ref):]
                                )

                            else:
                                alternate_sequence = "{}{}{}".format(
                                    alternate_sequence[:entry.pos - region_start],
                                    alt,
                                    alternate_sequence[entry.pos - region_start + len(alt):]
                                )

                        transcript_id = row["name"].split("__")[0]

                        # Retrieve all the genomic regions in the BED file corresponding to the current transcript ID
                        transcript_regions = bed.loc[bed["name"].str.startswith("{}__".format(transcript_id))]
                        
                        # Sort regions based on the start position in ascending order
                        transcript_regions.sort_values(by="start", ascending=True, inplace=True)

                        merged_reference_sequence = ""
                        merged_alternate_sequence = ""

                        absolute_start = -1
                        absolute_end = -1

                        for _, transcript_row in transcript_regions.iterrows():
                            if absolute_start < 0:
                                absolute_start = transcript_row["start"]
                            
                            absolute_end = transcript_row["end"]

                            if transcript_row["name"] == row["name"]:
                                merged_reference_sequence += reference_sequence
                                merged_alternate_sequence += alternate_sequence
                            
                            else:
                                partial_seq = genome.fetch("chr{}".format(entry.chrom), transcript_row["start"] - 1, transcript_row["end"]).upper()

                                merged_reference_sequence += partial_seq
                                merged_alternate_sequence += partial_seq

                        if row["strand"] == "-":
                            # Complement
                            comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "_": "_"}
                            rev_comp_reference_sequence = "".join([comp[n] for n in merged_reference_sequence.upper()])
                            rev_comp_alternate_sequence = "".join([comp[n] for n in merged_alternate_sequence.upper()])

                            # Reverse
                            rev_comp_reference_sequence = rev_comp_reference_sequence[::-1]
                            rev_comp_alternate_sequence = rev_comp_alternate_sequence[::-1]

                            merged_reference_sequence = rev_comp_reference_sequence
                            merged_alternate_sequence = rev_comp_alternate_sequence

                        bed_entry = {
                            "chromosome": row["chromosome"],
                            "start": absolute_start,
                            "end": absolute_end,
                            "strand": row["strand"],
                            "gene_symbol": row["name"].split("__")[1],  # To use with the selection_with_protein_coding_genes.bed file
                            "transcript_id": row["name"].split("__")[0],  # Take track of the transcript ID
                        }

                        jobs.append(
                            pool.apply_async(
                                process_vcf_entry_partial,
                                args=(
                                    entry.id,
                                    merged_reference_sequence,
                                    merged_alternate_sequence,
                                    bed_entry,
                                    {"chrom": entry.chrom, "pos": entry.pos}
                                ),
                                callback=progress,
                            )
                        )

                        bed_count += 1

        # Get results from jobs
        for job in jobs:
            entry_id, partial_out_table_row, partial_out_index_row = job.get()

            if entry_id not in processed:
                if partial_out_table_row:
                    out_table_rows.append(partial_out_table_row)

                if partial_out_index_row:
                    out_index_rows.append(partial_out_index_row)

                processed.add(entry_id)

    genome.close()

    vcf.close()

    if args.out_table and out_table_rows:
        for row in out_table_rows:
            with open(args.out_table, "a+") as table_file:
                table_file.write("{}\n".format("\t".join([str(val) for val in row])))

    if args.out_html:
        if out_index_rows:
            for row in out_index_rows:
                with open(out_data_index, "a+") as index_file:
                    index_file.write("{}\n".format(row))

        # Finalize forna index page
        forna_index(close=True, header=TABLE_COLS, filepath=out_data_index)


if __name__ == "__main__":
    main()