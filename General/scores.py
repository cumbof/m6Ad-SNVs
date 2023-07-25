import os
import re

tablepath = "./targets-ClinVar-m6A-gencode-20230724/targets.tsv"

headerline = [
    "Free DRA sites (reference)",
    "Free DRA sites (alternate)",
    "Free RAC sites (reference)",
    "Free RAC sites (alternate)",
    "Free ACH sites (reference)",
    "Free ACH sites (alternate)",
    "Free DRAC sites (reference)",
    "Free DRAC sites (alternate)",
    "Free RACH sites (reference)",
    "Free RACH sites (alternate)",
    "Free DRACH sites (reference)",
    "Free DRACH sites (alternate)"
]

print(os.path.abspath(tablepath))

print("<th>{}</th>".format("</th><th>".join(headerline)))

with open(tablepath) as targets:
    header = targets.readline().strip().split("\t")
    for line in targets:
        line = line.strip()
        if line:
            line_split = line.split("\t")
            
            reference_sequence = line_split[header.index("reference_sequence")]
            reference_structure = line_split[header.index("reference_structure")]
            alternate_sequence = line_split[header.index("alternate_sequence")]
            alternate_structure = line_split[header.index("alternate_structure")]

            regex = re.compile("[AGT][AG]AC[ACT]")

            drach_ref_free = 0

            rach_ref_free = 0
            drac_ref_free = 0

            ach_ref_free = 0
            rac_ref_free = 0
            dra_ref_free = 0

            for match_ref in regex.finditer(reference_sequence):
                # DRACH
                drach_ref_structure = reference_structure[match_ref.start():match_ref.start() + 5]
                if drach_ref_structure == ".....":
                    drach_ref_free += 1
                
                else:
                    check_dra = True
                    check_rac = True
                    check_ach = True

                    # DRAC
                    drac_ref_structure = reference_structure[match_ref.start():match_ref.start() + 4]
                    # RACH
                    rach_ref_structure = reference_structure[match_ref.start() + 1:match_ref.start() + 5]

                    if drac_ref_structure == "....":
                        drac_ref_free += 1
                        check_dra = False
                        check_rac = False
                    
                    elif rach_ref_structure == "....":
                        rach_ref_free += 1
                        check_rac = False
                        check_ach = False

                    if check_dra:
                        # DRA
                        dra_ref_structure = reference_structure[match_ref.start():match_ref.start() + 3]
                        if dra_ref_structure == "...":
                            dra_ref_free += 1
                    
                    if check_rac:
                        # RAC
                        rac_ref_structure = reference_structure[match_ref.start() + 1:match_ref.start() + 4]
                        if rac_ref_structure == "...":
                            rac_ref_free += 1

                    if check_ach:
                        # ACH
                        ach_ref_structure = reference_structure[match_ref.start() + 2:match_ref.start() + 5]
                        if ach_ref_structure == "...":
                            ach_ref_free += 1
            
            drach_alt_free = 0

            rach_alt_free = 0
            drac_alt_free = 0

            ach_alt_free = 0
            rac_alt_free = 0
            dra_alt_free = 0

            for match_alt in regex.finditer(alternate_sequence):
                # DRACH
                drach_alt_structure = alternate_structure[match_alt.start():match_alt.start() + 5]
                if drach_alt_structure == ".....":
                    drach_alt_free += 1
                
                else:
                    check_dra = True
                    check_rac = True
                    check_ach = True

                    # DRAC
                    drac_alt_structure = alternate_structure[match_alt.start():match_alt.start() + 4]
                    # RACH
                    rach_alt_structure = alternate_structure[match_alt.start() + 1:match_alt.start() + 5]

                    if drac_alt_structure == "....":
                        drac_alt_free += 1
                        check_dra = False
                        check_rac = False
                    
                    elif rach_alt_structure == "....":
                        rach_alt_free += 1
                        check_rac = False
                        check_ach = False

                    if check_dra:
                        # DRA
                        dra_alt_structure = alternate_structure[match_alt.start():match_alt.start() + 3]
                        if dra_alt_structure == "...":
                            dra_alt_free += 1
                    
                    if check_rac:
                        # RAC
                        rac_alt_structure = alternate_structure[match_alt.start() + 1:match_alt.start() + 4]
                        if rac_alt_structure == "...":
                            rac_alt_free += 1

                    if check_ach:
                        # ACH
                        ach_alt_structure = alternate_structure[match_alt.start() + 2:match_alt.start() + 5]
                        if ach_alt_structure == "...":
                            ach_alt_free += 1

            print(
                "<td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>".format(
                    dra_ref_free,
                    dra_alt_free,
                    rac_ref_free,
                    rac_alt_free,
                    ach_ref_free,
                    ach_alt_free,
                    drac_ref_free,
                    drac_alt_free,
                    rach_ref_free,
                    rach_alt_free,
                    drach_ref_free,
                    drach_alt_free
                )
            )
