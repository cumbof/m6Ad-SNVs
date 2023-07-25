#!/usr/bin/env python3
"""
Extract the last 2 CDS + 100bp of adjacent UTR
Consider protein_coding genes only
"""

import copy
import os

def extract_exons_utr_coordinates(bed_entry):
   # Split the BED entry into its components
   chrom, chrom_start, chrom_end, gene_name, score, strand, thickStart, thickEnd, _, block_count, block_sizes, block_starts = bed_entry.strip().split("\t")

   # Convert values to integers and lists
   chrom_start, chrom_end, thickStart, thickEnd = (
      int(chrom_start),
      int(chrom_end),
      int(thickStart),
      int(thickEnd),
   )
   
   block_count = int(block_count)
   block_sizes = list(map(int, block_sizes.strip(",").split(",")))
   block_starts = list(map(int, block_starts.strip(",").split(",")))

   # Initialize lists to store exon and UTR coordinates
   coding_exon_coordinates = list()
   utr_coordinates5 = list()
   utr_coordinates3 = list()

   if thickStart != thickEnd:
      exon_starts = [x + chrom_start for x in block_starts]
      exon_ends = [x + y for x, y in zip(exon_starts, block_sizes)]

      for start, end in zip(exon_starts, exon_ends):
         if start < thickStart:
            if strand == "+":
               utr_coordinates5.append((start, min(end, thickStart)))
            
            elif strand == "-":
               utr_coordinates3.append((start, min(end, thickStart)))
                
            # split block
            if end > thickStart:
               if end > thickEnd:
                  coding_exon_coordinates.append((thickStart, thickEnd))
                  
                  if strand == "+":
                     utr_coordinates3.append((thickEnd, end))
                  
                  elif strand == "-":
                     utr_coordinates5.append((thickEnd, end))
                    
               else:
                  coding_exon_coordinates.append((thickStart, end))
            
         elif start < thickEnd:
            coding_exon_coordinates.append((start, min(end, thickEnd)))
            
            if end > thickEnd:
               if strand == "+":
                  utr_coordinates3.append((thickEnd, end))
               
               elif strand == "-":
                  utr_coordinates5.append((thickEnd, end))
            
         else:
            if strand == "+":
               utr_coordinates3.append((start, end))
            
            elif strand == "-":
               utr_coordinates5.append((start, end))

   # These coordinates are for a "+" stranded region, will need to manipulated to represent a "-" stranded region
   return coding_exon_coordinates, utr_coordinates5, utr_coordinates3

"""
bed_data = [
   "chr1",
   "67092164",
   "67134970",
   "NM_001276352.2",
   "0",
   "-",
   "67093579",
   "67127240",
   "0",
   "9",
   "1440,70,145,68,113,158,92,86,41,",
   "0,4087,11073,19412,23187,33587,35001,38977,42765,",
]

coding_exon_coordinates, utr_coordinates5, utr_coordinates3 = extract_exons_utr_coordinates("\t".join(bed_data))

print(utr_coordinates5)
print(coding_exon_coordinates)
print(utr_coordinates3)
"""

last_cds = 2
utr_size = 100

regions = dict()

#whole_gene_filepath = "whole_gene__ncbiRefSeq.bed"
whole_gene_filepath = "whole_gene__gencode.bed"

transcripts_db = os.path.splitext(whole_gene_filepath)[0].split("__")[-1]

selection_filepath = "selection__last_{}_CDS_plus_{}_UTR__{}.bed".format(last_cds, utr_size, transcripts_db)

with open(selection_filepath, "w+") as selection_bed:
   with open(whole_gene_filepath) as whole_gene_bed:
      for line in whole_gene_bed:
         line = line.strip()
         if line:
            if not line.startswith("chr"):
               selection_bed.write("{}\n".format(line))
            else:
               line_split = line.split("\t")
               # Select RefSeq transcript identifiers that start with NM_ (coding)
               if (line_split[3].startswith("NM_") and transcripts_db == "ncbiRefSeq") or transcripts_db != "ncbiRefSeq":
                  coding_exon_coordinates, utr_coordinates5, utr_coordinates3 = extract_exons_utr_coordinates(line)
                  if coding_exon_coordinates:
                     # Define the new line
                     newline = line_split[:6]
                     entry_id = newline[3]

                     strand = line_split[5]

                     utrs = list()

                     if strand == "+":
                        coding_exon_coordinates = coding_exon_coordinates[-last_cds:]
                        
                        utr_len = utr_size
                        for utr_start, utr_end in utr_coordinates3:
                           utr_len -= (utr_end - utr_start)

                           if utr_len > 0:
                              utrs.append((utr_start, utr_end))

                           else:
                              utrs.append((utr_start, utr_end - abs(utr_len)))

                              break

                     elif strand == "-":
                        coding_exon_coordinates = coding_exon_coordinates[:last_cds]

                        utr_len = utr_size
                        for utr_start, utr_end in reversed(utr_coordinates3):
                           utr_len -= (utr_end - utr_start)

                           if utr_len > 0:
                              utrs.append((utr_start, utr_end))
                           
                           else:
                              utrs.append((utr_start + abs(utr_len), utr_end))

                              break

                     if entry_id not in regions:
                        regions[entry_id] = list()

                     for exon_count, (exon_start, exon_end) in enumerate(coding_exon_coordinates):
                        newline_list = copy.deepcopy(newline)
                        newline_list[1] = str(exon_start)
                        newline_list[2] = str(exon_end)
                        newline_list[3] = "{}__exon_CDS_{}".format(entry_id, exon_count + 1)
                        selection_bed.write("{}\n".format("\t".join(newline_list)))

                        regions[entry_id].append(newline_list)
                     
                     for utr_count, (utr_start, utr_end) in enumerate(utrs):
                        newline_list = copy.deepcopy(newline)
                        newline_list[1] = str(utr_start)
                        newline_list[2] = str(utr_end)
                        newline_list[3] = "{}__UTR_{}".format(entry_id, utr_count + 1)
                        selection_bed.write("{}\n".format("\t".join(newline_list)))

                        regions[entry_id].append(newline_list)

# Load the GENCODE.v43 annotation
# Focus on protein_coding only
genes = dict()

with open("gencode.v43.genes.bed") as m:
   for line in m:
      line = line.strip()
      if line:
         line_split = line.split("\t")
         if line_split[0] not in genes:
            genes[line_split[0]] = {}
         info = [x.strip() for x in line_split[9].strip().split(" ")]
         gene_type = info[info.index("gene_type")+1].replace("\"", "").replace(";", "").strip()
         if gene_type == "protein_coding":
            gene = info[info.index("gene_name")+1].replace("\"", "").replace(";", "").strip()
            genes[line_split[0]][gene] = {
               "start": int(line_split[1]),
               "end": int(line_split[2]),
               "strand": line_split[5]
            }

unique_regions = set()

limit_size = 250

selection_genes_filepath = "selection__last_{}_CDS_plus_{}_UTR__{}__with_protein_coding_genes.bed".format(last_cds, utr_size, transcripts_db)
selection_genes_limit_filepath = "selection__last_{}_CDS_plus_{}_UTR__{}__with_protein_coding_genes__{}nt.bed".format(last_cds, utr_size, transcripts_db, limit_size)

# Apply the annotation to the selection
# Add the Ensembl/RefSeq transcript ID with the gene symbol
with open(selection_genes_filepath, "w+") as selection_with_protein_coding_genes, open(selection_genes_limit_filepath, "w+") as selection_limit_with_protein_coding_genes:
   for entry_id in regions:
      chrom = regions[entry_id][0][0]
      
      if chrom in genes:
         start = min([int(bed_data[1]) for bed_data in regions[entry_id]])
         end = max([int(bed_data[2]) for bed_data in regions[entry_id]])
         strand = regions[entry_id][0][5]

         for gene in genes[chrom]:
            if genes[chrom][gene]["start"] <= start and genes[chrom][gene]["end"] >= end and genes[chrom][gene]["strand"] == strand:
               region_list = list()
               for bed_data in regions[entry_id]:
                  region_list.append(bed_data[1])
                  region_list.append(bed_data[2])

               region_list = tuple(sorted(region_list))

               if region_list not in unique_regions:
                  bed_data_sorted = sorted(regions[entry_id], key=lambda data: data[1], reverse=True)
                  if bed_data_sorted[0][5] == "-":
                     bed_data_sorted.reverse()

                  bed_data_len = limit_size

                  stop_processing_regions = False
                  last_of_limited = False

                  for bed_data in bed_data_sorted:
                     bed_data_list = copy.deepcopy(bed_data)
                     bed_data_list_limit = copy.deepcopy(bed_data)
                     
                     bed_data_len -= (int(bed_data[2]) - int(bed_data[1]))

                     if bed_data_len <= 0:
                        if strand == "-":
                           bed_data_list_limit[2] = str(int(bed_data_list_limit[2]) - abs(bed_data_len))
                        elif strand == "+":
                           bed_data_list_limit[1] = str(int(bed_data_list_limit[1]) + abs(bed_data_len))

                        last_of_limited = True

                     name = bed_data[3].split("__")

                     bed_data_list[3] = "{}__{}__{}".format(name[0], gene, name[1])
                     bed_data_list_limit[3] = "{}__{}__{}".format(name[0], gene, name[1])

                     selection_with_protein_coding_genes.write("{}\n".format("\t".join(bed_data_list)))

                     if not stop_processing_regions:
                        selection_limit_with_protein_coding_genes.write("{}\n".format("\t".join(bed_data_list_limit)))
                        
                        if last_of_limited:
                           stop_processing_regions = True

                  unique_regions.add(region_list)

                  break
