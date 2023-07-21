#!/usr/bin/env python3
"""
Extract the last 2 CDS + 100bp of adjacent UTR
Consider protein_coding genes only
"""

max_cds = 2
expand_utr = 100

region_max_size = 250
#region_max_size = 0  # disabled

regions = set()

with open("selection_250nt.bed", "w+") as m1:
   with open("whole_gene__ncbiRefSeq.bed") as m2:
      for line in m2:
         line = line.strip()
         if line:
            if not line.startswith("chr"):
               m1.write("{}\n".format(line))
            else:
               line_split = line.split("\t")
               # Select RefSeq transcript identifiers that start with NM_ (coding)
               if line_split[3].startswith("NM_"):
                  # https://github.com/BlankenbergLab/nAltORFs/blob/af1db6f23cda893a22d8c0b747c3c9af8dc68959/naltorfs/find_nested_alt_orfs.py#LL215C39-L215C39
                  # get_coding_exon_regions
                  chromStart = int(line_split[1])
                  blockStarts = [int(x) for x in line_split[11].split(",") if x.strip()]
                  blockSizes = [int(x) for x in line_split[10].split(",") if x.strip()]
                  thickStart = int(line_split[6])
                  thickEnd = int(line_split[7])
                  starts = list()
                  ends = list()
                  exon_starts = [x + chromStart for x in blockStarts]
                  exon_ends = [x + y for x, y in zip(exon_starts, blockSizes)]
                  for start, end in zip(exon_starts, exon_ends):
                     start = max(start, thickStart)
                     end = min(end, thickEnd)
                     if start < end:
                        starts.append(start)
                        ends.append(end)
                  coding_exon_regions = [(x, y) for x, y in zip(starts, ends)]
                  # Get the last (or first) two coding exon regions
                  strand = line_split[5]
                  if strand == "+":
                     selection = coding_exon_regions[-max_cds:]
                     reshape_selection = list()
                     for x, y in selection:
                        y = y + expand_utr
                        if region_max_size > 0:
                           x =  y - region_max_size
                        reshape_selection.append((x, y))
                  else:
                     selection = coding_exon_regions[:max_cds]
                     reshape_selection = list()
                     for x, y in selection:
                        x = x - expand_utr
                        if region_max_size > 0:
                           y =  x + region_max_size
                        reshape_selection.append((x, y))
                  if reshape_selection:
                     # Define the new line
                     newline = line_split[:6]
                     entry_id = newline[3]
                     if region_max_size > 0:
                        # Get the (last CDS + 100 UTR) that has been already limited to 250nt
                        x, y = reshape_selection[-1]
                        if (x, y) not in regions:
                           newline[1] = str(x)
                           newline[2] = str(y)
                           newline[3] = "{}__exon_CDS_{}_plus_{}_UTR".format(entry_id, 1, expand_utr)
                           m1.write("{}\n".format("\t".join(newline)))
                           regions.add((x, y))
                     else:
                        # Merge results
                        merged = (reshape_selection[0][0], reshape_selection[-1][1])
                        if merged not in regions:
                           newline[1] = str(merged[0])
                           newline[2] = str(merged[1])
                           newline[3] = "{}__exon_CDS_{}_plus_{}_UTR".format(entry_id, len(reshape_selection), expand_utr)
                           m1.write("{}\n".format("\t".join(newline)))
                           regions.add(merged)

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

# Apply the annotation to the selection
# Add the Ensembl transcript ID with the gene symbol
with open("selection_250nt_with_protein_coding_genes.bed", "w+") as m1:
   with open("selection_250nt.bed") as m2:
      for line in m2:
         line = line.strip()
         if line:
            if not line.startswith("chr"):
               m1.write("{}\n".format(line))
            else:
               line_split = line.split("\t")
               if line_split[0] in genes:
                  start = int(line_split[1])
                  end = int(line_split[2])
                  strand = line_split[5]
                  for gene in genes[line_split[0]]:
                     if genes[line_split[0]][gene]["start"] <= start and genes[line_split[0]][gene]["end"] >= end and genes[line_split[0]][gene]["strand"] == strand:
                        #line_split.append(gene)
                        line_split[3] = "{}__{}__{}".format(line_split[3].split("__")[0], gene, line_split[3].split("__")[1])
                        m1.write("{}\n".format("\t".join(line_split)))
                        break
