biglist = list()
smalllist = list()

with open("/Users/fabio/GitHub/m6A/General/targets-ClinVar-m6A-gencode-20230728/targets.tsv") as m:
    header = m.readline().strip().split("\t")
    for line in m:
        line = line.strip()
        if line:
            line_split = line.split("\t")
            biglist.append((line_split[header.index("entry_id")], line_split[header.index("transcript_id")]))

with open("/Users/fabio/GitHub/m6A/General/targets-ClinVar-m6A-gencode-20230728/targets-250nt.tsv") as m:
    header = m.readline().strip().split("\t")
    for line in m:
        line = line.strip()
        if line:
            line_split = line.split("\t")
            smalllist.append((line_split[header.index("entry_id")], line_split[header.index("transcript_id")]))

inter = list(set(biglist).intersection(set(smalllist)))

with open("/Users/fabio/GitHub/m6A/General/targets-ClinVar-m6A-gencode-20230728/index/index.html") as m:
    for line in m:
        line = line.strip()
        if line:
            for couple in inter:
                if ">{}</a>".format(couple[0]) in line and "<td>{}</td>".format(couple[1]) in line:
                    line_split = line.split("<td>")
                    line_split.insert(2, "no</td>")
                    print("<td>".join(line_split))
                    #newline = "<tr><td>{}</td><td>{}</td><td>no</td><td>".format(couple[0], couple[1]) + "<td>".join(line_split[-12:])
                    #print(newline)
                    break

with open("/Users/fabio/GitHub/m6A/General/targets-ClinVar-m6A-gencode-20230728/index/index-250nt.html") as m:
    for line in m:
        line = line.strip()
        if line:
            for couple in inter:
                if ">{}</a>".format(couple[0]) in line and "<td>{}</td>".format(couple[1]) in line:
                    line_split = line.split("<td>")
                    line_split.insert(2, "yes</td>")
                    print("<td>".join(line_split))
                    #newline = "<tr><td>{}</td><td>{}</td><td>yes</td><td>".format(couple[0], couple[1]) + "<td>".join(line_split[-12:])
                    #print(newline)
                    break

