rmvar_clinvar_rs_ids = list()

with open("/data/home/cumbof/isilon/cumbof/m6A/RMVar/RMVar_Human_Clinvar_info.txt") as m:
    header = m.readline().strip().split("\t")
    for line in m:
        line = line.strip()
        if line:
            rmvar_clinvar_rs_ids.append(line.split("\t")[1][2:])

rmvar_clinvar_rs_ids = list(set(rmvar_clinvar_rs_ids))

len(rmvar_clinvar_rs_ids)

rmvar_clinvar_rs_ids_m6A = list()

with open("/data/home/cumbof/isilon/cumbof/m6A/RMVar/RMVar_Human_basic_info_m6A.txt") as m:
    header = m.readline().strip().split("\t")
    for line in m:
        line = line.strip()
        if line:
            line_split = line.split("\t")
            rs_id = line_split[header.index("rs_id")][2:]
            if rs_id in rmvar_clinvar_rs_ids:
                rmvar_clinvar_rs_ids_m6A.append(rs_id)

rmvar_clinvar_rs_ids_m6A = list(set(rmvar_clinvar_rs_ids_m6A))

len(rmvar_clinvar_rs_ids_m6A)

rmvar_clinvar_rs_ids_m6A = {
    rs_id: 0 for rs_id in rmvar_clinvar_rs_ids_m6A
}

with open("/data/home/cumbof/isilon/cumbof/m6A/General/clinvar_m6A.vcf", "w+") as clinvar_m6A:
    with open("/data/home/cumbof/isilon/cumbof/m6A/General/clinvar.vcf") as clinvar:
        for line in clinvar:
            if line.strip():
                if line.startswith("#"):
                    clinvar_m6A.write(line)
                else:
                    line_split = line.strip().split("\t")
                    info = line_split[7].split(";")
                    for i in info:
                        if i.startswith("RS="):
                            rs_id = i[3:]
                            if rs_id in rmvar_clinvar_rs_ids_m6A:
                                clinvar_m6A.write(line)
                                break
