import numpy as np


def clear_name(name):
    name = name.replace('"', '')
    name = name.replace("'", "")
    return name


def gtf_parser(file_name, drop_n_first_lines=1, f_format='gtf'):

    exon_length = {}

    with open(file_name, 'r') as f:
        gtf = f.readlines()

    gtf = gtf[drop_n_first_lines:]

    for line in gtf:

        tmp = line.split(';')
        if f_format == 'gff':
            name = tmp[2].split('=')[1].upper() if tmp[2].startswith('transcript_id') else None
        elif f_format == 'gtf':
            name = tmp[1].split(' ')[2].upper() if tmp[1].split(' ')[1].startswith('transcript_id') else None
        else:
            print("Incorrect f_format parameter. Please, choose `gtf` or `gff`")
            return -1

        tmp = tmp[0].split('\t')
        exon = int(tmp[4]) - int(tmp[3]) if tmp[2] == 'exon' else None

        if name is None:
            print("Incorrect data format. Use gtf file!")
            return -1
        elif exon is None:
            continue
        else:
            name = clear_name(name)
            if name in exon_length:
                exon_length[name].append(exon)
            else:
                exon_length[name] = [exon]

    res_max = {key: max(val) for key, val in exon_length.items()}
    res_mean = {key: np.mean(val) for key, val in exon_length.items()}
    res_num = {key: len(val) for key, val in exon_length.items()}
    return res_max, res_mean, res_num