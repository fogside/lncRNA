import numpy as np
import pandas as pd

def read_FASTA(file_name):
    with open(file_name, "r") as fn:
        text = fn.read().split(">")
    text = [x.split("\n") for x in text if x != ""]
    text = [[x[0], "".join(x[1:]).upper()] for x in text]
    text_dict = {line[0].split('|')[0]: line[1] for line in text}
    return text_dict


def process_hmm_file(filename, const_for_zeros=-1000):
    data_lines = []
    with open(filename, 'r') as hmm:
        data_lines = hmm.readlines()

    ## because there's a metadata on those lines
    ## which is unusable
    data_lines = data_lines[3:-10]

    def process_line(line):
        tmp = line.split()

        name = tmp[2].split('|')[0]
        eVal, score, bias = float(tmp[4]), float(tmp[5]), float(tmp[6])
        if eVal != 0.0:
            eVal = np.log10(eVal)
        else:
            eVal = const_for_zeros  ## this const is necessary for replacing true zero values

        if bias < score:
            return (name, eVal, score)
        else:
            return None

    seq_dict = {}  ## making dict for DataFrame creation
    for i, line in enumerate(data_lines):
        res = process_line(line)
        if res:  ## if it is not None then add
            if res[0] not in seq_dict.keys():
                seq_dict[res[0]] = {}
                seq_dict[res[0]]['eVal'] = [res[1]]
                seq_dict[res[0]]['score'] = [res[2]]
            else:
                seq_dict[res[0]]['eVal'].append(res[1])
                seq_dict[res[0]]['score'].append(res[2])

    for key in seq_dict.keys():
        seq_dict[key]['num_hits'] = len(seq_dict[key]['eVal'])
        tmp_arg = np.argmin(seq_dict[key]['eVal'])
        seq_dict[key]['eVal'] = seq_dict[key]['eVal'][tmp_arg]
        seq_dict[key]['score'] = seq_dict[key]['score'][tmp_arg]

    result_df = pd.DataFrame(seq_dict).T
    return result_df