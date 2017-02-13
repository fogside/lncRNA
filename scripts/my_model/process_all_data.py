import numpy as np
import pandas as pd
from utils import *
from optparse import OptionParser
import os

parser = OptionParser()

parser.add_option("-c", "--code_fasta", action="store", dest="code_fasta",
                  help="FASTA file for coding")
parser.add_option("-n", "--noncode_fasta", action="store", dest="noncode_fasta",
                  help="FASTA file for coding")
parser.add_option("--uniref_code", action="store", dest="code_uniref_raw",
                  help="uniref for coding")
parser.add_option("--uniref_noncode", action="store", dest="noncode_uniref_raw",
                  help="uniref for noncoding")
parser.add_option("--my_cod_file", action="store", dest="my_cod_file",
                  help="my feats for coding")
parser.add_option("--my_noncod_file", action="store", dest="my_noncod_file",
                  help="my feats for noncoding")

parser.add_option("--out", action="store", dest="out_file",
                  help="path for out file")
parser.add_option("--folder_code", action="store", dest="folder_code",
                  help="folder of coding hmm 3 frames result")
parser.add_option("--folder_noncode", action="store", dest="folder_noncode",
                  help="folder of noncoding hmm 3 frames result")

uniref_cols = ["queryId", "subjectId", "percIdentity", "alnLength",
               "mismatchCount", "gapOpenCount", "queryStart",
               "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore"]

## this is necessary to transform eValue into log10 and used to fill 0 values
NULL_CONST = -1000
(options, args) = parser.parse_args()


## Uniref ###
#############


def process_uniref(refs, uniref_raw, verbose=True):
    ## fasta names should be unique
    ## refs is a df with all original fasta names

    if verbose:
        print("Processing uniref....\nReading data....")

    data = pd.read_csv(uniref_raw)
    data.columns = uniref_cols
    data = data.ix[:, ["queryId", "percIdentity", "alnLength", "eVal", "bitScore"]]
    not_in_blast = list(set(refs.all_names) - set(data.queryId))

    if verbose:
        print("Raw_uniref len: ", len(data))
        print('Len of names in uniref: ', len(set(data.queryId)))
        print('Num names in reference: ', len(refs))
        print('Num of names not in uniref: ', len(not_in_blast))

    ## choose only names from ref and maybe add some missed names from ref
    data = pd.merge(refs, data, left_on='all_names', right_on='queryId', how='left')
    data.drop("queryId", axis=1, inplace=True)
    ## rename all_names to queryId
    data.columns = ["queryId", "percIdentity", "alnLength", "eVal", "bitScore"]
    data.eVal = [np.log10(i) if i != 0 else NULL_CONST for i in data.eVal]
    data = data.fillna({"percIdentity": 0, "alnLength": 0, "eVal": 0, "bitScore": 0}, inplace=True)

    ### try to choose only rows with the lowest eVal,
    ### but there're still a lot of duplicates..
    ### So we continued to check the other params to find the better combination
    ### for each queryId

    eVal = pd.DataFrame(data.groupby('queryId')['eVal'].min())
    eVal['queryId'] = eVal.index
    merge_eVal = pd.merge(data, eVal)

    bscore = pd.DataFrame(merge_eVal.groupby('queryId')['bitScore'].max())
    bscore['queryId'] = bscore.index
    merge_bscore = pd.merge(merge_eVal, bscore)
    del merge_eVal

    pid = pd.DataFrame(merge_bscore.groupby('queryId')['percIdentity'].max())
    pid['queryId'] = pid.index
    merge_pid = pd.merge(merge_bscore, pid)
    del merge_bscore

    length = pd.DataFrame(merge_pid.groupby('queryId')['alnLength'].max())
    length['queryId'] = length.index
    merge_len = pd.merge(merge_pid, length)
    del merge_pid

    hits = pd.DataFrame(data.groupby('queryId')['percIdentity'].count())
    hits['queryId'] = hits.index
    ## just rename column
    hits.columns = ['blast_hits_count', 'queryId']
    hits.ix[not_in_blast, 'blast_hits_count'] = 0
    merge_hits = pd.merge(merge_len, hits)
    del merge_len

    merge_hits.drop_duplicates(inplace=True)

    if verbose:
        print("Result uniref size without duplicates: ", len(merge_hits))

    return merge_hits


## Hmmer ##
###########

def process_hmmer(refs, folder, verbose=True):

    if verbose:
        print("Starting reading from hmmer folder....")

    file_names = []
    for (dirpath, dirnames, filenames) in os.walk(folder):
        file_names.extend([os.path.join(dirpath,f) for f in filenames])
        break

    files_df = []
    for f in file_names:
        files_df.append(process_hmm_file(f, NULL_CONST))

    data = pd.concat(files_df)
    data['queryId'] = data.index
    del files_df

    if verbose:
        print("HMMER: found hits for {0} queryIds in {1} files".format(len(set(data.queryId)), len(file_names)))
        print("Ref length: {0} files".format(len(refs)))

    data = pd.merge(refs, data, left_on='all_names', right_on='queryId', how='left')
    data.drop('queryId', axis=1, inplace=True)
    data.columns = ['queryId', 'eVal_hmm', 'num_hits_hmm', 'score_hmm']
    data.fillna(0, inplace=True)

    eVal = pd.DataFrame(data.groupby('queryId')['eVal_hmm'].min())
    eVal['queryId'] = eVal.index
    merge_eVal = pd.merge(data, eVal)

    score = pd.DataFrame(merge_eVal.groupby('queryId')['score_hmm'].max())
    score['queryId'] = score.index
    merge_score = pd.merge(merge_eVal, score)
    del merge_eVal

    hits = pd.DataFrame(merge_score.groupby('queryId')['num_hits_hmm'].max())
    hits['queryId'] = hits.index
    # eVal.head()
    merge_hits = pd.merge(merge_score, hits)
    del merge_score

    merge_hits.drop_duplicates(inplace=True)

    if verbose:
        print("Result Hmm size without duplicates: ", len(merge_hits))

    return merge_hits


#########################
## Making all the data ##
#########################

# code_fasta =
# noncode_fasta =
# code_uniref_raw =
# noncode_uniref_raw =
#
# my_cod_file =
# my_noncod_file =
#
# out_file =
# folder_code =
# folder_noncode =

names_cod = list(read_FASTA(options.code_fasta).keys())
refs_cod = pd.DataFrame({'all_names': names_cod}, index=np.arange(len(names_cod)))

names_noncod = list(read_FASTA(options.noncode_fasta).keys())
refs_noncod = pd.DataFrame({'all_names': names_noncod}, index=np.arange(len(names_noncod)))

uniref_cod_data = process_uniref(refs_cod, options.code_uniref_raw, verbose=True)
uniref_noncod_data = process_uniref(refs_noncod, options.noncode_uniref_raw, verbose=True)
uniref_cod_data['TYPE'] = [1] * len(uniref_cod_data)
uniref_noncod_data['TYPE'] = [0] * len(uniref_noncod_data)

hmm_code_data = process_hmmer(refs_cod, options.folder_code, verbose=True)
hmm_noncode_data = process_hmmer(refs_noncod, options.folder_noncode, verbose=True)

my_coding = pd.read_csv(options.my_cod_file, sep='\t')
my_noncoding = pd.read_csv(options.my_noncod_file, sep='\t')

uniref_all = pd.concat((uniref_cod_data, uniref_noncod_data))
hmm_all = pd.concat((hmm_code_data, hmm_noncode_data))
my_feats = pd.concat((my_coding, my_noncoding))

merged = pd.merge(hmm_all, uniref_all, on='queryId')
merged = pd.merge(my_feats, merged, right_on='queryId', left_on='sname')
merged.drop('queryId', axis=1, inplace=True)


print("Finished. Merged shape: {}".format(merged.shape))

merged.to_csv(options.out_file, index=False)
