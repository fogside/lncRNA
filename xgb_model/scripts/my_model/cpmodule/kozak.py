import numpy as np

def kozak_feat(seq, find_num):
    # first_feature {-3,+4}
    # second -- {-2, -1}
    # third -- {-6}
    # method returns 3 boolean values,
    # that mean find or not find for every feature


    if find_num < 6:
        return -1

    else:
        feat34 = 0
        feat21 = 0
        feat6 = 0

        num6 = seq[find_num - 6]
        if len(seq) >= find_num + 4:
            num4 = seq[find_num + 3]
        else:
            num4 = -1

        num3 = seq[find_num - 3]
        num2 = seq[find_num - 2]
        num1 = seq[find_num - 1]

    if num6 == "G":
        feat6 = 1
    if num4 == "G" and (num3 == "A" or num3 == "G"):
        feat34 = 1
    if num2 == "C" and num1 == "C":
        feat21 = 1

    return feat34, feat21, feat6


def find_kozak_feat(seq, atg_addr):

    atg_feat = []
    feat_all = False

    if len(atg_addr) < 1:
        return (0, 0, 0)

    for atg in atg_addr:
        feat_all = kozak_feat(seq, atg)
        atg_feat.append(feat_all)
        if feat_all == -1:
            return (0, 0, 0)
        if feat_all == (1, 1, 1):
            feat_all = 1
            break

    if feat_all == 1:
        return (1, 1, 1)
    else:
        return atg_feat[np.argmax([sum(i) for i in atg_feat])]







