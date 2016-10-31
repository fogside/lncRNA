import re
import numpy as np


class ORFFinder:
    def __init__(self, seq, start_codons=['ATG'], stop_codons=['TAG', 'TAA', 'TGA']):

        self.start_codons = start_codons
        self.stop_codons = stop_codons
        self.seq = seq.upper()

    def _pairwise_distances_(self, start_vec, stop_vec):
        '''
        start_vec, stop_vec -- numpy 1-d arrays with possibly different length;
        function returns the matrix of pairwise distances between elements of vectors;
        stop_vec -- expected to be the vector of stop_codons' locations.
        Distance is positive value. Negative (and zero) distance is replaced by np.inf.

        '''
        mtx = np.array(stop_vec) - np.array(start_vec).reshape((len(start_vec), 1))
        sh = mtx.shape
        mtx = np.array([i if i > 0 else np.inf for i in mtx.flatten()]).reshape(sh)
        return mtx

    def _return_start_stop_(self, start_vec, stop_vec):

        '''
        return positions of start and stop codons
        with max distance between them without any stop 

        '''

        dist_mtx = self._pairwise_distances_(start_vec, stop_vec)

        ## just find raw with max element among min elements of the each row ##
        indx_start = np.argmax([i if i is not np.inf else -1 for i in np.min(dist_mtx, axis=1)])

        ## find column with min value for the row number [indx_start]
        indx_stop = np.argmin(dist_mtx[indx_start])
        return start_vec[indx_start], stop_vec[indx_stop]

    def _find_three_starts_(self, patterns=['ATG']):

        addrs = {0: [], 1: [], 2: []}
        bag_of_start = []

        for patt in patterns:
            bag_of_start.extend([n.start() for n in re.finditer(patt, self.seq)])

        if len(bag_of_start) == 0:
            print("Pattern hasn't been found")
            return -1

        for i in range(3):
            lst = [j for j in bag_of_start if j % 3 == i]
            if len(lst) == 0:
                addrs[i] = -1
                continue
            addrs[i] = lst

        return addrs

    def find_three_longest(self):

        start_addrs = self._find_three_starts_(patterns=self.start_codons)
        stop_addrs = self._find_three_starts_(patterns=self.stop_codons)

        self.longests = []
        for start_vec, stop_vec in zip(start_addrs.values(), stop_addrs.values()):
            if (start_vec != -1) and (stop_vec != -1):
                start, stop = self._return_start_stop_(start_vec, stop_vec)
                self.longests.append(self.seq[start:stop + 3])

        return self.longests

    def longest_orf(self):

        find_max = lambda lst: lst[np.argmax([len(k) for k in lst])]
        return find_max(self.find_three_longest())
