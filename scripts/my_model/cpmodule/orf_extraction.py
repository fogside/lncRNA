import re
import numpy as np


class ORFFinder:
    def __init__(self, seq, start_codons=['ATG'], stop_codons=['TAG', 'TAA', 'TGA']):

        self.start_codons = start_codons
        self.stop_codons = stop_codons
        self.longests = None
        self.seq = seq.upper()

    def _find_three_starts_(self, patterns=['ATG'], fun=np.min):

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
            addrs[i] = fun(lst)

        return addrs

    def find_three_longest(self):

        start_addrs = self._find_three_starts_(patterns=self.start_codons)
        stop_addrs = self._find_three_starts_(patterns=self.stop_codons, fun=np.max)

        if start_addrs == -1 or stop_addrs == -1:
            self.longests = -1
            return -1

        self.longests = []
        for start, stop in zip(start_addrs.values(), stop_addrs.values()):
            if (start >= 0) and (stop >= 0):
                self.longests.append(self.seq[start:stop + 3])

        return self.longests, start_addrs.values()

    def longest_orf(self):

        find_max = lambda lst: lst[np.argmax([len(k) for k in lst])] if lst!=-1 else -1


        if self.longests is not None:
            longest = find_max(self.longests)

        else:
            longest = find_max(self.find_three_longest()[0])

        if longest == -1:
            size = -1
        else:
            size = len(longest)

        return longest, size
