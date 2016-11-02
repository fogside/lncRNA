import numpy as np


class ORFFinder:
    def __init__(self, seq):

        self.seq = seq.upper()

    def _find_start(self, seq, start=0):
        for i in range(start, len(seq), 3):
            if seq[i:i + 3] == 'ATG':
                return i
        return -1

    def _find_stop(self, seq, start=0, stop_codons=['TAG', 'TAA', 'TGA']):
        for i in range(start, len(seq), 3):
            if seq[i:i + 3] in stop_codons:
                return i
        return -1

    def find_longest(self):

        starts = []
        stops = []

        for i in range(3):
            ret1 = ret2 = 0
            results = []
            while ret2 != -1 and ret1 != -1:
                ret1 = self._find_start(self.seq, start=(i if ret2 == 0 else ret2 + 3))
                if ret1 == -1:
                    break
                ret2 = self._find_stop(self.seq, ret1 + 3)
                if ret2 == -1:
                    break

                results.append((ret2 - ret1, ret1, ret2))

            if results != []:
                max_size_idx = np.argmax([m[0] for m in results])
                starts.append(results[max_size_idx][1])
                stops.append(results[max_size_idx][2])
                # print('max_size:', results[max_size_idx][0])

        if len(starts) == 0:
            return -1

        long3 = [self.seq[n:k + 3] for n, k in zip(starts, stops)]
        longest = long3[np.argmax([len(s) for s in long3])]

        mean_orf_length = np.mean([len(seqv) for seqv in long3])
        longest_size = len(longest)
        orf_coverage = (longest_size*1.0)/len(self.seq)

        return starts, longest, longest_size, mean_orf_length, orf_coverage