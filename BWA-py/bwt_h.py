# requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
OCC_INTV_SHIFT = 7
OCC_INTERVAL = 1 << OCC_INTV_SHIFT
OCC_INTV_MASK = OCC_INTERVAL - 1

class BWT:
    def __init__(self, primary, L2, seq_len, bwt_size, bwt, cnt_table, sa_intv, n_sa, sa):
        self.primary = primary  # S^{-1}(0), or the primary index of BWT
        self.L2 = L2  # C(), cumulative count
        self.seq_len = seq_len  # sequence length
        self.bwt_size = bwt_size  # size of bwt, about seq_len/4
        self.bwt = bwt  # BWT
        # occurrence array, separated into two parts
        self.cnt_table = cnt_table
        # suffix array
        self.sa_intv = sa_intv
        self.n_sa = n_sa
        self.sa = sa

class BWTInterval:
    def __init__(self, x, info):
        self.x = x  # x[0], x[1], x[2]
        self.info = info

class BWTIntervalVector:
    def __init__(self):
        self.n = 0
        self.m = 0
        self.a = []  # List of BWTInterval

# The following two functions are ONLY correct when OCC_INTERVAL == 0x80
def bwt_bwt(bwt, k):
    return bwt.bwt[((k >> 7) << 4) + 8 + ((k & 0x7f) >> 4)]

def bwt_occ_intv(bwt, k):
    return bwt.bwt + ((k >> 7) << 4)

def bwt_B0(bwt, k):
    return (bwt_bwt(bwt, k) >> ((~k & 0xf) << 1)) & 3

def bwt_set_intv(bwt, c, ik):
    ik.x[0] = bwt.L2[c] + 1
    ik.x[2] = bwt.L2[c + 1] - bwt.L2[c]
    ik.x[1] = bwt.L2[3 - c] + 1
    ik.info = 0
