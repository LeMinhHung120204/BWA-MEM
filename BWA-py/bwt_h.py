from dataclasses import dataclass, field
from typing import List


# requirement: (OCC_INTERVAL%16 == 0); please DO NOT change this line because some part of the code assume OCC_INTERVAL=0x80
OCC_INTV_SHIFT = 7
OCC_INTERVAL = 1 << OCC_INTV_SHIFT
OCC_INTV_MASK = OCC_INTERVAL - 1

bwtint_t = int  

@dataclass
class bwt_t:
    primary: bwtint_t = 0  # S^{-1}(0), or the primary index of BWT
    L2: List[bwtint_t] = field(default_factory=lambda: [0] * 5)  # C(), cumulative count
    seq_len: bwtint_t = 0  # Sequence length
    bwt_size: bwtint_t = 0  # Size of BWT, about seq_len/4
    bwt: List[int] = field(default_factory=list)  # BWT array

    # Occurrence array, separated into two parts
    cnt_table: List[int] = field(default_factory=lambda: [0] * 256)

    # Suffix array
    sa_intv: int = 0
    n_sa: bwtint_t = 0
    sa: List[bwtint_t] = field(default_factory=list)

@dataclass
class bwtintv_t:
    x: List[bwtint_t] = field(default_factory=lambda: [0] * 3)  # Array of 3 elements
    info: bwtint_t = 0

@dataclass
class bwtintv_v:
    n: int = 0  # Current number of elements
    m: int = 0  # Maximum capacity (not strictly needed in Python)
    a: List[bwtintv_t] = field(default_factory=list)  # Dynamic list of intervals

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
