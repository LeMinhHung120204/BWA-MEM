import struct
from dataclasses import dataclass, field
from typing import List

from kstring import *
from bwt import *
from bntseq import *

BWA_IDX_BWT = 0x1
BWA_IDX_BNS = 0x2
BWA_IDX_PAC = 0x4
BWA_IDX_ALL = 0x7

BWA_CTL_SIZE = 0x10000

BWTALGO_AUTO = 0
BWTALGO_RB2 = 1
BWTALGO_BWTSW = 2
BWTALGO_IS = 3

BWA_DBG_QNAME = 0x1

@dataclass
class bwaidx_t:
    bwt: bwt_t = field(default_factory=bwt_t)  # FM-index từ Bwt
    bns: bntseq_t = field(default_factory=bntseq_t)  # Thông tin về reference sequences
    pac: List[int] = field(default_factory=list)  # 2-bit encoded reference sequences

    is_shm: int = 0
    l_mem: bwtint_t = 0
    mem: List[int] = field(default_factory=list)  # Lưu dữ liệu memory

class bseq1_t:
    def __init__(self):
        self.l_seq = 0
        self.id = 0
        self.name = None
        self.comment = None
        self.seq = None
        self.qual = None
        self.sam = None

bwa_verbose = 3
bwa_dbg = 0
bwa_rg_id = "" * 256
bwa_pg = None

#===============================================================================
# Batch FASTA/Q reader

def trim_readno(s: kstring_t):
    if s.l > 2 and s.s[s.l - 2] == '/' and s.s[s.l - 1].isdigit():
        s.l -= 2
        s.s = s.s[:s.l]

def dupkstring(str_obj: kstring_t, dupempty: int) -> str:
    if str_obj.l > 0 or dupempty:
        return str_obj.s[:str_obj.l] if str_obj.s else ""
    return None

def kseq2bseq1(ks, s: bseq1_t):
    s.name = dupkstring(ks.name, 1)
    s.comment = dupkstring(ks.comment, 0)
    s.seq = dupkstring(ks.seq, 1)
    s.qual = dupkstring(ks.qual, 0)
    s.l_seq = ks.seq.l

def bseq_read(chunk_size, n_, ks1, ks2=None):
    ks = ks1
    size = 0
    seqs = []
    n = 0
    while ks.read() >= 0:
        if ks2 and ks2.read() < 0:
            print("[W::bseq_read] the 2nd file has fewer sequences.")
            break
        trim_readno(ks.name)
        seq = bseq1_t()
        kseq2bseq1(ks, seq)
        seq.id = n
        seqs.append(seq)
        size += seq.l_seq
        n += 1
        if ks2:
            trim_readno(ks2.name)
            seq2 = bseq1_t()
            kseq2bseq1(ks2, seq2)
            seq2.id = n
            seqs.append(seq2)
            size += seq2.l_seq
            n += 1
        if size >= chunk_size and (n & 1) == 0:
            break
    n_[0] = n
    return seqs

def bseq_classify(n, seqs):
    a0, a1 = [], []
    has_last = True
    for i in range(1, n):
        if has_last:
            if seqs[i].name == seqs[i-1].name:
                a1.extend([seqs[i-1], seqs[i]])
                has_last = False
            else:
                a0.append(seqs[i-1])
        else:
            has_last = True
    if has_last:
        a0.append(seqs[-1])
    return {0: a0, 1: a1}

#===============================================================================
# CIGAR related

def bwa_fill_scmat(a, b, mat):
    k = 0
    for i in range(4):
        for j in range(4):
            mat[k] = a if i == j else -b
            k += 1
        mat[k] = -1  # ambiguous base
        k += 1
    for j in range(5):
        mat[k] = -1
        k += 1

#===============================================================================
# Full index reader

def bwa_mem2idx(l_mem, mem):
    idx = bwaidx_t()
    k = 0

    # Tạo idx.bwt
    x = struct.calcsize("2q")  # Giả định struct bwt_t có 2 int64_t
    idx.bwt = bwt_t(*struct.unpack_from("2q", mem, k))
    k += x
    x = idx.bwt.bwt_size * 4
    idx.bwt.bwt = mem[k:k + x]
    k += x
    x = idx.bwt.n_sa * 8  # bwtint_t là int64_t
    idx.bwt.sa = mem[k:k + x]
    k += x

    # Tạo idx.bns
    x = struct.calcsize("3q")  # Giả định struct bntseq_t có 3 int64_t
    idx.bns = bseq1_t(*struct.unpack_from("3q", mem, k))
    k += x
    x = idx.bns.n_holes * 8  # bntamb1_t giả định là int64_t
    idx.bns.ambs = mem[k:k + x]
    k += x
    x = idx.bns.n_seqs * struct.calcsize("2q")  # Giả định struct bntann1_t có 2 int64_t
    idx.bns.anns = [bntann1_t() for _ in range(idx.bns.n_seqs)]
    k += x

    # Đọc tên và chú thích
    for ann in idx.bns.anns:
        end = mem.index(b'\0', k)
        ann.name = mem[k:end].decode()
        k = end + 1
        end = mem.index(b'\0', k)
        ann.anno = mem[k:end].decode()
        k = end + 1

    # Tạo idx.pac
    idx.pac = mem[k:k + idx.bns.l_pac // 4 + 1]
    k += len(idx.pac)

    assert k == l_mem, "Lỗi: không khớp kích thước bộ nhớ"

    idx.l_mem = k
    idx.mem = mem
    return idx

#===============================================================================
# SAM header routines

def bwa_escape(s: str) -> str:
    p = 0
    q = []
    while p < len(s):
        if s[p] == "\\":
            p += 1
            if p < len(s):
                if s[p] == 't':
                    q.append('\t')
                elif s[p] == 'n':
                    q.append('\n')
                elif s[p] == 'r':
                    q.append('\r')
                elif s[p] == '\\':
                    q.append('\\')
        else:
            q.append(s[p])
        p += 1
    return "".join(q)

def bwa_insert_header(s: str, hdr: str = None) -> str:
    if not s or s[0] != '@':
        return hdr
    if hdr:
        hdr += '\n' + s
    else:
        hdr = s
    return bwa_escape(hdr)

if __name__ == "__main__":
    pass