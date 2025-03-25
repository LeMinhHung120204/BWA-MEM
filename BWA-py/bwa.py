import os
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

def bwa_idx_infer_prefix(hint: str) -> str | None:
    l_hint = len(hint)
    
    # Kiểm tra file với đuôi ".64.bwt"
    prefix = hint + ".64.bwt"
    if os.path.isfile(prefix):
        return hint + ".64"  # Cắt bỏ phần ".bwt"

    # Kiểm tra file với đuôi ".bwt"
    prefix = hint + ".bwt"
    if os.path.isfile(prefix):
        return hint  # Giữ nguyên hint nếu tồn tại file .bwt

    return None  # Không tìm thấy file nào

def bwa_idx_load_bwt(hint: str):
    prefix = bwa_idx_infer_prefix(hint)
    if prefix is None:
        if bwa_verbose >= 1:
            print(f"[E::bwa_idx_load_bwt] fail to locate the index files")
        return None

    # Khởi tạo bwt object
    bwt = None

    # Khôi phục FM-index (BWT)
    bwt_file = prefix + ".bwt"
    if os.path.isfile(bwt_file):
        bwt = bwt_restore_bwt(bwt_file)
    
    # Khôi phục partial suffix array (SA)
    sa_file = prefix + ".sa"
    if os.path.isfile(sa_file) and bwt is not None:
        bwt_restore_sa(sa_file, bwt)
    
    return bwt

def bwa_idx_load_from_disk(hint: str, which: int) -> BwaIndex:
    prefix = bwa_idx_infer_prefix(hint)
    if prefix is None:
        print(f"[E::{bwa_idx_load_from_disk.__name__}] fail to locate the index files")
        return None

    idx = bwaidx_t()

    if which & BWA_IDX_BWT:
        idx.bwt = bwa_idx_load_bwt(hint)

    if which & BWA_IDX_BNS:
        idx.bns = bns_restore(prefix)

        alt_count = sum(1 for ann in idx.bns.anns if ann.is_alt)
        if bwa_verbose >= 3:
            print(f"[M::{bwa_idx_load_from_disk.__name__}] read {alt_count} ALT contigs")

        if which & BWA_IDX_PAC:
            idx.pac = bytearray(idx.bns.l_pac // 4 + 1)
            idx.bns.fp_pac.readinto(idx.pac)  # Đọc dữ liệu 2-bit từ file
            idx.bns.fp_pac.close()
            idx.bns.fp_pac = None  # Giải phóng file sau khi đọc

    return idx

def bwa_mem2idx(l_mem, mem, idx : bwaidx_t):
    k = 0

    # Generate idx.bwt
    idx.bwt = bwt_t()
    idx.bwt.bwt_size = int.from_bytes(mem[k:k+8], "little"); k += 8
    idx.bwt.n_sa = int.from_bytes(mem[k:k+8], "little"); k += 8

    x = idx.bwt.bwt_size * 4
    idx.bwt.bwt = list(mem[k:k+x]); k += x

    x = idx.bwt.n_sa * 8  # Giả sử bwtint_t là uint64_t (8 byte)
    idx.bwt.sa = list(mem[k:k+x]); k += x

    # Generate idx.bns and idx.pac
    idx.bns = bntseq_t()
    idx.bns.n_holes = int.from_bytes(mem[k:k+4], "little"); k += 4
    idx.bns.n_seqs = int.from_bytes(mem[k:k+4], "little"); k += 4
    idx.bns.l_pac = int.from_bytes(mem[k:k+8], "little"); k += 8

    x = idx.bns.n_holes * 8  # Giả sử bntamb1_t là 8 byte
    idx.bns.ambs = list(mem[k:k+x]); k += x

    x = idx.bns.n_seqs * 16  # Giả sử bntann1_t là 16 byte
    idx.bns.anns = list(mem[k:k+x]); k += x

    for i in range(idx.bns.n_seqs):
        name_len = mem[k:].index(0)  # Tìm null-terminated string
        idx.bns.anns[i] = {"name": mem[k:k+name_len].decode()}; k += name_len + 1

        anno_len = mem[k:].index(0)
        idx.bns.anns[i]["anno"] = mem[k:k+anno_len].decode(); k += anno_len + 1

    idx.pac = list(mem[k:k + (idx.bns.l_pac // 4) + 1])
    k += len(idx.pac)

    assert k == l_mem

    idx.l_mem = k
    idx.mem = mem
    return 0

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