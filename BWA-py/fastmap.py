import sys
import getopt
import math
import string

from bwa import *
from bwamem import *

class KtpAux:
    def __init__(self):
        self.ks = None
        self.ks2 = None
        self.opt = None
        self.pes0 = None
        self.n_processed = 0
        self.copy_comment = 0
        self.actual_chunk_size = 0
        self.idx = None

class KtpData:
    def __init__(self, aux, n_seqs, seqs):
        self.aux = aux
        self.n_seqs = n_seqs
        self.seqs = seqs

def process(shared, step, _data):
    if step == 0:
        size = 0
        seqs, n_seqs = bseq_read(shared.actual_chunk_size, shared.ks, shared.ks2)
        if not seqs:
            return None
        
        if not shared.copy_comment:
            for seq in seqs:
                seq.comment = None
                
        for seq in seqs:
            size += len(seq.seq)
        
        if bwa_verbose >= 3:
            print("[M::process] read {} sequences ({} bp)...".format(n_seqs, size), file=sys.stderr)
        
        return KtpData(shared, n_seqs, seqs)
    
    elif step == 1:
        opt = shared.opt
        idx = shared.idx
        
        if opt.flag & MEM_F_SMARTPE:
            sep = [[], []]
            n_sep = [0, 0]
            tmp_opt = opt
            
            bseq_classify(_data.n_seqs, _data.seqs, n_sep, sep)
            
            if bwa_verbose >= 3:
                print("[M::process] {} single-end sequences; {} paired-end sequences".format(n_sep[0], n_sep[1]), file=sys.stderr)
            
            if n_sep[0]:
                tmp_opt.flag &= ~MEM_F_PE
                mem_process_seqs(tmp_opt, idx.bwt, idx.bns, idx.pac, shared.n_processed, n_sep[0], sep[0], None)
                for seq in sep[0]:
                    _data.seqs[seq.id].sam = seq.sam
            
            if n_sep[1]:
                tmp_opt.flag |= MEM_F_PE
                mem_process_seqs(tmp_opt, idx.bwt, idx.bns, idx.pac, shared.n_processed + n_sep[0], n_sep[1], sep[1], shared.pes0)
                for seq in sep[1]:
                    _data.seqs[seq.id].sam = seq.sam
        
        else:
            mem_process_seqs(opt, idx.bwt, idx.bns, idx.pac, shared.n_processed, _data.n_seqs, _data.seqs, shared.pes0)
        
        shared.n_processed += _data.n_seqs
        return _data
    
    elif step == 2:
        for seq in _data.seqs:
            if seq.sam:
                sys.stdout.write(seq.sam)
            
            del seq.name, seq.comment, seq.seq, seq.qual, seq.sam
        
        return None
    
    return None

def update_a(opt, opt0):
    if opt0.a:  # matching score is changed
        if not opt0.b: opt.b *= opt.a
        if not opt0.T: opt.T *= opt.a
        if not opt0.o_del: opt.o_del *= opt.a
        if not opt0.e_del: opt.e_del *= opt.a
        if not opt0.o_ins: opt.o_ins *= opt.a
        if not opt0.e_ins: opt.e_ins *= opt.a
        if not opt0.zdrop: opt.zdrop *= opt.a
        if not opt0.pen_clip5: opt.pen_clip5 *= opt.a
        if not opt0.pen_clip3: opt.pen_clip3 *= opt.a
        if not opt0.pen_unpaired: opt.pen_unpaired *= opt.a

def main_mem(argv):
    opt = MemOpt()
    opt0 = MemOpt()
    ignore_alt = 0
    no_mt_io = 0
    fixed_chunk_size = -1
    fp = None
    fp2 = None
    rg_line = None
    hdr_line = None
    mode = None
    ko = None
    ko2 = None

    # Khởi tạo mảng pes gồm 4 phần tử kiểu MemPestat
    pes = [MemPestat() for _ in range(4)]

    # Tạo thể hiện của KtpAux và gán giá trị
    aux = KtpAux()
    
    aux.opt = opt = mem_opt_init()
    opts, args = getopt.getopt(sys.argv[1:], "51qpaMCSPVYjuk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:o:f:W:x:G:h:y:K:X:H:F:z:")

    for c, arg in opts:
        if c == '-k':
            opt.min_seed_len = int(arg)
            opt0.min_seed_len = 1
        elif c == '-1':
            no_mt_io = 1
        elif c == '-x':
            mode = arg
        elif c == '-w':
            opt.w = int(arg)
            opt0.w = 1
        elif c == '-A':
            opt.a = int(arg)
            opt0.a = 1
        elif c == '-B':
            opt.b = int(arg)
            opt0.b = 1
        elif c == '-T':
            opt.T = int(arg)
            opt0.T = 1
        elif c == '-U':
            opt.pen_unpaired = int(arg)
            opt0.pen_unpaired = 1
        elif c == '-t':
            opt.n_threads = max(int(arg), 1)
        elif c == '-P':
            opt.flag |= (1 << 0)  # MEM_F_NOPAIRING
        elif c == '-a':
            opt.flag |= (1 << 1)  # MEM_F_ALL
        elif c == '-p':
            opt.flag |= (1 << 2) | (1 << 3)  # MEM_F_PE | MEM_F_SMARTPE
        elif c == '-M':
            opt.flag |= (1 << 4)  # MEM_F_NO_MULTI
        elif c == '-S':
            opt.flag |= (1 << 5)  # MEM_F_NO_RESCUE
        elif c == '-Y':
            opt.flag |= (1 << 6)  # MEM_F_SOFTCLIP
        elif c == '-V':
            opt.flag |= (1 << 7)  # MEM_F_REF_HDR
        elif c == '-5':
            opt.flag |= (1 << 8) | (1 << 9)  # MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ
        elif c == '-q':
            opt.flag |= (1 << 9)  # MEM_F_KEEP_SUPP_MAPQ
        elif c == '-u':
            opt.flag |= (1 << 10)  # MEM_F_XB
        elif c == '-c':
            opt.max_occ = int(arg)
            opt0.max_occ = 1
        elif c == '-d':
            opt.zdrop = int(arg)
            opt0.zdrop = 1
        elif c == '-v':
            bwa_verbose = int(arg)
        elif c == '-j':
            ignore_alt = 1
        elif c == '-r':
            opt.split_factor = float(arg)
            opt0.split_factor = 1.0
        elif c == '-D':
            opt.drop_ratio = float(arg)
            opt0.drop_ratio = 1.0
        elif c == '-m':
            opt.max_matesw = int(arg)
            opt0.max_matesw = 1
        elif c == '-s':
            opt.split_width = int(arg)
            opt0.split_width = 1
        elif c == '-G':
            opt.max_chain_gap = int(arg)
            opt0.max_chain_gap = 1
        elif c == '-N':
            opt.max_chain_extend = int(arg)
            opt0.max_chain_extend = 1
        elif c in ('-o', '-f'):
            sys.stdout = open(arg, "wb")
        elif c == '-W':
            opt.min_chain_weight = int(arg)
            opt0.min_chain_weight = 1
        elif c == '-y':
            opt.max_mem_intv = int(arg)
            opt0.max_mem_intv = 1
        elif c == '-C':
            aux.copy_comment = 1
        elif c == '-K':
            fixed_chunk_size = int(arg)
        elif c == '-X':
            opt.mask_level = float(arg)
        elif c == '-F':
            bwa_dbg = int(arg)
        elif c == '-h':
            opt0.max_XA_hits = opt0.max_XA_hits_alt = 1
            parts = arg.split()
            opt.max_XA_hits = opt.max_XA_hits_alt = int(parts[0])
            if len(parts) > 1 and parts[1][0] in string.punctuation and parts[1][1:].isdigit():
                opt.max_XA_hits_alt = int(parts[1][1:])
        elif c == '-z':
            opt.XA_drop_ratio = float(arg)
        elif c == '-Q':
            opt0.mapQ_coef_len = 1
            opt.mapQ_coef_len = int(arg)
            opt.mapQ_coef_fac = math.log(opt.mapQ_coef_len) if opt.mapQ_coef_len > 0 else 0
        elif c == '-O':
            opt0.o_del = opt0.o_ins = 1
            parts = arg.split()
            opt.o_del = opt.o_ins = int(parts[0])
            if len(parts) > 1 and parts[1][0] in string.punctuation and parts[1][1:].isdigit():
                opt.o_ins = int(parts[1][1:])
        elif c == '-E':
            opt0.e_del = opt0.e_ins = 1
            parts = arg.split()
            opt.e_del = opt.e_ins = int(parts[0])
            if len(parts) > 1 and parts[1][0] in string.punctuation and parts[1][1:].isdigit():
                opt.e_ins = int(parts[1][1:])
        elif c == '-L':
            opt0.pen_clip5 = opt0.pen_clip3 = 1
            parts = arg.split()
            opt.pen_clip5 = opt.pen_clip3 = int(parts[0])
            if len(parts) > 1 and parts[1][0] in string.punctuation and parts[1][1:].isdigit():
                opt.pen_clip3 = int(parts[1][1:])
        elif c == '-R':
            rg_line = arg  # Giữ nguyên chuỗi dòng read group
        elif c == '-H':
            if not arg.startswith('@'):
                try:
                    with open(arg, "r") as fp:
                        for line in fp:
                            line = line.strip()
                            if line:
                                hdr_line = line  # Giữ nguyên nội dung header
                except FileNotFoundError:
                    print(f"File {arg} not found!", file=sys.stderr)
            else:
                hdr_line = arg
        elif c == '-I':  # Specify the insert size distribution
            pes[1].failed = 0
            parts = arg.split()
            pes[1].avg = float(parts[0])
            pes[1].std = pes[1].avg * 0.1
            if len(parts) > 1 and parts[1][0] in string.punctuation and parts[1][1:].isdigit():
                pes[1].std = float(parts[1][1:])
            pes[1].high = int(pes[1].avg + 4.0 * pes[1].std + 0.499)
            pes[1].low = int(pes[1].avg - 4.0 * pes[1].std + 0.499)
            if pes[1].low < 1:
                pes[1].low = 1
            if len(parts) > 2 and parts[2][0] in string.punctuation and parts[2][1:].isdigit():
                pes[1].high = int(float(parts[2][1:]) + 0.499)
            if len(parts) > 3 and parts[3][0] in string.punctuation and parts[3][1:].isdigit():
                pes[1].low = int(float(parts[3][1:]) + 0.499)
            if bwa_verbose >= 3:
                print(f"[M::parse_args] mean insert size: {pes[1].avg:.3f}, stddev: {pes[1].std:.3f}, max: {pes[1].high}, min: {pes[1].low}")
        else:
            sys.exit(1)

    if rg_line:
        hdr_line = rg_line if hdr_line is None else f"{hdr_line}\n{rg_line}"
        rg_line = None  

    if opt.n_threads < 1:
        opt.n_threads = 1
    
    if len(args) < 1 or len(args) > 2:
        print("\nUsage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n", file=sys.stderr)
        print("Algorithm options:\n", file=sys.stderr)
        print(f"       -t INT        number of threads [{opt.n_threads}]", file=sys.stderr)
        print(f"       -k INT        minimum seed length [{opt.min_seed_len}]", file=sys.stderr)
        print(f"       -w INT        band width for banded alignment [{opt.w}]", file=sys.stderr)
        print(f"       -d INT        off-diagonal X-dropoff [{opt.zdrop}]", file=sys.stderr)
        print(f"       -r FLOAT      look for internal seeds inside a seed longer than {opt.min_seed_len} * FLOAT [{opt.split_factor}]", file=sys.stderr)
        print(f"       -y INT        seed occurrence for the 3rd round seeding [{opt.max_mem_intv}]", file=sys.stderr)
        print(f"       -c INT        skip seeds with more than INT occurrences [{opt.max_occ}]", file=sys.stderr)
        print(f"       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [{opt.drop_ratio:.2f}]", file=sys.stderr)
        print("       -W INT        discard a chain if seeded bases shorter than INT [0]", file=sys.stderr)
        print(f"       -m INT        perform at most INT rounds of mate rescues for each read [{opt.max_matesw}]", file=sys.stderr)
        print("       -S            skip mate rescue", file=sys.stderr)
        print("       -P            skip pairing; mate rescue performed unless -S also in use\n", file=sys.stderr)

        print("Scoring options:\n", file=sys.stderr)
        print(f"       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [{opt.a}]", file=sys.stderr)
        print(f"       -B INT        penalty for a mismatch [{opt.b}]", file=sys.stderr)
        print(f"       -O INT[,INT]  gap open penalties for deletions and insertions [{opt.o_del},{opt.o_ins}]", file=sys.stderr)
        print(f"       -E INT[,INT]  gap extension penalty; a gap of size k costs '{{-O}} + {{-E}}*k' [{opt.e_del},{opt.e_ins}]", file=sys.stderr)
        print(f"       -L INT[,INT]  penalty for 5'- and 3'-end clipping [{opt.pen_clip5},{opt.pen_clip3}]", file=sys.stderr)
        print(f"       -U INT        penalty for an unpaired read pair [{opt.pen_unpaired}]\n", file=sys.stderr)

        print("Input/output options:\n", file=sys.stderr)
        print("       -p            smart pairing (ignoring in2.fq)", file=sys.stderr)
        print("       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]", file=sys.stderr)
        print("       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]", file=sys.stderr)
        print("       -o FILE       sam file to output results to [stdout]", file=sys.stderr)
        print("       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)", file=sys.stderr)
        print("       -5            for split alignment, take the alignment with the smallest query (not genomic) coordinate as primary", file=sys.stderr)
        print("       -q            don't modify mapQ of supplementary alignments", file=sys.stderr)
        print(f"       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [{bwa_verbose}]", file=sys.stderr)
        print(f"       -T INT        minimum score to output [{opt.T}]", file=sys.stderr)
        print(f"       -h INT[,INT]  if there are <INT hits with score >{opt.XA_drop_ratio * 100:.2f}% of the max score, output all in XA [{opt.max_XA_hits},{opt.max_XA_hits_alt}]", file=sys.stderr)
        print("       -a            output all alignments for SE or unpaired PE", file=sys.stderr)
        print("       -C            append FASTA/FASTQ comment to SAM output", file=sys.stderr)
        print("       -V            output the reference FASTA header in the XR tag", file=sys.stderr)
        print("       -Y            use soft clipping for supplementary alignments", file=sys.stderr)
        print("       -M            mark shorter split hits as secondary\n", file=sys.stderr)

        print("\nNote: Please read the man page for detailed description of the command line and options.\n", file=sys.stderr)

        sys.exit(1)
    if mode:
        if mode == "intractg":
            if not opt0.o_del: opt.o_del = 16
            if not opt0.o_ins: opt.o_ins = 16
            if not opt0.b: opt.b = 9
            if not opt0.pen_clip5: opt.pen_clip5 = 5
            if not opt0.pen_clip3: opt.pen_clip3 = 5
        elif mode in ["pacbio", "pbref", "ont2d"]:
            if not opt0.o_del: opt.o_del = 1
            if not opt0.e_del: opt.e_del = 1
            if not opt0.o_ins: opt.o_ins = 1
            if not opt0.e_ins: opt.e_ins = 1
            if not opt0.b: opt.b = 1
            if opt0.split_factor == 0.0: opt.split_factor = 10.0

            if mode == "ont2d":
                if not opt0.min_chain_weight: opt.min_chain_weight = 20
                if not opt0.min_seed_len: opt.min_seed_len = 14
                if not opt0.pen_clip5: opt.pen_clip5 = 0
                if not opt0.pen_clip3: opt.pen_clip3 = 0
            else:
                if not opt0.min_chain_weight: opt.min_chain_weight = 40
                if not opt0.min_seed_len: opt.min_seed_len = 17
                if not opt0.pen_clip5: opt.pen_clip5 = 0
                if not opt0.pen_clip3: opt.pen_clip3 = 0
        else:
            print(f"[E::{sys._getframe().f_code.co_name}] unknown read type '{mode}'", file=sys.stderr)
            sys.exit(1)  # FIXME: Cảnh báo về khả năng rò rỉ bộ nhớ trong C không cần thiết trong Python
    else:
        update_a(opt, opt0)

    # đang làm

if __name__ == "__main__":
    pass