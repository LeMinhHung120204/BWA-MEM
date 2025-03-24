import sys
import math
#from bwt import *  đang lỗi
from ksw import *
from bntseq import *
from bwa import * # bwa_verbose trong file bwa.c chưa hiểu cách xử lý 
from kvec import *
from bwt import *


MEM_MAPQ_COEF = 30.0
MEM_MAPQ_MAX = 60

MEM_F_PE = 0x2
MEM_F_NOPAIRING = 0x4
MEM_F_ALL = 0x8
MEM_F_NO_MULTI = 0x10
MEM_F_NO_RESCUE = 0x20
MEM_F_REF_HDR = 0x100
MEM_F_SOFTCLIP = 0x200
MEM_F_SMARTPE = 0x400
MEM_F_PRIMARY5 = 0x800
MEM_F_KEEP_SUPP_MAPQ = 0x1000
MEM_F_XB = 0x2000

class MemAlnReg:
    def __init__(self):
        self.rb = self.re = 0
        self.qb = self.qe = 0
        self.rid = 0
        self.score = 0
        self.truesc = 0
        self.sub = 0
        self.alt_sc = 0
        self.csub = 0
        self.sub_n = 0
        self.w = 0
        self.seedcov = 0
        self.secondary = 0
        self.secondary_all = 0
        self.seedlen0 = 0
        self.n_comp = 0
        self.is_alt = 0
        self.frac_rep = 0.0
        self.hash = 0

class MemAlnRegV:
    def __init__(self):
        self.n = self.m = 0
        self.a = []

class MemPestat:
    def __init__(self):
        self.low = self.high = 0
        self.failed = 0
        self.avg = 0.0
        self.std = 0.0

class MemAln:
    def __init__(self):
        self.pos = 0
        self.rid = 0
        self.flag = 0
        self.is_rev = 0
        self.is_alt = 0
        self.mapq = 0
        self.NM = 0
        self.n_cigar = 0
        self.cigar = []
        self.XA = None
        self.score = 0
        self.sub = 0
        self.alt_sc = 0


# Thiết lập các tham số
class MemOpt:
    def __init__(self):
        self.flag = 0
        self.a = 1
        self.b = 4
        self.o_del = self.o_ins = 6
        self.e_del = self.e_ins = 1
        self.w = 100
        self.T = 30
        self.zdrop = 100
        self.pen_unpaired = 17
        self.pen_clip5 = self.pen_clip3 = 5

        self.max_mem_intv = 20

        self.min_seed_len = 19
        self.split_width = 10
        self.max_occ = 500
        self.max_chain_gap = 10000
        self.max_ins = 10000
        self.mask_level = 0.50
        self.drop_ratio = 0.50
        self.XA_drop_ratio = 0.80
        self.split_factor = 1.5
        self.chunk_size = 10000000
        self.n_threads = 1
        self.max_XA_hits = 5
        self.max_XA_hits_alt = 200
        self.max_matesw = 50
        self.mask_level_redun = 0.95
        self.min_chain_weight = 0
        self.max_chain_extend = 1 << 30
        
        self.mapQ_coef_len = 50
        self.mapQ_coef_fac = math.log(self.mapQ_coef_len)
        
        self.mat = [[0] * 5 for _ in range(5)]  # Giả lập ma trận điểm số
        self.bwa_fill_scmat()

    def bwa_fill_scmat(self):
        # Đây là giả lập hàm fill_scmat từ C
        for i in range(5):
            for j in range(5):
                self.mat[i][j] = self.a if i == j else -self.b


def mem_opt_init():
    return MemOpt()

#===============================================================================
# Thu thập và xử lý các khoảng suffix array (SA)
class BwtIntv:
    #Lớp đại diện cho một khoảng bwtintv.
    def __init__(self, info, x):
        self.info = info  # Thông tin về khoảng
        self.x = x  # Danh sách chứa các giá trị liên quan

class SmemAux:
    # Lớp hỗ trợ xử lý các khoảng SA.
    def __init__(self):
        self.mem = []  # Danh sách lưu trữ các khoảng SMEM
        self.mem1 = []  # Danh sách tạm thời cho SMEM
        self.tmpv = [[], []]  # Hai danh sách tạm thời dùng để xử lý dữ liệu

    @staticmethod
    def smem_aux_init():
        # Khởi tạo đối tượng SmemAux mới.
        return SmemAux()

    def smem_aux_destroy(self):
        # Giải phóng bộ nhớ bằng cách làm sạch các danh sách.
        self.tmpv[0].clear()
        self.tmpv[1].clear()
        self.mem.clear()
        self.mem1.clear()

    def mem_collect_intv(self, opt, bwt, length, seq):
        #Thu thập các khoảng SA dựa trên thông tin trình tự và tham số đầu vào.
        x = 0
        start_width = 1  # Độ rộng bắt đầu
        split_len = int(opt.min_seed_len * opt.split_factor + 0.499)  # Độ dài tối thiểu để chia
        self.mem.clear()
        
        # Lượt đầu tiên: tìm tất cả SMEM
        while x < length:
            if seq[x] < 4:
                x = self.bwt_smem1(bwt, length, seq, x, start_width, self.mem1, self.tmpv)
                for p in self.mem1:
                    slen = (p.info & 0xFFFFFFFF) - (p.info >> 32)  # Tính độ dài seed
                    if slen >= opt.min_seed_len:
                        self.mem.append(p)
            else:
                x += 1

        # Lượt thứ hai: tìm MEM trong một SMEM dài
        old_n = len(self.mem)
        for k in range(old_n):
            p = self.mem[k]
            start = p.info >> 32
            end = p.info & 0xFFFFFFFF
            if end - start < split_len or p.x[2] > opt.split_width:
                continue
            self.bwt_smem1(bwt, length, seq, (start + end) >> 1, p.x[2] + 1, self.mem1, self.tmpv)
            for p in self.mem1:
                if (p.info & 0xFFFFFFFF) - (p.info >> 32) >= opt.min_seed_len:
                    self.mem.append(p)
        
        # Lượt thứ ba: Chiến lược giống LAST
        if opt.max_mem_intv > 0:
            x = 0
            while x < length:
                if seq[x] < 4:
                    m = BwtIntv(0, [0, 0, 0])
                    x = self.bwt_seed_strategy1(bwt, length, seq, x, opt.min_seed_len, opt.max_mem_intv, m)
                    if m.x[2] > 0:
                        self.mem.append(m)
                else:
                    x += 1
        
        # Sắp xếp kết quả theo thông tin khoảng
        self.mem.sort(key=lambda p: p.info)

    def bwt_smem1(self, bwt, length, seq, x, start_width, mem1, tmpv):
        """
        Placeholder: Hàm tìm kiếm SMEM dựa trên thuật toán BWT.
        """
        # Thực hiện tìm kiếm SMEM và cập nhật mem1 và tmpv
        # Trả về chỉ số tiếp theo (giả định)
        return x + 1
    
    def bwt_seed_strategy1(self, bwt, length, seq, x, min_seed_len, max_mem_intv, m):
        """
        Placeholder: Hàm chiến lược tìm seed BWT.
        """
        # Thực hiện chiến lược tìm seed và cập nhật m
        # Trả về chỉ số tiếp theo (giả định)
        return x + 1

#===============================================================================
# Chaining

#===============================================================================
# Filtering chains

#===============================================================================
# De-overlap single-end hits

#===============================================================================
# Test if a seed is good enough

#===============================================================================
def cal_max_gap(opt, qlen):
    """
    Tính toán khoảng cách tối đa (max_gap) giữa hai đoạn căn chỉnh trong chuỗi truy vấn (query sequence).
    
    Tham Số:
        - opt: dict	Tham số căn chỉnh.
        - qlen: int	Độ dài của chuỗi truy vấn.

    Trả về:
        - int	Khoảng cách tối đa.
    """
    # Tính khoảng cách tối đa khi có thao tác xóa (deletion)
    l_del = int((qlen * opt['a'] - opt['o_del']) / opt['e_del'] + 1)

    # Tính khoảng cách tối đa khi có thao tác chèn (insertion)
    l_ins = int((qlen * opt['a'] - opt['o_ins']) / opt['e_ins'] + 1)\

    l = max(l_del, l_ins)
    l = max(l, 1) # Đảm bảo giá trị tối thiểu của khoảng cách là 1
    return min(l, opt['w'] << 1) # Giới hạn khoảng cách tối đa

MAX_BAND_TRY = 2

def mem_chain2aln(opt, bns, pac, l_query, query, c, av):
    rmax = [0, 0]
    max_off = [0, 0]
    aw = [0, 0]
    max_len = 0

    if c['n'] == 0:
        return

    # get the max possible span
    rmax[0] = float('inf')
    rmax[1] = 0
    for seed in c['seeds']:
        b = seed['rbeg'] - (seed['qbeg'] + cal_max_gap(opt, seed['qbeg']))
        e = seed['rbeg'] + seed['len'] + ((l_query - seed['qbeg'] - seed['len']) + cal_max_gap(opt, l_query - seed['qbeg'] - seed['len']))
        rmax[0] = min(rmax[0], b)
        rmax[1] = max(rmax[1], e)
        max_len = max(max_len, seed['len'])

    rmax[0] = max(rmax[0], 0)
    rmax[1] = min(rmax[1], len(pac) * 2)

    if rmax[0] < len(pac) and len(pac) < rmax[1]: # crossing the forward-reverse boundary; then choose one side
        if c['seeds'][0]['rbeg'] < len(pac): # this works because all seeds are guaranteed to be on the same strand
            rmax[1] = len(pac)
        else:
            rmax[0] = len(pac)

    # retrieve the reference sequence
    rmax_beg, rmax_end = rmax  # Giữ giá trị ban đầu của rmax[0] và rmax[1]
    rseq, rmax_beg, rmax_end, rid = bns_fetch_seq(bns, pac, rmax_beg, c.seeds[0].rbeg, rmax_end)

    assert c.rid == rid

    srt = [(c.seeds[i].score << 32) | i for i in range(c.n)]
    # "ks_introsort_64(c->n, srt);" khúc này chưa biết hàm ks_introsort_64 nằm ở đâu
    srt.sort()

    for k in range(c.n - 1, -1, -1):  # Duyệt từ c.n - 1 về 0
        s = c.seeds[srt[k] & 0xFFFFFFFF] # (uint32_t) để lấy phần 32-bit cuối cùng

        for i in range(av.n):  # Kiểm tra xem có mở rộng trước đó không
            p = av.a[i]

            # Kiểm tra nếu seed không nằm hoàn toàn trong p
            if s.rbeg < p.rb or s.rbeg + s.len > p.re or s.qbeg < p.qb or s.qbeg + s.len > p.qe:
                continue

            # Nếu seed có thể tạo ra căn chỉnh tốt hơn
            if s.len - p.seedlen0 > 0.1 * l_query:  
                continue

            # Tính toán khoảng cách trước seed trên query/reference
            qd = s.qbeg - p.qb
            rd = s.rbeg - p.rb
            max_gap = cal_max_gap(opt, min(qd, rd))  # Giới hạn khoảng trống tối đa
            w = min(max_gap, p.w)  # Giới hạn bởi băng thông

            if abs(qd - rd) < w:
                break  # Seed nằm gần với một lần hit trước đó

            # Tương tự nhưng kiểm tra vùng phía sau
            qd = p.qe - (s.qbeg + s.len)
            rd = p.re - (s.rbeg + s.len)
            max_gap = cal_max_gap(opt, min(qd, rd))
            w = min(max_gap, p.w)

            if abs(qd - rd) < w:
                break  # Seed nằm gần với một lần hit trước đó

        if i < len(av):  # The seed is (almost) contained in an existing alignment; further testing is needed
            if bwa_verbose >= 4:
                    print(f"** Seed({k}) [{s.len};{s.qbeg},{s.rbeg}] is almost contained in an existing alignment [{av[i].qb},{av[i].qe}) <=> [{av[i].rb},{av[i].re})")

            for i in range(k + 1, c.n):  # Check overlapping seeds in the same chain
                if srt[i] == 0:
                    continue # Skip if already marked
            
                t = c.seeds[srt[i] & 0xFFFFFFFF]  # Equivalent to (uint32_t)srt[i] in C

                if t.len < s.len * 0.95:
                    continue # Only check overlapping if t is long enough

                if (s.qbeg <= t.qbeg and (s.qbeg + s.len - t.qbeg) >= (s.len >> 2) and 
                    (t.qbeg - s.qbeg) != (t.rbeg - s.rbeg)):
                    break

                if (t.qbeg <= s.qbeg and (t.qbeg + t.len - s.qbeg) >= (s.len >> 2) and 
                    (s.qbeg - t.qbeg) != (s.rbeg - t.rbeg)):
                    break
            
            if i == c.n:  # No overlapping seeds; then skip extension
                srt[k] = 0  # Mark that seed extension has not been performed
                continue

            if bwa_verbose >= 4:
                print(f"** Seed({k}) might lead to a different alignment even though it is contained. Extension will be performed.")

        a = av.kv_pushp()  # Thêm phần tử mới vào danh sách và lấy tham chiếu đến nó
        a.w = aw[0] = aw[1] = opt.w
        a.score = a.truesc = -1
        a.rid = c.rid

        if bwa_verbose >= 4:
            print(f"** ---> Extending from seed({k}) [{s.len};{s.qbeg},{s.rbeg}] @ {bns.anns[c.rid].name} <---", 
                file=sys.stderr)
        if s.qbeg:  # left extension
            qs = [query[s.qbeg - 1 - i] for i in range(s.qbeg)]
            tmp = s.rbeg - rmax[0]
            rs = [rseq[tmp - 1 - i] for i in range(tmp)]  # Tạo mảng rs (đảo ngược chuỗi tham chiếu)

            for i in range(MAX_BAND_TRY):
                prev = a.score
                aw[0] = opt.w << i  # Nhân đôi bandwidth sau mỗi vòng lặp

                if bwa_verbose >= 4:
                    ref_str = "".join("ACGTN"[r] for r in rs)
                    query_str = "".join("ACGTN"[q] for q in qs)
                    print(f"*** Left ref:   {ref_str}")
                    print(f"*** Left query: {query_str}")
                
                # Gọi ksw_extend2 để mở rộng
                qle, tle, gtle, gscore, max_off = [0], [0], [0], [0], [0]  # Khởi tạo danh sách để truyền vào hàm
                a.score = ksw_extend2(
                    s.qbeg, qs, tmp, rs, 5, opt.mat, opt.o_del, opt.e_del,
                    opt.o_ins, opt.e_ins, aw[0], opt.pen_clip5, opt.zdrop,
                    s.len * opt.a, qle, tle, gtle, gscore, max_off
                )
                if bwa_verbose >= 4:
                    print(f"*** Left extension: prev_score={prev}; score={a.score}; bandwidth={aw[0]}; max_off_diagonal_dist={max_off[0]}")

                # Nếu điểm số không cải thiện hoặc max_off[0] nhỏ hơn ngưỡng, thoát vòng lặp
                if a.score == prev or max_off[0] < (aw[0] >> 1) + (aw[0] >> 2):
                    break

            # Kiểm tra xem có nên mở rộng về cuối query không
            if gscore[0] <= 0 or gscore[0] <= a.score - opt.pen_clip5:
                # Local extension
                a.qb = s.qbeg - qle[0]
                a.rb = s.rbeg - tle[0]
                a.truesc = a.score
            else:
                # To-end extension
                a.qb = 0
                a.rb = s.rbeg - gtle[0]
                a.truesc = gscore[0]
            # Giải phóng bộ nhớ nếu cần
            del qs, rs
        else:
            a.score = a.truesc = s.len * opt.a
            a.qb = 0
            a.rb = s.rbeg

        if s.qbeg + s.len != l_query:  # right extension
            qle = tle = qe = re = gtle = gscore = 0
            sc0 = a.score
            qe = s.qbeg + s.len
            re = s.rbeg + s.len - rmax[0]
            assert re >= 0

            for i in range(MAX_BAND_TRY):
                prev = a.score
                aw[1] = opt.w << i  # dịch bit sang trái i lần  

                if bwa_verbose >= 4:
                    print("*** Right ref:   ", "".join("ACGTN"[rseq[re + j]] for j in range(rmax[1] - rmax[0] - re)))
                    print("*** Right query: ", "".join("ACGTN"[query[qe + j]] for j in range(l_query - qe)))

                a.score = ksw_extend2(
                    l_query - qe, query[qe:], rmax[1] - rmax[0] - re, rseq[re:], 5, opt.mat,
                    opt.o_del, opt.e_del, opt.o_ins, opt.e_ins, aw[1], opt.pen_clip3, opt.zdrop,
                    sc0, qle, tle, gtle, gscore, max_off[1]
                )

                if bwa_verbose >= 4:
                    print(f"*** Right extension: prev_score={prev}; score={a.score}; bandwidth={aw[1]}; max_off_diagonal_dist={max_off[1]}")
                if a.score == prev or max_off[1] < (aw[1] >> 1) + (aw[1] >> 2):
                    break

            if gscore <= 0 or gscore <= a.score - opt.pen_clip3:  # local extension
                a.qe = qe + qle
                a.re = rmax[0] + re + tle
                a.truesc += a.score - sc0
            else:  # to-end extension
                a.qe = l_query
                a.re = rmax[0] + re + gtle
                a.truesc += gscore - sc0
        else:
            a.qe = l_query
            a.re = s.rbeg + s.len
        
        if bwa_verbose >= 4:
            print(f"*** Added alignment region: [{a.qb},{a.qe}) <=> [{a.rb},{a.re}); score={a.score}; {{left,right}}_bandwidth={{ {aw[0]},{aw[1]} }}")
        
        # Compute seed coverage
        a.seedcov = sum(
            t.len for t in c.seeds
            if a.qb <= t.qbeg < a.qe and a.rb <= t.rbeg < a.re and t.qbeg + t.len <= a.qe and t.rbeg + t.len <= a.re
        )  # Seed fully contained

        # Set other attributes
        a.w = max(aw[0], aw[1])
        a.seedlen0 = s.len
        a.frac_rep = c.frac_rep
    del srt, rseq


            
        

#===============================================================================
# Basic hit->SAM conversion

#===============================================================================
# Integrated interface



def mem_process_seqs(opt, bwt, bns, pac, n_processed, n, seqs, pes0):
    pass

# Sử dụng hàm
if __name__ == "__main__":
    opt = mem_opt_init()