from bwt_h import *
import struct
from utils import *
import os

def __occ_aux4(bwt, b):
    """Tính tổng số lần xuất hiện của các ký tự trong một đoạn của chuỗi BWT"""
    return (bwt.cnt_table[b & 0xff] +
            bwt.cnt_table[(b >> 8) & 0xff] +
            bwt.cnt_table[(b >> 16) & 0xff] +
            bwt.cnt_table[(b >> 24)])

def bwt_occ4(bwt, k, cnt):
    """Đếm số lần xuất hiện của từng ký tự (A, C, G, T) trong chuỗi BWT từ vị trí 0 đến một vị trí nhất định k"""
    if k == -1:
        cnt[:] = [0] * 4
        return
    k -= int(k >= bwt.primary)
    p = bwt_occ_intv(bwt, k)
    cnt[:] = p[:4]
    p = p[4:]
    end = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4))
    x = sum(__occ_aux4(bwt, val) for val in p[:end])
    tmp = p[end] & ~((1 << ((~k & 15) << 1)) - 1)
    x += __occ_aux4(bwt, tmp) - (~k & 15)
    cnt[0] += x & 0xff
    cnt[1] += (x >> 8) & 0xff
    cnt[2] += (x >> 16) & 0xff
    cnt[3] += x >> 24

def bwt_2occ4(bwt, k, l, cntk, cntl):
    _k = k - int(k >= bwt.primary)
    _l = l - int(l >= bwt.primary)
    if (_l >> OCC_INTV_SHIFT) != (_k >> OCC_INTV_SHIFT) or k == -1 or l == -1:
        bwt_occ4(bwt, k, cntk)
        bwt_occ4(bwt, l, cntl)
    else:
        k -= int(k >= bwt.primary)
        l -= int(l >= bwt.primary)
        p = bwt_occ_intv(bwt, k)
        cntk[:] = p[:4]
        p = p[4:]
        endk = p + ((k >> 4) - ((k & ~OCC_INTV_MASK) >> 4))
        endl = p + ((l >> 4) - ((l & ~OCC_INTV_MASK) >> 4))
        x = sum(__occ_aux4(bwt, val) for val in p[:endk])
        y = x
        tmp = p[endk] & ~((1 << ((~k & 15) << 1)) - 1)
        x += __occ_aux4(bwt, tmp) - (~k & 15)
        for val in p[endk:endl]:
            y += __occ_aux4(bwt, val)
        tmp = p[endl] & ~((1 << ((~l & 15) << 1)) - 1)
        y += __occ_aux4(bwt, tmp) - (~l & 15)
        cntl[:] = cntk[:]
        cntk[0] += x & 0xff
        cntk[1] += (x >> 8) & 0xff
        cntk[2] += (x >> 16) & 0xff
        cntk[3] += x >> 24
        cntl[0] += y & 0xff
        cntl[1] += (y >> 8) & 0xff
        cntl[2] += (y >> 16) & 0xff
        cntl[3] += y >> 24

def bwt_extend(bwt, ik, ok, is_back):
    """Mở rộng một khoảng (interval) trong BWT bằng cách tính toán sự xuất hiện của từng ký tự trong khoảng đó."""
    tk = [0] * 4
    tl = [0] * 4

    # Tính số lần xuất hiện của các ký tự từ ik.x[not is_back] - 1 đến ik.x[not is_back] - 1 + ik.x[2]
    bwt_2occ4(bwt, ik.x[not is_back] - 1, ik.x[not is_back] - 1 + ik.x[2], tk, tl)
    for i in range(4):
        ok[i].x[not is_back] = bwt.L2[i] + 1 + tk[i] # Bắt đầu khoảng mới
        ok[i].x[2] = tl[i] - tk[i] # Độ dài khoảng mới

    # Điều chỉnh vị trí của khoảng mở rộng
    ok[3].x[is_back] = ik.x[is_back] + (ik.x[not is_back] <= bwt.primary and ik.x[not is_back] + ik.x[2] - 1 >= bwt.primary)
    ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2]
    ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2]
    ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2]

def bwt_reverse_intvs(p):
    """Đảo ngược thứ tự của một danh sách các khoảng (intervals) trong BWT"""
    if p.n > 1:
        for j in range(p.n // 2):
            p.a[j], p.a[p.n - 1 - j] = p.a[p.n - 1 - j], p.a[j]

def bwt_smem1a(bwt, length, q, x, min_intv, max_intv, mem, tmpvec):
    """
    Tham số:
    - bwt       Cấu trúc dữ liệu BWT chứa thông tin đã nén của chuỗi gốc
    - length	Độ dài của chuỗi truy vấn q
    - q	        Chuỗi truy vấn (dưới dạng danh sách số nguyên, với A=0, C=1, G=2, T=3)
    - x	        Vị trí bắt đầu tìm kiếm trong q
    - min_intv	Kích thước khoảng nhỏ nhất được chấp nhận (thường là 1)
    - max_intv	Kích thước khoảng lớn nhất (giới hạn số lần xuất hiện của SMEM)
    - mem	    Danh sách lưu các SMEM tìm được
    - tmpvec	Bộ đệm tạm thời để tối ưu hiệu suất
    """
    mem.n = 0
    if q[x] > 3:
        return x + 1
    if min_intv < 1:
        min_intv = 1
    
    a = [bwtintv_v(), bwtintv_v()]
    prev = tmpvec[0] if tmpvec and tmpvec[0] else a[0]
    curr = tmpvec[1] if tmpvec and tmpvec[1] else a[1]
    
    ik = bwtintv_t([0, 0, 0], 0)
    bwt_set_intv(bwt, q[x], ik)
    ik.info = x + 1
    
    i = x + 1
    curr.n = 0
    while i < length:
        if ik.x[2] < max_intv:
            curr.a.append(ik)
            curr.n += 1
            break
        elif q[i] < 4:
            c = 3 - q[i]
            ok = [bwtintv_t([0, 0, 0], 0) for _ in range(4)]
            bwt_extend(bwt, ik, ok, 0)
            if ok[c].x[2] != ik.x[2]:
                curr.a.append(ik)
                curr.n += 1
                if ok[c].x[2] < min_intv:
                    break
            ik = ok[c]
            ik.info = i + 1
        else:
            curr.a.append(ik)
            curr.n += 1
            break
        i += 1
    
    if i == length:
        curr.a.append(ik)
        curr.n += 1
    
    bwt_reverse_intvs(curr)
    ret = curr.a[0].info
    prev, curr = curr, prev
    
    for i in range(x - 1, -2, -1):
        c = -1 if i < 0 or q[i] >= 4 else q[i]
        curr.n = 0
        for j in range(prev.n):
            p = prev.a[j]
            if c >= 0 and ik.x[2] >= max_intv:
                ok = [bwtintv_t([0, 0, 0], 0) for _ in range(4)]
                bwt_extend(bwt, p, ok, 1)
            if c < 0 or ik.x[2] < max_intv or ok[c].x[2] < min_intv:
                if curr.n == 0:
                    if mem.n == 0 or i + 1 < (mem.a[mem.n - 1].info >> 32):
                        ik = p
                        ik.info |= (i + 1) << 32
                        mem.a.append(ik)
                        mem.n += 1
            elif curr.n == 0 or ok[c].x[2] != curr.a[curr.n - 1].x[2]:
                ok[c].info = p.info
                curr.a.append(ok[c])
                curr.n += 1
        
        if curr.n == 0:
            break
        prev, curr = curr, prev
    
    bwt_reverse_intvs(mem)
    return ret

def bwt_gen_cnt_table(bwt):
    """Sinh bảng đếm cho BWT."""
    bwt.cnt_table = [0] * 256

    for i in range(256):
        x = 0
        for j in range(4):
            x |= (((i & 3) == j) +
                  ((i >> 2 & 3) == j) +
                  ((i >> 4 & 3) == j) +
                  ((i >> 6) == j)) << (j << 3)
        bwt.cnt_table[i] = x

def bwt_restore_sa(filename, bwt):
    """Khôi phục mảng suffix array (SA) từ tệp chỉ mục .sa"""
    with open(filename, "rb") as fp:
        primary = struct.unpack("Q", err_fread_noeof(fp, 8, 1))[0]
        assert primary == bwt.primary, "SA-BWT inconsistency: primary is not the same."

        # Bỏ qua 4 giá trị bwtint_t (32 byte)
        err_fread_noeof(fp, 8, 4)

        bwt.sa_intv = struct.unpack("Q", err_fread_noeof(fp, 8, 1))[0]
        primary = struct.unpack("Q", err_fread_noeof(fp, 8, 1))[0]
        assert primary == bwt.seq_len, "SA-BWT inconsistency: seq_len is not the same."

        # Tính số lượng phần tử trong SA
        bwt.n_sa = (bwt.seq_len + bwt.sa_intv) // bwt.sa_intv
        bwt.sa = [-1] + list(struct.unpack(f"{bwt.n_sa - 1}Q", err_fread_noeof(fp, 8, bwt.n_sa - 1)))

def bwt_restore_bwt(fn):
    """Khôi phục BWT từ tệp .bwt"""
    bwt = bwt_t()

    # Mở file ở chế độ nhị phân
    with open(fn, "rb") as fp:
        # Tìm kích thước tệp
        fp.seek(0, os.SEEK_END)
        file_size = fp.tell()
        bwt.bwt_size = (file_size - 5 * 8) >> 2  # sizeof(bwtint_t) = 8 bytes trong C

        # Cấp phát bộ nhớ cho bwt
        bwt.bwt = [0] * bwt.bwt_size

        # Đọc dữ liệu
        fp.seek(0, os.SEEK_SET)
        bwt.primary = struct.unpack("<Q", fp.read(8))[0]  # bwtint_t (8 bytes)
        bwt.L2[1:5] = struct.unpack("<4Q", fp.read(32))  # 4 * bwtint_t (32 bytes)
        bwt.bwt = list(struct.unpack(f"<{bwt.bwt_size}I", fp.read(bwt.bwt_size * 4)))  # uint32_t (4 bytes mỗi phần tử)
        bwt.seq_len = bwt.L2[4]

    # Sinh bảng đếm
    bwt_gen_cnt_table(bwt)

    return bwt