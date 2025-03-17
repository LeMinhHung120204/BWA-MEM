KSW_XBYTE = 0x10000
KSW_XSTOP = 0x20000
KSW_XSUBO = 0x40000
KSW_XSTART = 0x80000

class KSWR:
    """Lưu trữ kết quả của thuật toán mở rộng Smith-Waterman."""
    def __init__(self, score=-1, te=-1, qe=-1, score2=-1, te2=-1, tb=-1, qb=-1):
        self.score = score      # Điểm số tối đa tìm thấy trong quá trình mở rộng
        self.te = te            # Chỉ số trên chuỗi target (trình tự tham chiếu) nơi xảy ra điểm tối đa
        self.qe = qe            # Chỉ số trên chuỗi query (trình tự truy vấn) nơi xảy ra điểm tối đa
        self.score2 = score2    # Điểm số tối đa thứ hai
        self.te2 = te2          # Chỉ số target tương ứng với score2
        self.tb = tb            # Vị trí bắt đầu của đoạn căn chỉnh trên target
        self.qb = qb            # Vị trí bắt đầu của đoạn căn chỉnh trên query

class KSWQ:
    """Lưu trữ thông tin về query để tối ưu hóa quá trình tính toán động (DP)."""
    def __init__(self, qlen, slen, shift, mdiff, max, size, qp, H0, H1, E, Hmax):
        self.qlen = qlen    # Độ dài chuỗi query
        self.slen = slen    # Độ dài chuỗi target
        self.shift = shift  # Dịch bit trong ma trận điểm số (không thường dùng)
        self.mdiff = mdiff  # Sự khác biệt tối đa có thể có trong điểm số
        self.max = max      # Giá trị điểm số tối đa có thể có
        self.size = size    # Kích thước lưu trữ điểm số (1 byte hoặc 2 byte)
        self.qp = qp        # Query profile (một danh sách lưu điểm số so với ma trận)
        self.H0 = H0        # Mảng điểm số của hàng trước trong DP
        self.H1 = H1        # Mảng điểm số của hàng hiện tại trong DP
        self.E = E          # Mảng lưu điểm số cho khoảng cách mở rộng
        self.Hmax = Hmax    # Điểm số lớn nhất trong DP


"""
size    Number of bytes used to store a score; valid valures are 1 or 2
qlen    Length of the query sequence
query   Query sequence
m       Size of the alphabet
mat     Scoring matrix in a one-dimension array
return  Query data structure
"""

# SW extension

class EHT:
    """Giúp lưu trữ trạng thái trong quá trình tính toán DP."""
    def __init__(self, h=0, e=0):
        self.h = h  # Điểm số tối ưu nhất tại vị trí (i, j) trong DP
        self.e = e  # Giá trị mở rộng gap (gap extension)


def ksw_extend2(qlen, query, tlen, target, m, mat, o_del, e_del, o_ins, e_ins, w, end_bonus, zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off):
    """
    Căn chỉnh một chuỗi query (query) với một chuỗi target (target) bằng cách sử dụng thuật toán Smith-Waterman có ràng buộc dải (banded Smith-Waterman).

    Tham số:
        - qlen	    int	        Độ dài chuỗi query.
        - query	    list[int]	Danh sách chứa mã ASCII của các nucleotide trong chuỗi query.
        - tlen	    int	        Độ dài chuỗi target.
        - target    list[int]	Danh sách chứa mã ASCII của các nucleotide trong chuỗi target.
        - m	        int	        Kích thước bảng điểm (số lượng nucleotide).
        - mat	    list[int]	Mảng chứa điểm số của các cặp nucleotide.
        - o_del	    int	        Phạt cho việc xóa một nucleotide (gap open).
        - e_del	    int	        Phạt cho việc xóa mở rộng một nucleotide (gap extension).
        - o_ins	    int	        Phạt cho việc chèn một nucleotide (gap open).
        - e_ins	    int	        Phạt cho việc chèn mở rộng một nucleotide (gap extension).
        - w	        int	        Chiều rộng của dải (band).
        - end_bonus	int	        Phần thưởng cho việc kết thúc chuỗi.
        - zdrop	    int	        Phần thưởng cho việc dừng quá trình mở rộng.
        - h0	    int	        Điểm số ban đầu.
        - _qle	    list[int]	Chỉ số cuối cùng của chuỗi query được căn chỉnh.
        - _tle	    list[int]	Chỉ số cuối cùng của chuỗi target được căn chỉnh.
        - _gtle	    list[int]	Chỉ số cuối cùng của chuỗi target được căn chỉnh (có điểm số cao nhất).
        - _gscore	list[int]	Điểm số cao nhất tìm thấy trong quá trình mở rộng.
        - _max_off	list[int]	Khoảng cách lớn nhất giữa chỉ số của chuỗi target và query.

    Trả về:
        - int	Điểm số tối đa tìm thấy trong quá trình m
    """
    eh = [EHT() for _ in range(qlen + 1)] # score array
    qp = [0] * (qlen * m) # query profile
    # Tính toán giới hạn tối đa của phép chèn (max_ins) và xóa (max_del) dựa trên ma trận điểm.
    oe_del = o_del + e_del
    oe_ins = o_ins + e_ins

    # adjust $w if it is too large
    # Xác định số lượng chèn (insertions) tối đa có thể xảy
    max_ins = int((qlen * max(mat) + end_bonus - o_ins) / e_ins + 1)
    max_ins = max(max_ins, 1)
    w = min(w, max_ins)

    # Xác định số lượng xóa (deletions) tối đa có thể xảy
    max_del = int((qlen * max(mat) + end_bonus - o_del) / e_del + 1)
    max_del = max(max_del, 1)
    w = min(w, max_del)

    # generate the query profile
    for k in range(m):
        p = mat[k * m:(k + 1) * m]
        for j in range(qlen):
            qp[k * qlen + j] = p[query[j]]

    # for i in range(m):
    #     print(qp[i * qlen: (i + 1) * qlen])
    # print("\n")

    # fill the first row
    eh[0].h = h0 # Khởi tạo phần tử đầu tiên của hàng với điểm số ban đầu h0
    eh[1].h = max(h0 - oe_ins, 0)
    for j in range(2, qlen + 1):
        eh[j].h = max(eh[j - 1].h - e_ins, 0)

    # DP loop
    # Khởi tạo các biến theo dõi điểm số và vị trí tốt nhất
    max_score = h0      # Điểm số cao nhất, ban đầu bằng điểm khởi tạo h0
    max_i = max_j = -1  # Vị trí (i,j) của điểm số cao nhất trong ma trận
    max_ie = -1         # Vị trí i của điểm cuối cùng trên target khi đạt đến cuối query
    gscore = -1         # Điểm số cao nhất khi đạt đến cuối query
    max_off = 0         # Độ lệch tối đa giữa vị trí trên query và target
    beg = 0             # Vị trí bắt đầu của dải trên query
    end = qlen          # Vị trí kết thúc của dải trên query (ban đầu là toàn bộ query)

    # thực hiện quét qua chuỗi target để tính toán ma trận động Smith-Waterman với ràng buộc dải (banding)
    for i in range(tlen):
        f = 0   # Giá trị gap mở rộng theo chiều ngang (F)
        h1 = 0 if beg == 0 else -1 # Giá trị H(i, j-1), nếu ở cột đầu tiên thì đặt về 0, ngược lại đặt -1
        q = qp[target[i] * qlen:(target[i] + 1) * qlen] # Truy xuất query profile của ký tự target[i]
        # chứa điểm số so sánh giữa ký tự target[i] với tất cả các ký tự trong query

        # apply the band and the constraint (if provided)
        # Điều chỉnh biên (beg, end) để chỉ xét trong phạm vi dải (band)
        if beg < i - w:
            beg = i - w
        if end > i + w + 1:
            end = i + w + 1
        if end > qlen:
            end = qlen

        # compute the first column
        # Tính giá trị H(i, 0) cho cột đầu tiên
        if beg == 0:
            h1 = h0 - (o_del + e_del * (i + 1)) # H(i,0) = H(0,0) - gap mở rộng
            if h1 < 0:
                h1 = 0
        else:
            h1 = 0

        m = 0  # Điểm số cao nhất trong hàng hiện tại
        mj = beg  # Vị trí của điểm số cao nhất
        # Vòng lặp chính tính toán ma trận DP
        for j in range(beg, end):
            # At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
            # Similar to SSE2-SW, cells are computed in the following order:
            # H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
            # E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
            # F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
            p = eh[j]
            M, e = p.h, p.e  # Lấy H(i-1,j-1) và E(i-1,j)
            p.h = h1  # Cập nhật H(i,j-1) cho dòng tiếp theo

            # Cập nhật điểm số căn chỉnh H(i, j)
            M = (M + q[j]) if M else 0  # Nếu M != 0 thì cộng với q[j], ngược lại đặt bằng 0
            h = max(M, e)  # So sánh với e
            h = max(h, f)  # So sánh với f
            h1 = h  # Lưu H(i,j) vào h1 cho cột tiếp theo

            mj = mj if m > h else j  # Cập nhật vị trí có điểm số cao nhất
            m = max(m, h)  # Lưu giá trị điểm số lớn nhất

            t = max(M - oe_del, 0)  # t = M - oe_del, nếu <0 thì đặt = 0
            e = max(e - e_del, t)  # Tính E(i+1,j)
            p.e = e  # Lưu E(i+1,j) vào p.e

            t = max(M - oe_ins, 0)  # t = M - oe_ins, nếu <0 thì đặt = 0
            f = max(f - e_ins, t)  # Tính F(i,j+1)
            
        # Lưu giá trị H(i, end) và reset E(i, end).
        eh[end].h = h1
        eh[end].e = 0

        # Cập nhật gscore nếu đạt đến cuối query
        if end == qlen:
            max_ie = i if gscore < h1 else max_ie   # Cập nhật vị trí i của điểm số cao nhất
            gscore = max(gscore, h1)                # Cập nhật điểm số cao nhất

        # Kiểm tra điểm số lớn nhất (max_score)
        if max_score == 0:
            break

        # Cập nhật điểm số cao nhất và vị trí tương ứng
        if m > max_score:
            # m là điểm số cao nhất trong hàng hiện tại
            # i là vị trí hiện tại trên target
            # mj là vị trí trên query có điểm số cao nhất trong hàng
            max_score, max_i, max_j = m, i, mj

            # Cập nhật độ lệch tối đa giữa vị trí trên query và target
            max_off = max(max_off, abs(mj - i))
        elif zdrop > 0: # Điều kiện dừng sớm (Z-drop)
            # So sánh độ chênh lệch giữa vị trí hiện tại và vị trí có điểm số cao nhất
            if i - max_i > mj - max_j:
                # Trường hợp deletion
                # Tính penalty dựa trên khoảng cách và gap extension
                penalty = ((i - max_i) - (mj - max_j)) * e_del
                if max_score - m - penalty > zdrop:
                    break
            else:
                # Trường hợp insertion (chèn)
                # Tính penalty dựa trên khoảng cách và gap extension
                penalty = ((mj - max_j) - (i - max_i)) * e_ins
                if max_score - m - penalty > zdrop:
                    break

        # update beg and end for the next round
        # Cập nhật vị trí bắt đầu (beg)
        """
        - Tìm vị trí đầu tiên mà cả h (điểm số) và e (gap extension) đều bằng 0
        - Điểm số 0 nghĩa là không có căn chỉnh có ý nghĩa tại vị trí đó
        - Vị trí này sẽ là điểm bắt đầu mới cho dải trong lần lặp tiếp theo
        """
        for j in range(beg, end):
            if eh[j].h == 0 and eh[j].e == 0: # Tìm vị trí đầu tiên có điểm số = 0
                beg = j
                break
        
        # Cập nhật vị trí kết thúc (end)
        """
        - Tìm vị trí cuối cùng có điểm số = 0
        - Đặt end mới = vị trí đó + 2 (để có khoảng đệm)
        - Nếu end mới vượt quá qlen thì đặt end = qlen
        """
        for j in range(end - 1, beg - 1, -1):   # Duyệt ngược từ cuối lên đầu
            if eh[j].h == 0 and eh[j].e == 0:   # Tìm vị trí cuối cùng có điểm số = 0
                end = j + 2 if j + 2 < qlen else qlen # Đặt end = j+2 hoặc qlen
                break

    # Giải phóng bộ nhớ
    del eh, qp
    if _qle is not None:
        _qle[0] = max_j + 1     # Vị trí kết thúc trên query (0-based -> 1-based)
    if _tle is not None:
        _tle[0] = max_i + 1     # Vị trí kết thúc trên target (0-based -> 1-based)
    if _gtle is not None:
        _gtle[0] = max_ie + 1   # Vị trí kết thúc trên target khi đạt điểm cao nhất tại cuối query
    if _gscore is not None:
        _gscore[0] = gscore     # Điểm số cao nhất khi đạt đến cuối query
    if _max_off is not None:
        _max_off[0] = max_off   # Độ lệch tối đa giữa vị trí trên query và target

    return max_score

def ksw_extend(qlen, query, tlen, target, m, mat, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off):
    return ksw_extend2(qlen, query, tlen, target, m, mat, gapo, gape, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off)

# Global alignment

# Main function (not compiled by default)

if __name__ == "__main__":
    # Query và Target dưới dạng danh sách các số nguyên (mã ASCII)
    nucleotide_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    query = [nucleotide_to_idx[c] for c in "ACGT"]  # Chuyển thành [0,1,2,3]
    target = [nucleotide_to_idx[c] for c in "AGCT"]
    
    # Ma trận điểm đơn giản: match = 2, mismatch = -1
    mat = [
        2, -1, -1, -1,  # A
        -1, 2, -1, -1,  # C
        -1, -1, 2, -1,  # G
        -1, -1, -1, 2   # T
    ]
    
    # Các tham số khác
    qlen = len(query)
    tlen = len(target)
    m = 4  # Số loại nucleotide (A, C, G, T)
    o_del = 2  # Gap opening penalty cho deletion
    e_del = 1  # Gap extension penalty cho deletion
    o_ins = 2  # Gap opening penalty cho insertion
    e_ins = 1  # Gap extension penalty cho insertion
    w = 2  # Chiều rộng dải
    end_bonus = 0  # Không có phần thưởng khi kết thúc
    zdrop = 5  # Ngưỡng dừng sớm
    h0 = 0  # Điểm bắt đầu
    
    # Biến kết quả
    _qle = [0]
    _tle = [0]
    _gtle = [0]
    _gscore = [0]
    _max_off = [0]
    
    max_score = ksw_extend2(
        qlen, query, tlen, target, m, mat, 
        o_del, e_del, o_ins, e_ins, w, end_bonus, 
        zdrop, h0, _qle, _tle, _gtle, _gscore, _max_off
    )
    
    # In kết quả
    print("Max Score:", max_score)
    print("Query End:", _qle[0])
    print("Target End:", _tle[0])
    print("GScore:", _gscore[0])
    print("Max Offset:", _max_off[0])
