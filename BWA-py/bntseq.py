# khai báo các hằng số
class bntann1_t:
    """
    Chứa thông tin về một chuỗi trong bộ genome.
    """
    def __init__(self, offset, length, n_ambs, gi, is_alt, name, anno):
        self.offset = offset       # Vị trí bắt đầu của chuỗi
        self.length = length       # Độ dài chuỗi
        self.n_ambs = n_ambs       # Số lượng nucleotide không xác định (N)
        self.gi = gi               # GenBank identifier
        self.is_alt = is_alt       # Cờ đánh dấu bản thay thế
        self.name = name           # Tên chuỗi (str hoặc None)
        self.anno = anno           # Thông tin chú thích bổ sung (str hoặc None)

class bntamb1_t:
    """ 
    Chứa thông tin về một khoảng nucleotide không xác định (N).
    """
    def __init__(self, offset, length, amb):
        self.offset = offset       # Vị trí bắt đầu của đoạn không xác định
        self.length = length       # Độ dài của đoạn không xác định
        self.amb = amb             # Ký tự đại diện cho đoạn không xác định ('N', 'X', ...)

class bntseq_t:
    """ 
    Chứa toàn bộ thông tin về bộ genome.
    """
    def __init__(self, l_pac, n_seqs, seed, anns, n_holes, ambs, fp_pac):
        self.l_pac = l_pac         # Tổng chiều dài bộ genome (ở dạng packed)
        self.n_seqs = n_seqs       # Số lượng chuỗi trong bộ genome
        self.seed = seed           # Giá trị seed ngẫu nhiên
        self.anns = anns           # Danh sách thông tin chuỗi (mảng BntAnn1)
        self.n_holes = n_holes     # Số lượng đoạn không xác định (N)
        self.ambs = ambs           # Danh sách đoạn không xác định (mảng BntAmb1)
        self.fp_pac = fp_pac       # Đường dẫn đến file chứa genome (str hoặc None)



def bns_pos2rid(bns: bntseq_t, pos_f):
    """
    Xác định ID của chuỗi (rid) chứa vị trí `pos_f` trong bộ genome.

    Parameters:
    - bns (BntSeq): Đối tượng chứa thông tin về bộ genome.
    - pos_f (int): Vị trí cần tìm trong bộ genome.

    Returns:
    - (int): ID của chuỗi chứa `pos_f`, hoặc -1 nếu `pos_f` nằm ngoài phạm vi genome.
    """
    # Nếu vị trí `pos_f` lớn hơn hoặc bằng độ dài genome, trả về -1 (nằm ngoài phạm vi hợp lệ)
    if pos_f >= bns.l_pac:
        return -1

    # Khởi tạo các biến cho tìm kiếm nhị phân
    left, mid, right = 0, 0, bns.n_seqs

    # Thực hiện tìm kiếm nhị phân để tìm chuỗi chứa `pos_f`
    while left < right:
        mid = (left + right) // 2

        # Kiểm tra xem `pos_f` có nằm trong chuỗi `mid` không
        if pos_f >= bns.anns[mid].offset:
            if mid == bns.n_seqs - 1 or pos_f < bns.anns[mid + 1].offset:
                break
            left = mid + 1
        else:
            right = mid

    # Trả về ID của chuỗi tìm được
    return mid


def bns_depos(bns: bntseq_t, pos):
    """
    Chuyển đổi vị trí pos trong bộ genome thành vị trí tương ứng theo hướng xuôi hoặc ngược.

    Tham số:
    - bns (bntseq_t): Đối tượng chứa thông tin về bộ genome.
    - pos (int): Vị trí trong bộ genome (đơn vị base).

    Kết quả trả về:
    - int: Vị trí đã được chuyển đổi.
    - bool: Giá trị True nếu vị trí nằm ở nửa sau (tức là hướng ngược), False nếu không.
    """
    # Kiểm tra xem vị trí pos có nằm trong nửa sau của genome không
    is_rev = pos >= bns.l_pac  

    # Nếu is_rev = True (pos nằm trong nửa sau), chuyển đổi vị trí theo công thức:
    # (bns.l_pac * 2) - 1 - pos
    # Nếu is_rev = False, giữ nguyên vị trí pos
    return (bns.l_pac * 2) - 1 - pos if is_rev else pos, is_rev

def _get_pac(pac, l):
    index = l >> 2
    if index >= len(pac):
        return 0  # Trả về giá trị mặc định thay vì gây lỗi
    return (pac[index] >> ((~l & 3) << 1)) & 3


def bns_get_seq(l_pac, pac, beg, end):
    """
    Lấy trình tự nucleotide từ bộ genome, có thể ở chiều xuôi hoặc ngược.

    Tham số:
    - l_pac: Tổng chiều dài genome (ở dạng packed).
    - pac: Mảng chứa dữ liệu genome đã nén.
    - beg, end: Chỉ số bắt đầu và kết thúc.

    Trả về:
    - tuple (seq, length): `seq` là bytearray chứa trình tự, `length` là độ dài chuỗi.
    """
    if end < beg:
        beg, end = end, beg  # Swap nếu end nhỏ hơn beg

    # Giới hạn beg và end trong phạm vi hợp lệ
    end = min(end, l_pac << 1)
    beg = max(beg, 0)

    # Kiểm tra xem beg và end có hợp lệ hay không
    if beg >= l_pac or end <= l_pac:
        # Tạo bộ nhớ để lưu trình tự genome
        length = end - beg
        seq = bytearray(length)  # Thay cho malloc()

        # Xử lý trường hợp beg nằm ở reverse strand
        if beg >= l_pac:  # reverse strand
            # Chuyển đổi chỉ số từ reverse về forward
            beg_f = (l_pac << 1) - 1 - end
            end_f = (l_pac << 1) - 1 - beg

            # Lặp qua đoạn này theo chiều ngược lại và đảo giá trị nucleotide
            for l, k in enumerate(range(end_f - 1, beg_f - 1, -1)):  # Sửa range()
                seq[l] = 3 - _get_pac(pac, k)
        
        # Xử lý trường hợp forward strand
        else:  # forward strand
            for l, k in enumerate(range(beg, end)):
                seq[l] = _get_pac(pac, k)

        return seq, length  # Trả về cả chuỗi và độ dài
    return bytearray(), 0  # Trả về mảng rỗng và độ dài 0 nếu không hợp lệ


def bns_fetch_seq(bns: bntseq_t, pac, beg, mid, end):
    """
    Lấy một đoạn trình tự nucleotide (seq) từ bộ genome.
    
    Tham số:    
        - bns	bntseq_t	Đối tượng chứa thông tin về bộ genome.
        - pac	bytearray hoặc list	Mảng chứa dữ liệu genome đã nén.
        - beg	int	Vị trí bắt đầu lấy trình tự genome.
        - mid	int	Vị trí giữa (điểm tham chiếu để xác định rid).
        - end	int	Vị trí kết thúc lấy trình tự genome.

    Trả về:
        - seq	bytearray	Trình tự nucleotide lấy được.
        - beg	int	Vị trí bắt đầu đã được điều chỉnh theo phạm vi hợp lệ.
        - end	int	Vị trí kết thúc đã được điều chỉnh theo phạm vi hợp lệ.
        - rid	int	ID của chuỗi chứa vị trí mid.
    """
    if end < beg:
        beg, end = end, beg # Swap nếu end nhỏ hơn beg

    # Kiểm tra mid nằm trong khoảng [beg, end]
    assert beg <= mid < end

    # Xác định rid và kiểm tra strand
    pos_f, is_rev = bns_depos(bns, mid)
    rid = bns_pos2rid(bns, pos_f)

    # Xác định far_beg và far_end
    far_beg = bns.anns[rid].offset
    far_end = far_beg + bns.anns[rid].length

    # Nếu là reverse strand, chuyển đổi far_beg và far_end
    if is_rev:  # Flip to the reverse strand
        far_beg, far_end = (bns.l_pac << 1) - far_end, (bns.l_pac << 1) - far_beg

    # Giới hạn beg và end trong phạm vi far_beg và far_end
    beg = max(beg, far_beg)
    end = min(end, far_end)

    # Lấy trình tự genome
    seq, length = bns_get_seq(bns.l_pac, pac, beg, end)

    # Kiểm tra lỗi
    if seq is None or length != end - beg:
        raise AssertionError(f"Error in {__name__}: begin={beg}, mid={mid}, end={end}, len={len(seq)}, seq={seq}, rid={rid}, far_beg={far_beg}, far_end={far_end}")
    
    # Kiểm tra lại bằng assert
    assert seq is not None and length == end - beg  # Assertion failure should never happen

    return seq, beg, end, rid # Return sequence, begin, end, and reference ID

#===============================================================================
# Test
if __name__ == "__main__":
    # Tạo một số đối tượng bntann1_t đại diện cho các chuỗi genome
    anns = [
        bntann1_t(offset=0, length=100, n_ambs=0, gi=1, is_alt=False, name="seq1", anno=None),
        bntann1_t(offset=100, length=150, n_ambs=0, gi=2, is_alt=False, name="seq2", anno=None),
        bntann1_t(offset=250, length=200, n_ambs=0, gi=3, is_alt=False, name="seq3", anno=None)
    ]
    # Khởi tạo đối tượng bntseq_t
    bns = bntseq_t(l_pac=450, n_seqs=len(anns), seed=42, anns=anns, n_holes=0, ambs=[], fp_pac=None)
    # Danh sách các vị trí cần kiểm tra
    test_positions = [50, 120, 260, 400, 500]  # Vị trí 500 nằm ngoài phạm vi

    print(bns)
    # Kiểm tra hàm bns_pos2rid
    for pos in test_positions:
        rid = bns_pos2rid(bns, pos)
        print(f"Position {pos} is in sequence ID: {rid}")

    # Kiểm tra hàm bns_depos
    test_positions_depos = [100, 200, 450, 500, 700]
    for pos in test_positions_depos:
        pos_f, is_rev = bns_depos(bns, pos)
        # Nếu is_rev = False, nghĩa là vị trí nằm trong nửa đầu của genome, tức là theo chiều xuôi (forward strand).
        # Nếu is_rev = True, nghĩa là vị trí nằm trong nửa sau của genome, tức là theo chiều ngược (reverse strand). 
        print(f"Original position: {pos}, Forward position: {pos_f}, Is reverse: {is_rev}")

    # Kiểm tra hàm bns_get_seq
    # pac là một mảng byte (bytearray) chứa dữ liệu genome đã được mã hóa (packed).
    # Mỗi nucleotide được lưu trong 2 bit (vì có 4 loại: A, T, C, G → 2-bit là đủ).
    # Mỗi byte chứa 4 nucleotide, do đó pac có thể lưu tối đa 2 x len(pac) nucleotide.
    # pac = [0b00011011, 0b11001000] có 2 byte, tức là lưu được 2 × 4 = 8 nucleotide.
    # 00 01 10 11 → A C G T
    # 11 00 10 00 → T A G A
    pac = bytearray([0b00011011, 0b11001000])  # Dữ liệu giả lập
    beg, end = 0, 8 # Lấy toàn bộ 8 nucleotide
    print("-" * 10, "\n", bns.l_pac, pac, beg, end, "-" * 10)
    seq, length = bns_get_seq(bns.l_pac, pac, beg, end)
    print(f"Extracted sequence (bns_get_seq): {list(seq)}, Length: {length}")


    # Kiểm tra hàm bns_fetch_seq
    # lấy trình tự từ vị trí 10 đến 90
    # mid = 50 là điểm giữa, dùng để xác định trình tự (rid) nào chứa đoạn này.
    beg, mid, end = 10, 50, 90
    seq, new_beg, new_end, rid = bns_fetch_seq(bns, pac, beg, mid, end)
    print(f"Fetched sequence (bns_fetch_seq): {list(seq)}, Length: {len(seq)}, Beg: {new_beg}, End: {new_end}, RID: {rid}")