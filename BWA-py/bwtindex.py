import sys

def bwt_bwtupdate_core(bwt):
    n_occ = (bwt['seq_len'] + bwt['occ_intv'] - 1) // bwt['occ_intv'] + 1
    bwt['bwt_size'] += n_occ * 4  # the new size
    buf = [0] * bwt['bwt_size']  # will be the new bwt
    c = [0, 0, 0, 0]
    k = 0
    for i in range(bwt['seq_len']):
        if i % bwt['occ_intv'] == 0:
            buf[k:k+4] = c
            k += 4
        if i % 16 == 0:
            buf[k] = bwt['bwt'][i // 16]
            k += 1
        c[bwt_B00(bwt, i)] += 1
    # the last element
    buf[k:k+4] = c
    assert k + 4 == bwt['bwt_size'], "inconsistent bwt_size"
    # update bwt
    bwt['bwt'] = buf

def bwt_B00(bwt, k):
    return (bwt['bwt'][k >> 4] >> ((~k & 0xf) << 1)) & 3

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: bwtgen <in.pac> <out.bwt>")
        sys.exit(1)