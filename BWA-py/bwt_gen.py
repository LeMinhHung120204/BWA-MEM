import sys

def bwt_bwtgen2(fn_pac, fn_bwt, block_size):
    bwtInc = BWTIncConstructFromPacked(fn_pac, block_size, block_size)
    print(f"[bwt_gen] Finished constructing BWT in {bwtInc['numberOfIterationDone']} iterations.")
    BWTSaveBwtCodeAndOcc(bwtInc['bwt'], fn_bwt, 0)
    BWTIncFree(bwtInc)

def bwt_bwtgen(fn_pac, fn_bwt):
    bwt_bwtgen2(fn_pac, fn_bwt, 10000000)

def BWTIncConstructFromPacked(fn_pac, initial_max_build_size, inc_max_build_size):
    pass

def BWTSaveBwtCodeAndOcc(bwt, fn_bwt, _):
    pass

def BWTIncFree(bwtInc):
    pass

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: bwtgen <in.pac> <out.bwt>")
        sys.exit(1)
    bwt_bwtgen(sys.argv[1], sys.argv[2])