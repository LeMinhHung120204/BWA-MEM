import os
import mmap
import struct
from bwa import *

BWA_CTL_SIZE = 1024  # Giả sử kích thước bộ nhớ điều khiển

def bwa_idx_load_from_shm(hint: str):
    if not hint:
        return None

    # Lấy `name` từ `hint`
    name = hint.split("/")[-1]

    try:
        # Mở shared memory "/bwactl" để đọc
        fd = os.open("/dev/shm/bwactl", os.O_RDONLY)
        shm = mmap.mmap(fd, BWA_CTL_SIZE, mmap.PROT_READ, mmap.MAP_SHARED)
    except OSError:
        return None

    # Đọc số lượng index trong shared memory
    cnt = struct.unpack("H", shm[:2])[0]  # uint16_t (2 bytes)
    if cnt == 0:
        return None

    # Tìm `name` trong shared memory
    offset = 4
    l_mem = 0
    for _ in range(cnt):
        l_mem = struct.unpack("q", shm[offset:offset+8])[0]  # int64_t (8 bytes)
        offset += 8
        end = shm.find(b'\x00', offset)  # Tìm byte NULL kết thúc chuỗi
        idx_name = shm[offset:end].decode()
        offset = end + 1
        if idx_name == name:
            break
    else:
        return None  # Không tìm thấy

    # Tạo đường dẫn shared memory cho index
    path = f"/dev/shm/bwaidx-{name}"

    try:
        # Mở shared memory index
        fd = os.open(path, os.O_RDONLY)
        shm_idx = mmap.mmap(fd, l_mem, mmap.PROT_READ, mmap.MAP_SHARED)
    except OSError:
        return None

    # Cấp phát và khởi tạo `bwaidx_t`
    idx = bwaidx_t()
    bwa_mem2idx(l_mem, shm_idx, idx)
    idx.is_shm = 1

    return idx
