import sys

def err_fread_noeof(fp, size, count):
    """Đọc `count` phần tử, mỗi phần tử có kích thước `size`, từ tệp fp."""
    data = fp.read(size * count)
    if len(data) != size * count:
        raise IOError("Unexpected end of file or read error")
    return data