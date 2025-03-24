def kroundup32(x):
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    return x + 1

class kstring_t:
    def __init__(self):
        self.l = 0
        self.m = 0
        self.s = None

if __name__ == "__main__":
    pass