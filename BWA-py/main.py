import sys
from fastmap import *

def bwa_index(args):
    # hàm này bên fastmap chưa viết
    pass

def main():
    if len(sys.argv) < 2:
        print("Usage: script.py <command>")
        sys.exit(1)
    
    command = sys.argv[1]
    args = sys.argv[2:]
    
    if command == "index":
        ret = bwa_index(args)
    elif command == "mem":
        ret = main_mem(args)
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
    
    sys.exit(ret if ret is not None else 0)

if __name__ == "__main__":
    main()
