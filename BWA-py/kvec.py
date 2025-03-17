class kvec_t:
    def __init__(self, type=int):
        self.n = 0  # Số phần tử hiện tại
        self.m = 0  # Kích thước bộ nhớ được cấp phát
        self.a = []  # Danh sách lưu trữ nội bộ
        self.type = type  # Kiểu dữ liệu

    def kv_init(self):
        """Khởi tạo vector."""
        self.n = 0
        self.m = 0
        self.a = []

    def kv_roundup32(self, x):
        """Làm tròn lên x thành số nguyên gần nhất là lũy thừa của 2."""
        x -= 1
        x |= x >> 1
        x |= x >> 2
        x |= x >> 4
        x |= x >> 8
        x |= x >> 16
        x += 1
        return x

    def kv_destroy(self):
        """Giải phóng bộ nhớ."""
        self.a = []
        self.n = 0
        self.m = 0

    def kv_A(self, i):
        """Lấy phần tử tại chỉ mục cụ thể."""
        if i >= self.n:
            raise IndexError("Index out of bounds")
        return self.a[i]
    
    def kv_pop(self):
        """Xóa phần tử cuối cùng khỏi mảng."""
        if self.n == 0:
            raise IndexError("Pop from empty kvec_t")
        self.n -= 1
        return self.a[self.n]
    
    def kv_size(self):
        """Trả về số lượng phần tử hiện tại."""
        return self.n

    def kv_max(self):
        """Trả về kích thước bộ nhớ tối đa được cấp phát."""
        return self.m

    def kv_resize(self, new_size):
        """Thay đổi kích thước bộ nhớ để chứa ít nhất new_size phần tử."""
        if new_size <= self.m:
            return
        self.m = self.kv_roundup32(new_size)  # Làm tròn kích thước lên lũy thừa của 2
        new_array = [self.type() for _ in range(self.m)]  # Tạo danh sách mới với kích thước tăng lên
        for i in range(self.n):
            new_array[i] = self.a[i]  # Sao chép phần tử cũ vào danh sách mới
        self.a = new_array  # Thay thế danh sách cũ bằng danh sách mới

    def kv_copy(self, other):
        """Sao chép phần tử từ một kvec_t khác."""
        if not isinstance(other, kvec_t):
            raise TypeError("Incompatible kvec_t types")
        self.kv_resize(other.n)
        for i in range(other.n):
            self.a[i] = other.a[i]  # Sao chép phần tử
        self.n = other.n

    def kv_push(self, value):
        """Thêm một phần tử vào mảng."""
        if self.n == self.m:
            self.kv_resize(self.n + 1)
        if self.n < len(self.a):
            self.a[self.n] = value
        else:
            self.a.append(value)
        self.n += 1

    def kv_pushp(self):
        """Thêm một phần tử và trả về con trỏ đến phần tử mới."""
        if self.n == self.m:
            self.kv_resize(self.n + 1)
        if self.n < len(self.a):
            self.a[self.n] = self.type()
        else:
            self.a.append(self.type())
        self.n += 1
        return self.a[self.n - 1]

    def kv_a(self, i):
        """Truy cập hoặc mở rộng mảng tại vị trí i."""
        if self.m <= i:
            self.m = self.n = i + 1
            self.m = self.kv_roundup32(self.m)
            new_array = [self.type() for _ in range(self.m)]
            for j in range(len(self.a)):
                new_array[j] = self.a[j]
            self.a = new_array
        elif self.n <= i:
            self.n = i + 1
        return self.a[i]