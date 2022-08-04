from BigInt import BigInt


# secp224k1
class Point:
    """p = BigInt(2) ** 224 - 2 ** 32 - 2 ** 12 - 2 ** 11 - 2 ** 9 - 2 ** 7 - 2 ** 4 - 3
    a = BigInt(0)
    b = BigInt(5)
    """
    p = BigInt(17)
    a = BigInt(2)
    b = BigInt(2)

    @classmethod
    def inv(cls, n):
        return pow(n, cls.p - 2, cls.p)

    def __init__(self, x, y):
        self.x = x % Point.p
        self.y = y % Point.p

    def __repr__(self):
        return f'Point({self.x}, {self.y})'

    def __neg__(self):
        return Point(self.x, -self.y)

    def __sub__(self, other):
        return self + (-other)

    def __add__(self, other):
        if self == Point_E:
            return other.copy()
        if other == Point_E:
            return self.copy()
        if self == -other:
            return Point_E.copy()

        if self == other:
            k = (self.x * self.x * 3 + Point.a) * Point.inv(self.y * 2)
        else:
            k = (other.y - self.y) * Point.inv(other.x - self.x)

        k %= Point.p
        x = k * k - self.x - other.x
        y = k * (self.x - x) - self.y
        return Point(x, y)

    def __mul__(self, n):
        result = copy.deepcopy(Point_E)
        point = copy.deepcopy(self)
        n_tmp = copy.deepcopy(n)
        while n_tmp > 0:
            if n_tmp % 2 == 1:
                result += point
            point += point
            n_tmp //= 2
        return result

    def __rmul__(self, n):
        return self * n

    def __eq__(self, other):
        return (self.x, self.y) == (other.x, other.y)

    def copy(self):
        return Point(self.x, self.y)


#Point_G = Point(BigInt('16983810465656793445178183341822322175883642221536626637512293983324'),
#                BigInt('13272896753306862154536785447615077600479862871316829862783613755813'))
Point_G = Point(BigInt(5), BigInt(1))
Point_E = Point(BigInt(0), BigInt(0))


