import copy
import random


class BigInt:
    BASE = 10
    # working with numbers with base BASE
    def __init__(self, init_value=None):

        self.is_pos = True  # 0 is positive too
        self.digits = []

        if init_value is not None:

            assert type(init_value) == str or type(init_value) == int, 'Initial_value must be int or str'
            # should be changed if BASE != 10
            init_value = str(init_value)
            right_border = 0
            if init_value[0] == '-':
                self.is_pos = False
                right_border = 1
            elif init_value[0] == '+':
                right_border = 1

            for x in range(len(init_value) - 1, right_border - 1, -1):
                self.digits.append(int(init_value[x]))

    def __str__(self):
        res = ''
        if not self.is_pos:
            res += '-'
        for i in range(len(self.digits) - 1, -1, -1):
            res += str(self.digits[i])
        return res

    def __eq__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        return self.is_pos == other.is_pos and self.digits == other.digits

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if self.is_pos != other.is_pos:
            return not self.is_pos

        if len(self.digits) != len(other.digits):
            if self.is_pos:
                return len(self.digits) < len(other.digits)
            else:
                return len(self.digits) > len(other.digits)

        for i in range(len(self.digits) - 1, -1, -1):
            if self.digits[i] != other.digits[i]:
                if self.is_pos:
                    return self.digits[i] < other.digits[i]
                else:
                    return self.digits[i] > other.digits[i]

        return False

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def __neg__(self):
        tmp = copy.deepcopy(self)
        tmp.is_pos = not tmp.is_pos
        return tmp

    def __pos__(self):
        return self

    def __add__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if self.is_pos == other.is_pos:
            res = BigInt()
            res.is_pos = self.is_pos
            extra = 0
            for i in range(max(len(self.digits), len(other.digits))):
                if i >= len(self.digits):
                    a = 0
                else:
                    a = self.digits[i]

                if i >= len(other.digits):
                    b = 0
                else:
                    b = other.digits[i]

                fullSum = a + b + extra
                res.digits.append(fullSum % self.BASE)
                extra = fullSum // self.BASE

            if extra > 0:
                res.digits.append(extra)
            return res

        if self.is_pos:
            return self - (- other)
        return other - (- self)

    def __sub__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if self == other:
            return BigInt(0)

        if self.is_pos == other.is_pos:
            res = BigInt()
            res.is_pos = self.is_pos
            extra = 0
            if self < other:
                for i in range(max(len(self.digits), len(other.digits))):
                    if len(self.digits) <= i:
                        a = 0
                    else:
                        a = self.digits[i]

                    if len(other.digits) <= i:
                        b = 0
                    else:
                        b = other.digits[i]
                    if self.is_pos:
                        fullDif = b - a - extra
                    else:
                        fullDif = a - b - extra

                    if fullDif < 0:
                        fullDif += self.BASE
                        extra = 1
                    else:
                        extra = 0

                    res.digits.append(fullDif)
                res.is_pos = False
            else:
                res = -(other - self)
            res.__clearZeros()
            return res
            # clearZeros

        if self.is_pos:
            return self + (-other)
        return -(-self + other)

    def __mul__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if self == 0 or other == 0:
            return BigInt(0)

        if len(self.digits) < len(other.digits):
            return other * self

        blocks = []
        for i in range(len(other.digits)):
            pCarry = 0
            blocks.append(BigInt())

            for j in range(i):
                blocks[i].digits.append(0)

            for j in range(len(self.digits)):
                fullProd = self.digits[j] * other.digits[i] + pCarry
                blocks[i].digits.append(fullProd % self.BASE)
                pCarry = fullProd//self.BASE

            if pCarry > 0:
                blocks[i].digits.append(pCarry)

        res = BigInt(0)

        for j in range(len(blocks)):
            res += blocks[j]
        res.is_pos = (self.is_pos == other.is_pos)
        return res

    def __floordiv__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if other == 0:
            raise ZeroDivisionError
        if self == 0:
            return BigInt(0)
        if other == 1:
            return self
        if other == -1:
            return -self

        divRes = BigInt()
        divRes.is_pos = (self.is_pos == other.is_pos)
        subst = BigInt()
        subst.is_pos = True

        absOther = abs(other)
        if absOther is other:
            return 0

        i = len(self.digits) - 1
        firstStep = True

        while i >= 0:
            added = 0
            cond = True
            while cond:
                subst.digits.insert(0, self.digits[i])
                i -= 1
                added += 1

                subst.__clearZeros()
                if added > 1 and (not firstStep):
                    divRes.digits.insert(0, 0)
                cond = subst < absOther and i >= 0

            firstStep = False
            quot = subst.__simpleDiv(absOther)
            # divRes.digits.insert(0, quot)
            for j in range(len(quot.digits) -1, -1, -1):
                divRes.digits.insert(0, quot.digits[j])

            subst = subst - absOther*quot

            if subst == 0:
                subst.digits.remove(0)

        if (not self.is_pos) and (not other.is_pos) and \
            (subst != 0) and len(subst.digits) != 0:
            divRes = divRes + 1
        if (not self.is_pos) and other.is_pos and \
            (subst != 0) and len(subst.digits) != 0:
            divRes = divRes - 1
        if divRes.digits[0] == 0:
            divRes.is_pos = True
        return divRes

    def __mod__(self, other):
        other = BigInt.__fromIntToBigInt(other)
        if other == 0:
            raise ZeroDivisionError

        res = self - other * (self // other)
        if res == other:
            return BigInt(0)
        return res

    def __truediv__(self, other):
        return self // other

    def __pow__(self, power, modulo=None):
        if self == 0:
            return BigInt(1)
        if power == 0:
            return BigInt(1)
        if power == 1:
            return self
        if power % 2 == 0:
            temp = self**(power//2)
            return temp*temp
        return self**(power - 1)*self

    def powWithMod(self, power, mod):
        if self == 0:
            return BigInt(1)
        if power == 0:
            return BigInt(1)
        if power == 1:
            return self % mod
        if power % 2 == 0:
            temp = self.powWithMod(power//2, mod)
            return (temp*temp) % mod
        return (self.powWithMod(power - 1, mod)*self) % mod

    def __abs__(self):
        if self < 0:
            return -copy.deepcopy(self)
        return copy.deepcopy(self)

    def __int__(self):
        return int(self.__str__())



    def addWithMod(self, other, mod):
        return (self + other) % mod

    def subsWithMod(self, other, mod):
        return (self - other) % mod

    def multWithMod(self, other, mod):
        return (self * other) % mod

    def divWithMod(self, other, mod):
        return (self // other) % mod

    @staticmethod
    def sqrt(obj):
        obj = BigInt.__fromIntToBigInt(obj)
        if not obj.is_pos:
            raise ValueError("The number should be positive")
        if obj == 0:
            return BigInt(0)
        if obj < 4:
            return BigInt(1)

        res = BigInt.sqrt((obj - obj % 4) / 4)*2
        if obj < (res + 1)**2:
            return res
        return res + 1

    @staticmethod
    def gcd(a, b):
        a = BigInt.__fromIntToBigInt(a)
        b = BigInt.__fromIntToBigInt(b)
        a = abs(a)
        b = abs(b)

        while b != 0:
            a %= b
            a, b = b, a

        return a

    @staticmethod
    def rand(upper):
        upper = BigInt.__fromIntToBigInt(upper)
        max_int = (1 << 63) - 1
        if upper < max_int:
            res = BigInt(random.randint(1, upper))
        else:
            res = BigInt(random.randint(1, max_int))*upper / max_int
        res.__clearZeros()
        return res

    @staticmethod
    def __fromIntToBigInt(obj):
        assert type(obj) == BigInt or type(obj) == int, 'The element should be int or BigInt'
        if type(obj) == int:
            obj = BigInt(obj)
        return obj


    def __clearZeros(self):
        while len(self.digits) > 0 and self.digits[-1] == 0:
            self.digits.pop()

    def __simpleDiv(self, other):
        res = BigInt(0)
        a = self
        b = other
        while a >= b:
            a = a - b
            res = res + 1
        return res



