from ElGamal import *


def log_bsgs(a, g, mod):
    # solving for a^res = g % mod
    a = BigInt(a) if type(a) == int else a
    g = BigInt(g) if type(g) == int else g
    mod = BigInt(mod) if type(mod) == int else mod
    a %= mod
    g %= mod
    m = BigInt.sqrt(mod) + 1

    a0 = a.powWithMod(m, mod)
    tmp = a0
    l1 = {1: a0}
    for j in range(2, int(m + 1)):
        tmp = tmp.multWithMod(a0, mod)
        l1.update({j: tmp})

    tmp = g.multWithMod(a.powWithMod(mod - m, mod), mod)

    for i in range(int(m) - 1, -1, -1):
        for key, value in l1.items():
            if value == tmp:
                return (m * key + i) % mod
        tmp = tmp.multWithMod(a, mod)

    return None


def factorize_pollard(num):
    num = BigInt(num) if type(num) == int else num
    func = lambda a: (a * a + 1) % num
    x = BigInt(2)
    y = BigInt(2)

    res = []
    attempts = 1000
    while num != 1 and attempts > 0:
        attempts -= 1
        x = func(x)
        y = func(func(y))
        d = BigInt.gcd(num, abs(x - y))
        if d != 1:
            res.append(d)
            num //= d
            x = BigInt(2)
            y = BigInt(2)
            attempts = 1000

    if num > 1:
        res.append(num)
    return res


# euler's totient function
def phi(num):
    res = BigInt(num) if type(num) == int else num
    for i in range(2, int(BigInt.sqrt(num)) + 1):
        if num % i == 0:
            while num % i == 0:
                num //= i
            res -= res // i

    if num > 1:
        res -= res // num
    return res


# mobius function
def mu(num):
    num = BigInt(num) if type(num) == int else num
    if num == 1:
        return 1

    res = BigInt(0)
    for i in range(2, int(num + 1)):
        if num % i == 0 and is_prime(i):
            if num % (BigInt(i) * BigInt(i)) == 0:
                return 0
            res += 1

    if res % 2 == 0:
        return 1
    return -1


def legendre(a, p):
    a = BigInt(a) if type(a) == int else a
    p = BigInt(p) if type(p) == int else p

    if p <= 2 or not is_prime(p):
        return None
    if a % p == 0:
        return 0
    if a.powWithMod((p - 1)//2, p) == 1:
        return 1
    return -1


def jacobi(a, n):
    a = BigInt(a) if type(a) == int else a
    n = BigInt(n) if type(n) == int else n

    if n <= 0 or n % 2 == 0:
        return None

    if BigInt.gcd(a, n) != 1:
        return 0

    a %= n
    t = 1
    while a != 0:
        while a % 2 == 0:
            a /= 2
            r = n % 8
            if r == 3 or r == 5:
                t = -t

        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            t = -t
        a %= n
    if n == 1:
        return t
    return 0


def is_prime(num):
    num = BigInt(num) if type(num) == int else num
    if num == 2:
        return True
    if num % 2 == 0:
        return False
    for i in range(3, int(BigInt.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True


def sqrt_cipolla(num, p):
    num = BigInt(num) if type(num) == int else num
    p = BigInt(p) if type(p) == int else p

    if legendre(num, p) != 1:
        return None

    a = BigInt(0)
    while legendre((a*a - num) % p, p) != -1:
        a += 1
    w = (a*a - num) % p

    def mult_ext(a, b):
        return (a[0]*b[0] % p + a[1]*b[1]*w % p) % p, (a[0]*b[1] % p + a[1]*b[0] % p) % p

    res = (BigInt(1), BigInt(0))
    a0 = (a, BigInt(1))
    pw = (p + 1) / 2
    while pw != 0:
        if pw % 2 != 0:
            res = mult_ext(res, a0)
        a0 = mult_ext(a0, a0)
        pw //= 2

    if res[1] != 0 or res[0]*res[0] % p != num:
        return None

    return res[0], p - res[0]


def is_prime_solovay_strassen(num, accuracy):
    if num < 2:
        return False
    if num == 2:
        return True

    if num % 2 == 0:
        return False

    for i in range(1, accuracy):
        a = BigInt.rand(num - 1 )
        if BigInt.gcd(a, num) > 1:
            return False
        if a.powWithMod((num - 1) / 2, num) != jacobi(a, num) % num:
            return False

    return True

# Decomment any block to test the functions

"""
# 11^x = 13 mod 127
b = log_bsgs(7, 156, 1237)
print(b)

# 87619876495133
s = factorize_pollard(12374869128373)
for i in s:
    print(i)

s = phi(7197107)
print(s)

s = mu(19)
print(s)
"""
a = BigInt(-68) % 821
print(a)
s = legendre(a,821)
print(s)
"""

s = jacobi(1, 2)
print(s)

print(BigInt(-6) % 13)


s = sqrt_cipolla(157, 1237)
print(s[0], s[1])

s = is_prime_solovay_strassen(12374869128373, 100)
print(s)

message = 10*Point_G

print(f"message: {message}")

sk, pk = ElGamal.generate_key()
ciphertext = ElGamal.encrypt(pk, message)
received_message = ElGamal.decrypt(sk, ciphertext)

print(f"received_message: {received_message}")
"""
