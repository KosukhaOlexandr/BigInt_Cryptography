import copy

from ElGamal import *


def log_bsgs(a_log, g, mod):
    # solving for a^res = g % mod
    a_tmp = copy.deepcopy(a_log)
    g_tmp = copy.deepcopy(g)
    mod_tmp = copy.deepcopy(mod)
    a_tmp = BigInt(a_tmp) if type(a_tmp) == int else a_tmp
    g_tmp = BigInt(g_tmp) if type(g_tmp) == int else g_tmp
    mod_tmp = BigInt(mod_tmp) if type(mod_tmp) == int else mod_tmp
    a_tmp %= mod_tmp
    g_tmp %= mod_tmp
    m = BigInt.sqrt(mod_tmp) + 1

    a0 = a_tmp.pow_with_mod(m, mod_tmp)
    tmp = a0
    l1 = {1: a0}
    for j in range(2, int(m + 1)):
        tmp = tmp.mult_with_mod(a0, mod_tmp)
        l1.update({j: tmp})

    tmp = g_tmp.mult_with_mod(a_tmp.pow_with_mod(mod_tmp - m, mod_tmp), mod_tmp)

    for i in range(int(m) - 1, -1, -1):
        for key, value in l1.items():
            if value == tmp:
                return (m * key + i) % mod_tmp
        tmp = tmp.mult_with_mod(a_tmp, mod_tmp)

    return None


def factorize_pollard(num):
    num_tmp = copy.deepcopy(num)
    num_tmp = BigInt(num_tmp) if type(num_tmp) == int else num_tmp
    func = lambda arg: (arg * arg + 1) % num_tmp
    x = BigInt(2)
    y = BigInt(2)

    res = []
    attempts = 1000
    while num_tmp != 1 and attempts > 0:
        attempts -= 1
        x = func(x)
        y = func(func(y))
        d = BigInt.gcd(num_tmp, abs(x - y))
        if d != 1:
            res.append(d)
            num_tmp //= d
            x = BigInt(2)
            y = BigInt(2)
            attempts = 1000

    if num_tmp > 1:
        res.append(num_tmp)
    return res


# euler's totient function
def phi(num):
    res = BigInt(num) if type(num) == int else copy.deepcopy(num)
    num_tmp = copy.deepcopy(num)
    for i in range(2, int(BigInt.sqrt(num_tmp)) + 1):
        if num_tmp % i == 0:
            while num_tmp % i == 0:
                num_tmp //= i
            res -= res // i

    if num_tmp > 1:
        res -= res // num_tmp
    return res


# mobius function
def mu(num):
    num_tmp = BigInt(num) if type(num) == int else copy.deepcopy(num)

    if num_tmp == 1:
        return 1

    res = BigInt(0)
    for i in range(2, int(num_tmp + 1)):
        if num_tmp % i == 0 and is_prime(i):
            if num_tmp % (BigInt(i) * BigInt(i)) == 0:
                return 0
            res += 1

    if res % 2 == 0:
        return 1
    return -1


def legendre(a_l, p_l):
    a_tmp = BigInt(a_l) if type(a_l) == int else copy.deepcopy(a_l)
    p_tmp = BigInt(p_l) if type(p_l) == int else copy.deepcopy(p_l)

    if p_tmp <= 2 or not is_prime(p_tmp):
        return None
    if a_tmp % p_tmp == 0:
        return 0
    if a_tmp.pow_with_mod((p_tmp - 1) // 2, p_tmp) == 1:
        return 1
    return -1


def jacobi(a_j, n_j):
    a_tmp = BigInt(a_j) if type(a_j) == int else a_j
    n_tmp = BigInt(n_j) if type(n_j) == int else n_j

    if n_tmp <= 0 or n_tmp % 2 == 0:
        return None

    if BigInt.gcd(a_tmp, n_tmp) != 1:
        return 0

    a_tmp %= n_tmp
    t = 1
    while a_tmp != 0:
        while a_tmp % 2 == 0:
            a_tmp /= 2
            r = n_tmp % 8
            if r == 3 or r == 5:
                t = -t

        a_tmp, n_tmp = n_tmp, a_tmp
        if a_tmp % 4 == 3 and n_tmp % 4 == 3:
            t = -t
        a_tmp %= n_tmp
    if n_tmp == 1:
        return t
    return 0


def is_prime(num):
    num_tmp = BigInt(num) if type(num) == int else copy.deepcopy(num)
    if num_tmp == 2:
        return True
    if num_tmp % 2 == 0:
        return False
    for i in range(3, int(BigInt.sqrt(num_tmp)) + 1):
        if num_tmp % i == 0:
            return False
    return True


def sqrt_cipolla(num, p):
    num_tmp = BigInt(num) if type(num) == int else copy.deepcopy(num)
    p_tmp = BigInt(p) if type(p) == int else copy.deepcopy(p)

    if legendre(num_tmp, p_tmp) != 1:
        return None

    a_s = BigInt(0)
    while legendre((a_s * a_s - num_tmp) % p_tmp, p_tmp) != -1:
        a_s += 1
    w = (a_s * a_s - num_tmp) % p_tmp

    def mult_ext(a, b):
        return (a[0] * b[0] % p_tmp + a[1] * b[1] * w % p_tmp) % p_tmp, \
               (a[0] * b[1] % p_tmp + a[1] * b[0] % p_tmp) % p_tmp

    res = (BigInt(1), BigInt(0))
    a0 = (a, BigInt(1))
    pw = (p_tmp + 1) / 2
    while pw != 0:
        if pw % 2 != 0:
            res = mult_ext(res, a0)
        a0 = mult_ext(a0, a0)
        pw //= 2

    if res[1] != 0 or res[0] * res[0] % p_tmp != num_tmp:
        return None

    return res[0], p_tmp - res[0]


def is_prime_solovay_strassen(num, accuracy):
    if num < 2:
        return False
    if num == 2:
        return True

    if num % 2 == 0:
        return False

    for i in range(1, accuracy):
        a_p = BigInt.rand(num - 1)
        if BigInt.gcd(a_p, num) > 1:
            return False
        if a.pow_with_mod((num - 1) / 2, num) != jacobi(a, num) % num:
            return False

    return True


# Uncomment any block to test the functions


# 11^x = 13 mod 127
b = log_bsgs(11, 13, 127)
a = BigInt(11)
print(a.pow_with_mod(107, 127))
print(b)

# 87619876495133
s = factorize_pollard(12374869128373)
for i in s:
    print(i)

s = phi(7197107)
print(s)

s = mu(19)
print(s)

a = BigInt(-68) % 821
print(a)
s = legendre(a, 821)
print(s)


s = jacobi(1, 2)
print(s)

print(BigInt(-6) % 13)


s = sqrt_cipolla(157, 1237)
print(s[0], s[1])

s = is_prime_solovay_strassen(12374869128373, 100)
print(s)
"""
message = 10*Point_G

print(f"message: {message}")

sk, pk = ElGamal.generate_key()
ciphertext = ElGamal.encrypt(pk, message)
received_message = ElGamal.decrypt(sk, ciphertext)

print(f"received_message: {received_message}")
"""
