from BigInt import BigInt
from ElCurve import Point, Point_G


class ElGamal:
    @staticmethod
    def encrypt(pk, message: Point):
        y = BigInt.rand(Point.p)
        secret = pk * y
        c1 = Point_G * y
        c2 = message + secret
        return c1, c2

    @staticmethod
    def decrypt(sk, ciphertext) -> Point:
        c1, c2 = ciphertext
        secret = c1 * sk
        message = c2 - secret
        return message

    @staticmethod
    def generate_key():
        sk = BigInt.rand(Point.p)
        pk = Point_G * sk
        return sk, pk
