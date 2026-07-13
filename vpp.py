def pp_verify(p, A, x0):
    from math import gcd, isqrt

    if p < 5 or p%2 == 0:
        return False

    q = isqrt(p)                            # every proper prime divisor of p is bounded by q
    k = (q + 1 + isqrt(4*q)).bit_length()   # least integer k such that 2^k > q + 1 + 2*sqrt(q)

    # check that Montgomery curve By^2 = x^3 + A*x^2 + x is non-singular
    if gcd(A * A - 4, p) != 1:
        return False

    # Repeated doubling in projective Montgomery coordinates (X:Z)
    X, Z = x0 % p, 1
    Zprev = None

    # precompute C = (A + 2)/4 mod p
    C = ((A + 2) * ((p + 1) // 4 if p % 4 == 3 else (3*p +1) //4)) % p
    for i in range(k):
        Zprev = Z
        U,V = (X+Z)*(X+Z) % p, (X-Z)*(X-Z) % p
        XZ4 = U - V
        X,Z = U*V % p, XZ4*(V+C*XZ4) % p

    # Zprev prime to p and Z = 0 mod p proves a point of order 2^k modulo every prime q|p
    # we chose 2^k > q + 1 + 2*sqrt(q), so this is possible only when p is prime
    return gcd(Zprev, p) == 1 and Z % p == 0

a = pp_verify(127,103,51)
print(a)