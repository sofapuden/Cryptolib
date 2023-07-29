def tonelli_shank(n, p):
    '''
    Find a squareroot x of n modulo p satisfying x^2 = n (mod p). Assumes p is an odd prime, or p is of the form q^k for prime q
    returns -1 if no such x exists
    https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    '''

    # n = Q * 2^S
    Q = p - 1
    S = 0
    while Q & 1 == 0:
        Q >>= 1
        S += 1

    # Find a quadratic non-residue
    z = 2
    while pow(z, (p-1) // 2, p) != p - 1:
        z += 1

    M = S
    c = pow(z, Q, p)
    t = pow(n, Q, p)
    R = pow(n, (Q + 1) // 2, p)

    # Do the magic
    while t:
        if t == 1:
            return R
        i = 0
        for i in range(1, M + 1):
            if pow(t, pow(2, i, p-1), p) == 1:
                break

        if i == M:
            return -1

        b = pow(c, pow(2, M-i-1, p-1), p)
        M = i
        c = b * b % p
        t = t * c % p
        R = R * b % p

    return -1

