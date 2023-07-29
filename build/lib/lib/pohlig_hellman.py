def extended_euclidean(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = extended_euclidean(b % a, a)
        return (g, x - (b // a) * y, y)
 
def modinv(a, m):
    g, x, y = extended_euclidean(a, m)
    return x % m
 
def chinese_remainder(m, x):
    '''
    finds value X so X = x_i (mod m_i) for all i
    '''
    while len(x) > 1:
        tmp1 = modinv(m[-2],m[-1]) * x[-1] * m[-2] + modinv(m[-1],m[-2]) * x[-2] * m[-1]
        tmp2 = m[-1] * m[-2]
        x.pop()
        x.pop()
        x.append(x[tmp1 % tmp2])
        m.pop()
        m.pop()
        m.append(m[tmp2])
    return x[0]


def baby_giant_step(g, h, n):
    '''
    Finds value x satisfying g^x = h, where |<g>| = n
    O(sqrt{n})
    '''
    m = int(1 + n**0.5)
    dic = {}
    a = 1
    for j in range(m):
        dic[a] = j
        a = a * g % n

    c = pow(g, -m, n)
    l = h
    for i in range(m):
        if l in dic:
            return i * m + dic[l]
        l = l * c % n

    return -1


def pohlig_hellman_prime(g, h, p, e):
    n = p**e
    x = 0
    l = pow(g, n // p, n)
    for k in range(e):
        h_k = pow(pow(g,-x,n)*h,p**(e-1-k),n)
        d_k = baby_giant_step(l,h_k,p)
        if d_k == -1:
            return -1
        x = (x + pow(p,k,n) * d_k) % n
    return x


def pohlig_hellman_general(g, h, n, factors):
    '''
    Given a cyclic group G of order n and generator <g>, an element h \in G and a primefactorization n = \prod{p_i^e_i}
    In other word find x such that g^x = h (mod n)
    factors should be list of tuples (prime, exponent)
    O(sum{e_i(log n + sqrt{p_i})})
    '''
    r = len(factors)
    X = []
    P = []
    ans = 0
    for i in range(r):
        G = (pow(g,n//(factors[i][0]**factors[i][1]),n))
        H = (pow(h,n//(factors[i][0]**factors[i][1]),n))
        X.append(pohlig_hellman_prime(G, H, factors[i][0], factors[i][1]))
        if X[-1] == -1:
            return -1
        P.append(factors[i][0]**factors[i][1])
    return chinese_remainder(P,X)


