class Crypto:

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

    def chinese_remainder(m, x):
        '''
        finds value X so X = x_i (mod m_i) for all i
        '''
        assert len(m) == len(x)
        n = len(m)
        prod = 1
        for nums in m:
            prod *= nums

        ret = 0
        for i in range(n):
            pp = prod // m[i]
            ret = ret + x[i] * Crypto.modinv(pp, m[i]) * pp

        return ret % prod

    def determinant(M):
        '''
        given Matrix M returns the determinant of the Matrix, if the matrix isn't a N x N matrix it will return 0
        '''
        n = len(M)

        if n != len(M[0]):
            return 0

        AM = M[:]
        for fd in range(n):
            for i in range(fd+1, n):
                if AM[fd][fd] == 0:
                    AM[fd][fd] == 1.0e-18
                c = AM[i][fd] / AM[fd][fd]
                AM[i] = Crypto.vecsub(AM[i], Crypto.scalar(AM[fd],c))

        ret = 1.0
        for i in range(n):
            ret *= AM[i][i]

        return ret


    def dot_product(a, b):
        '''
        Takes two vectors of same size and produces the dot product
        '''
        assert len(a) == len(b)
        n = len(a)
        su = 0
        for i in range(n):
            su += a[i] * b[i]
        return su


    def extended_euclidean(a, b):
        '''
        given numbers a, b returns (g, x, y), where g is the gcd, and x, y satisfies ax + by = gcd(a, b)
        Special if b is prime y will be modinv of a
        Same if a is prime x will be modinv of b
        '''
        if a == 0:
            return (b, 0, 1)
        else:
            g, y, x = Crypto.extended_euclidean(b % a, a)
            return (g, x - (b // a) * y, y)

    def gauss_lattice(v1, v2):
        '''
        Given two vectors v1, v2. Finds an optimal basis for a two-dimensional lattice. 
        The algorithm finds u1, u2, so u1 solves SVP
        '''
        while True:
            if Crypto.dot_product(v2, v2) < Crypto.dot_product(v1, v1):
                v1,v2=v2,v1
            m = Crypto.dot_product(v1,v2) // Crypto.dot_product(v1,v1)
            if m == 0:
                return v1, v2
            v2 = Crypto.vecsub(v2, Crypto.scalar(v1,m))

    def gauss_elimination_mod(v, p):
        '''
        Given a matrix V with dimensions n x n, returns the reduced echelon form of the matrix
        '''
        n = len(v)
        m = len(v[0])
        v_copy = v[:]

        for i in range(n):
            for j in range(m):
                v_copy[i][j] %= p
        
        h = 0
        k = 0

        while h < n and k < m:
            if v_copy[h][k] == 0:
                for i in range(h,n):
                    if v_copy[i][k] != 0:
                        v_copy[i], v_copy[h] = v_copy[h], v_copy[i]
                        break
            if v_copy[h][k] == 0:
                k += 1
            else:
                inv = Crypto.modinv(v_copy[h][k], p)
                for i in range(k,m):
                    v_copy[h][i] = (v_copy[h][i] * inv) % p
                for i in range(n):
                    if i == h:
                        continue
                    f = v_copy[i][k]
                    for j in range(k,m):
                        v_copy[i][j] -= f * v_copy[h][j]
                        v_copy[i][j] %= p
                h += 1
                k += 1

        return v_copy
     
    def gram_schmidt(v):
        '''
        Given a basis v_1, v_2, ..., v_n for a vector space V computes an orthogonal basis u_1, u_2, ..., u_n for the same vector space V
        If v_1, v_2, ..., v_n is not linearly indepent this will return nonsense or throw an error
        vectors should be of the form [v_i_0, v_i_1, ..., v_i_n], using tuples resolves in errors
        '''
        n = len(v)
        u = [v[0]]
        for i in range(1, n):
            mu = []
            for j in range(i):
                mu.append(Crypto.dot_product(v[i], u[j]) / Crypto.dot_product(u[j],u[j]))

            u.append(v[i])

            for j in range(i):
                u[i] = Crypto.vecsub(u[i], Crypto.scalar(u[j], mu[j]))

        return u

    def LLL(B, delta):
        '''
        https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm
        Lenstra-Lenstra-Lovász lattice basis reduction algorithm
        Given a lattice basis b1, b2, ..., bn in Z^m, and a delta: 1/4 < delta < 1, most commonly delta = 3/4, greater delta lead to stronger reductions
        of the basis. For delta = 1, not guaranteed polynomial time
        O(n^5 * m * log^3(max_i(||b_i||)))
        '''
        n = len(B)
        k = 1
        _B = Crypto.gram_schmidt(B)

        def mu(i, j):
            return Crypto.dot_product(B[i], _B[j]) / Crypto.dot_product(_B[j], _B[j])

        while k < n:
            for j in range(k - 1, -1, -1):
                mu_kj = mu(k,j)
                if abs(mu_kj) > 0.5:
                    B[k] = Crypto.vecsub(B[k], Crypto.scalar(B[j], round(mu_kj)))
                    _B = Crypto.gram_schmidt(B)

            if Crypto.dot_product(_B[k], _B[k]) >= (delta - mu(k,k-1)**2) * Crypto.dot_product(_B[k-1], _B[k-1]):
                k += 1
            else:
                B[k], B[k-1] = B[k-1], B[k]
                _B = Crypto.gram_schmidt(B)
                k = max(k-1, 1)

        return B

    def modinv(a, m):
        '''
        Calculates inverse of a mod m using extended_euclidean
        '''
        g, x, y = Crypto.extended_euclidean(a, m)
        return x % m
     

    def pohlig_hellman_prime(g, h, p, e):
        n = p**e
        x = 0
        l = pow(g, n // p, n)
        for k in range(e):
            h_k = pow(pow(g,-x,n)*h,p**(e-1-k),n)
            d_k = Crypto.baby_giant_step(l,h_k,p)
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
            X.append(Crypto.pohlig_hellman_prime(G, H, factors[i][0], factors[i][1]))
            if X[-1] == -1:
                return -1
            P.append(factors[i][0]**factors[i][1])
        return Crypto.chinese_remainder(P,X)

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
    
    def scalar(v, c):
        '''
        multiply a vector v with a scalar c
        '''
        v_copy = v[:]
        for i in range(len(v)):
            v_copy[i] *= c
        return v_copy

    def vecsub(v, u):
        '''
        substracts a vector u from vector v
        '''
        assert len(v) == len(u)
        n = len(v)
        v_copy = v[:]
        for i in range(n):
            v_copy[i] -= u[i]
        return v_copy



