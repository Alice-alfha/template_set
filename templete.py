# エラトステネスの篩（戻り値: 素数の昇順list）
# 10^7で 80msくらい
def eratosthenes(maxN: int):
    sieve = bytearray([0, 1]) * ((maxN+1)//2)
    sieve.extend([0] * ((maxN+1) % 2))

    sieve[1] = 0
    sieve[2] = 1

    for i in range(3, int(pow(maxN+1, 0.5))+1, 2):
        if sieve[i]:
            sieve[i*i::i*2] = bytearray(len(range(i*i, maxN+1, i*2)))

    return [x for x in range(maxN+1) if sieve[x]]


# 線形篩 (戻り値: 最小素因数のlist)
# 10^7で 70msくらい
def linear_sieve(maxN: int):
    spf = [x for x in range(maxN+1)]
    p = []

    for i in range(2, (maxN+1)//2):
        if spf[i] == i:
            p.append(i)

        for j in p:
            if i*j > maxN:
                break
            spf[i*j] = j

            if j == spf[i]:
                break

    return spf


# BIT
# add: bit[p]にxを加算する
# sum: l <= x < rの和
class Fenwick_Tree:
    def __init__(self, n):
        self._n = n
        self.data = [0] * n

    def add(self, p, x):
        assert 0 <= p < self._n
        p += 1
        while p <= self._n:
            self.data[p - 1] += x
            p += p & -p

    def sum(self, l, r):
        assert 0 <= l <= r <= self._n
        return self._sum(r) - self._sum(l)

    def _sum(self, r):
        s = 0
        while r > 0:
            s += self.data[r - 1]
            r -= r & -r
        return s
