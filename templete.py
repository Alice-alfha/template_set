# エラトステネスの篩
# 10^7で 600msくらい
def eratosthenes(maxN: int):
    prime = [2]
    m = set()

    for i in range(3, maxN, 2):
        if i in m:
            continue

        prime.append(i)

        t = pow(i, 2)
        while maxN > t:
            m.add(t)
            t += i * 2

    return prime


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
