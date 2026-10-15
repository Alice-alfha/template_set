# Union-find (毎回書くのが面倒)
# find: rootノード探索, union: 2ノードを接続, connected: 2ノードの接続判定
# rank: rootノードにするノードの優先度 (高いほどより優先的にrootになる)
class Unionfind:
    def __init__(self, size):
        self.root = [x for x in range(size + 1)]
        self.rank = [1] * (size + 1)

    def find(self, x):
        if x != self.root[x]:
            self.root[x] = self.find(self.root[x])

        return self.root[x]

    def union(self, x, y):
        rootx = self.find(x)
        rooty = self.find(y)

        if rootx != rooty:
            if self.rank[rootx] > self.rank[rooty]:
                self.root[rooty] = rootx
            elif self.rank[rootx] < self.rank[rooty]:
                self.root[rootx] = rooty
            else:
                self.root[rooty] = rootx
                self.rank[rootx] += 1

    def connected(self, x, y):
        return self.find(x) == self.find(y)


# Segment Tree
# 0-indexed
# class_function: set, get, prod, all_prod
class SegTree:
    from typing import Union

    def __init__(self, values: int | list, op, ide):

        if isinstance(values, int):
            rank = 1 << (values-1).bit_length()
        elif isinstance(values, list):
            rank = 1 << (len(values)-1).bit_length()
        else:
            raise NotImplementedError()

        self._op = op
        self._ide = ide
        self._tree_size = rank
        tree = [ide for _ in range(rank * 2)]

        if isinstance(values, list):
            for i in range(len(values)):
                tree[rank+i] = values[i]

            for i in range(rank-1, 0, -1):
                tree[i] = op(tree[i << 1], tree[(i << 1) + 1])

        self.tree = tree

    def set(self, index: int, element) -> None:
        tree = self.tree
        index += self._tree_size
        tree[index] = element
        op = self._op

        while index > 1:
            # update関数実行時間の9割はこの1行
            index >>= 1
            tree[index] = op(tree[index << 1], tree[index << 1 | 1])

    def get(self, index: int) -> any:
        return self.tree[self._tree_size + index]

    def prod(self, l_: int, r_: int) -> any:
        l_res = self._ide
        r_res = self._ide
        tree = self.tree
        op_ = self._op
        l_ += self._tree_size
        r_ += self._tree_size - 1

        if l_ == r_:
            return tree[l_]

        while l_ <= r_:
            if l_ & 1:
                l_res = op_(l_res, tree[l_])
                l_ += 1
            if not r_ & 1:
                r_res = op_(tree[r_], r_res)
                r_ -= 1

            l_ >>= 1
            r_ >>= 1

        return op_(l_res, r_res)

    def all_prod(self) -> any:
        return self.tree[1]


# Lazy Segment Tree（完成していない）
# class_function: set, get, prod
class LazySegTree:
    def __init__(self, values: list, segfunc, ide):
        n = len(values)
        r = 1 << (n-1).bit_length()

        self.tree_size = r
        self.tree = [ide] * 2 * r
        self.lazy = [ide] * 2 * r

        for i in range(n):
            self.tree[r + i] = values[i]

        self.ide = ide
        self.seg_func = segfunc
        self.update_all()

    def set(self, value, l, r):
        lazytree = self.lazy
        l += self.tree_size
        r += self.tree_size-1

        while l <= r:
            if l & 1:
                lazytree[l] = self.seg_func(lazytree[l], value)
                l += 1
            if not r & 1:
                lazytree[r] = self.seg_func(lazytree[r], value)
                r -= 1

            l >>= 1
            r >>= 1

    def get(self, x):
        return self.tree[self.tree_size + x]

    def prod(self, l, r):
        if l >= r:
            return None

        size = self.tree_size
        tree = self.tree
        lazy = self.lazy
        ide = self.ide
        ind: int = size + l

        self.__lazy_update(l)
        tree[ind] = self.seg_func(tree[ind], lazy[ind])
        lazy[ind] = ide

        for i in range(l+1, r):
            ind += 1
            self.__lazy_update(i, (ind & -ind).bit_length() - 1)
            tree[ind] = self.seg_func(tree[ind], lazy[ind])
            lazy[ind] = ide

        self.update_all()

        res = 0
        l += size
        r += size - 1

        while l <= r:
            if l & 1:
                res = self.seg_func(tree[l], res)
                l += 1
            if not r & 1:
                res = self.seg_func(tree[r], res)
                r -= 1
            l >>= 1
            r >>= 1
        return res

    def update_all(self):
        nsize = self.tree_size >> 1
        tree = self.tree

        while nsize >= 1:
            compare = nsize << 1
            for i in range(nsize):
                tree[nsize+i] = self.seg_func(tree[compare], tree[compare+1])
                compare += 2

            nsize >>= 1

    def __lazy_update(self, index, loop=None):
        lazy = self.lazy
        index += self.tree_size
        nindex = 0
        l = 0
        r = 0

        if loop is None:
            loop = self.tree_size.bit_length()-1

        for i in range(loop, 0, -1):
            nindex = index >> i
            if lazy[nindex] == 0:
                continue

            l, r = nindex << 1, 1 + (nindex << 1)
            lazy[l] = self.seg_func(lazy[nindex], lazy[l])
            lazy[r] = self.seg_func(lazy[nindex], lazy[r])
            lazy[nindex] = 0


# エラトステネスの篩
# 10 ** 7で 600msくらい
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
