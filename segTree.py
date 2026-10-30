# Generate AI

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


# Lazy Segment Tree
# 0-indexed
# op: segTree_function, e: segTree_e, mapping: f(lazy, seg)
# composition: lazy_function, id_: lazy_e, v: segTree_init
class LazySegTree:
    def __init__(self, op, e, mapping, composition, id_, v):
        self.op = op
        self.e = e
        self.mapping = mapping
        self.composition = composition
        self.id = id_

        if isinstance(v, int):
            v = [e] * v

        n = len(v)
        self.n = n
        self.log = max(n - 1, 0).bit_length()
        self.size = 1 << self.log
        self.d = [e] * (2 * self.size)
        self.lz = [id_] * self.size

        d = self.d
        size = self.size
        for i in range(n):
            d[size + i] = v[i]

        op_ = self.op
        for i in range(size - 1, 0, -1):
            d[i] = op_(d[2 * i], d[2 * i + 1])

    def set(self, p, x):
        p += self.size
        log = self.log
        d = self.d
        lz = self.lz
        mapping = self.mapping
        composition = self.composition
        id_ = self.id
        op_ = self.op
        size = self.size

        for i in range(log, 0, -1):
            k = p >> i
            f = lz[k]
            lk, rk = 2 * k, 2 * k + 1
            d[lk] = mapping(f, d[lk])
            if lk < size:
                lz[lk] = composition(f, lz[lk])
            d[rk] = mapping(f, d[rk])
            if rk < size:
                lz[rk] = composition(f, lz[rk])
            lz[k] = id_

        d[p] = x

        for i in range(1, log + 1):
            k = p >> i
            d[k] = op_(d[2 * k], d[2 * k + 1])

    def get(self, p):
        p += self.size
        log = self.log
        d = self.d
        lz = self.lz
        mapping = self.mapping
        composition = self.composition
        id_ = self.id
        size = self.size

        for i in range(log, 0, -1):
            k = p >> i
            f = lz[k]
            lk, rk = 2 * k, 2 * k + 1
            d[lk] = mapping(f, d[lk])
            if lk < size:
                lz[lk] = composition(f, lz[lk])
            d[rk] = mapping(f, d[rk])
            if rk < size:
                lz[rk] = composition(f, lz[rk])
            lz[k] = id_
        return d[p]

    def prod(self, l, r):
        if l == r:
            return self.e

        l += self.size
        r += self.size
        log = self.log
        d = self.d
        lz = self.lz
        op_ = self.op
        mapping = self.mapping
        composition = self.composition
        id_ = self.id
        size = self.size

        for i in range(log, 0, -1):
            if ((l >> i) << i) != l:
                k = l >> i
                f = lz[k]
                lk, rk = 2 * k, 2 * k + 1
                d[lk] = mapping(f, d[lk])
                if lk < size:
                    lz[lk] = composition(f, lz[lk])
                d[rk] = mapping(f, d[rk])
                if rk < size:
                    lz[rk] = composition(f, lz[rk])
                lz[k] = id_
            if ((r >> i) << i) != r:
                k = (r - 1) >> i
                f = lz[k]
                lk, rk = 2 * k, 2 * k + 1
                d[lk] = mapping(f, d[lk])
                if lk < size:
                    lz[lk] = composition(f, lz[lk])
                d[rk] = mapping(f, d[rk])
                if rk < size:
                    lz[rk] = composition(f, lz[rk])
                lz[k] = id_

        sml = smr = self.e
        while l < r:
            if l & 1:
                sml = op_(sml, d[l])
                l += 1
            if r & 1:
                r -= 1
                smr = op_(d[r], smr)
            l >>= 1
            r >>= 1
        return op_(sml, smr)

    def all_prod(self):
        return self.d[1]

    def apply(self, p, f):
        p += self.size
        log = self.log
        d = self.d
        lz = self.lz
        mapping = self.mapping
        composition = self.composition
        id_ = self.id
        op_ = self.op
        size = self.size

        for i in range(log, 0, -1):
            k = p >> i
            fx = lz[k]
            lk, rk = 2 * k, 2 * k + 1
            d[lk] = mapping(fx, d[lk])
            if lk < size:
                lz[lk] = composition(fx, lz[lk])
            d[rk] = mapping(fx, d[rk])
            if rk < size:
                lz[rk] = composition(fx, lz[rk])
            lz[k] = id_

        d[p] = mapping(f, d[p])

        for i in range(1, log + 1):
            k = p >> i
            d[k] = op_(d[2 * k], d[2 * k + 1])

    def range_apply(self, l, r, f):
        if l == r:
            return
        l += self.size
        r += self.size
        log = self.log
        d = self.d
        lz = self.lz
        mapping = self.mapping
        composition = self.composition
        id_ = self.id
        op_ = self.op
        size = self.size

        for i in range(log, 0, -1):
            if ((l >> i) << i) != l:
                k = l >> i
                fx = lz[k]
                lk, rk = 2 * k, 2 * k + 1
                d[lk] = mapping(fx, d[lk])
                if lk < size:
                    lz[lk] = composition(fx, lz[lk])
                d[rk] = mapping(fx, d[rk])
                if rk < size:
                    lz[rk] = composition(fx, lz[rk])
                lz[k] = id_
            if ((r >> i) << i) != r:
                k = (r - 1) >> i
                fx = lz[k]
                lk, rk = 2 * k, 2 * k + 1
                d[lk] = mapping(fx, d[lk])
                if lk < size:
                    lz[lk] = composition(fx, lz[lk])
                d[rk] = mapping(fx, d[rk])
                if rk < size:
                    lz[rk] = composition(fx, lz[rk])
                lz[k] = id_

        l2, r2 = l, r
        while l < r:
            if l & 1:
                d[l] = mapping(f, d[l])
                if l < size:
                    lz[l] = composition(f, lz[l])
                l += 1
            if r & 1:
                r -= 1
                d[r] = mapping(f, d[r])
                if r < size:
                    lz[r] = composition(f, lz[r])
            l >>= 1
            r >>= 1

        l, r = l2, r2
        for i in range(1, log + 1):
            if ((l >> i) << i) != l:
                k = l >> i
                d[k] = op_(d[2 * k], d[2 * k + 1])
            if ((r >> i) << i) != r:
                k = (r - 1) >> i
                d[k] = op_(d[2 * k], d[2 * k + 1])
