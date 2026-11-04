# AI製セグメントツリー

# Segment Tree
# 0-indexed
# max_right: f(op(l, l+1, ..., r))=Trueとなる最大のrを返す
class SegTree:
    def __init__(self, op, e, n_or_a):
        if isinstance(n_or_a, int):
            a = []
            n = n_or_a
        else:
            a = list(n_or_a)
            n = len(a)
        size = 1
        while size < n:
            size *= 2
        d = [e] * (2 * size)
        if a:
            d[size:size + n] = a
            for i in range(size - 1, 0, -1):
                d[i] = op(d[2 * i], d[2 * i + 1])
        self.op = op
        self.e = e
        self.n = n
        self.size = size
        self.d = d

    def set(self, p, x):
        d = self.d
        op = self.op
        p += self.size
        d[p] = x
        while p > 1:
            p >>= 1
            d[p] = op(d[2 * p], d[2 * p + 1])

    def get(self, p):
        return self.d[p + self.size]

    def prod(self, l, r):
        d = self.d
        op = self.op
        sml = smr = self.e
        l += self.size
        r += self.size
        while l < r:
            if l & 1:
                sml = op(sml, d[l])
                l += 1
            if r & 1:
                r -= 1
                smr = op(d[r], smr)
            l >>= 1
            r >>= 1
        return op(sml, smr)

    def all_prod(self):
        return self.d[1]

    def max_right(self, l, f):
        n = self.n
        if l >= n:
            return n
        d = self.d
        op = self.op
        size = self.size
        l += size
        sm = self.e
        while True:
            while l % 2 == 0:
                l >>= 1
            if not f(op(sm, d[l])):
                while l < size:
                    l *= 2
                    if f(op(sm, d[l])):
                        sm = op(sm, d[l])
                        l += 1
                return l - size
            sm = op(sm, d[l])
            l += 1
            if (l & -l) == l:
                return n

    def min_left(self, r, f):
        if r <= 0:
            return 0
        d = self.d
        op = self.op
        size = self.size
        r += size
        sm = self.e
        while True:
            r -= 1
            while r > 1 and r % 2:
                r >>= 1
            if not f(op(d[r], sm)):
                while r < size:
                    r = 2 * r + 1
                    if f(op(d[r], sm)):
                        sm = op(d[r], sm)
                        r -= 1
                return r + 1 - size
            sm = op(d[r], sm)
            if (r & -r) == r:
                return 0


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
