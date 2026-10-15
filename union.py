# https://github.com/Alice-alfha/template_set/blob/main/union.py
# AI製()
class UnionFind:
    def __init__(self, n: int):
        self.parent = [-1] * (n+1)
        self.groups = n

    def find(self, x: int) -> int:
        parent = self.parent
        while parent[x] >= 0:
            if parent[parent[x]] >= 0:
                parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(self, a: int, b: int) -> bool:
        parent = self.parent
        a = self.find(a)
        b = self.find(b)

        if a == b:
            return False

        # 要素数が多い方を親にする
        if parent[a] > parent[b]:
            a, b = b, a

        parent[a] += parent[b]
        parent[b] = a
        self.groups -= 1
        return True

    def connected(self, a: int, b: int) -> bool:
        return self.find(a) == self.find(b)

    def size(self, x: int) -> int:
        return -self.parent[self.find(x)]

    def count(self) -> int:
        return self.groups
