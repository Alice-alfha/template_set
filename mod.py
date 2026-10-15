# https://strangerxxx.hateblo.jp/entry/20230419/1681873929の記事参考（というかほぼそのまま）

# 有理数のmodをとるやつ
def invmod(n, m, mod=998244353):
    r = pow(m, -1, mod)
    return (n * r) % mod
